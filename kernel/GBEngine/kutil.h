#ifndef KUTIL_H
#define KUTIL_H
/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/*
* ABSTRACT: kernel: utils for kStd
*/


#include <string.h>
#include <stdlib.h>

#include "writable_set.h"
#include <vector>
#include <unordered_map>
#include <utility>

#include "omalloc/omalloc.h"
#ifdef HAVE_OMALLOC
#include "omalloc/omallocClass.h"
#endif

#include "misc/mylimits.h"

#include "kernel/polys.h"
#include "polys/operations/pShallowCopyDelete.h"

#include "kernel/structs.h"
#include "kernel/GBEngine/kstd1.h"   /* for s_poly_proc_t */
#include "coeffs/bigintmat.h"   /* for s_poly_proc_t */

// define if tailrings should be used
#define HAVE_TAIL_RING

#define setmax 128
#define setmaxT ((int)((4096-12)/sizeof(TObject)))
#define setmaxTinc ((int)((4096)/sizeof(TObject)))

#define RED_CANONICALIZE 200
#define REDNF_CANONICALIZE 60
#define REDTAIL_CANONICALIZE 100

// if you want std computations as in Singular version < 2:
// This disables RedThrough, tailReductions against T (bba),
// sets posInT = posInT15 (bba, strat->honey), and enables redFirst with LDeg
// NOTE: can be achieved with option(oldStd)

#undef NO_KINLINE
#if !defined(KDEBUG) && !defined(NO_INLINE)
#define KINLINE inline
#else
#define KINLINE
#define NO_KINLINE 1
#endif

typedef int* intset;
typedef int64  wlen_type;
typedef wlen_type* wlen_set;

// Block-allocated array: elements are stored in fixed-size blocks.
// Existing blocks never move when new blocks are appended, so pointers
// to elements (&arr[i]) remain stable across growth.
template<typename Elem, int BLOCK_SHIFT = 10>
class BlockArray {
  static const int BLOCK_SIZE = 1 << BLOCK_SHIFT;
  static const int BLOCK_MASK = BLOCK_SIZE - 1;
  Elem **blocks;        // directory of block pointers
  int num_blocks;       // current number of allocated blocks
  int dir_capacity;     // allocated directory slots
protected:
  int count;            // number of elements in use
public:
  BlockArray() : blocks(NULL), num_blocks(0), dir_capacity(0), count(0) {}

  Elem& operator[](int i) {
    return blocks[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }
  const Elem& operator[](int i) const {
    return blocks[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }

  // Return pointer to element i (stable across growth)
  Elem* addr(int i) {
    return &blocks[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }

  // Number of elements in use
  int size() const { return count; }
  bool empty() const { return count == 0; }

  // Set the count directly (for migration from external tl/sl counters)
  void setsize(int n) { count = n; }

  // Ensure at least n elements are allocated (indices 0..n-1)
  void ensure_capacity(int n) {
    int needed_blocks = (n + BLOCK_SIZE - 1) >> BLOCK_SHIFT;
    if (needed_blocks <= num_blocks) return;
    // Grow directory if needed
    if (needed_blocks > dir_capacity) {
      int new_cap = dir_capacity == 0 ? 4 : dir_capacity;
      while (new_cap < needed_blocks) new_cap *= 2;
      Elem **new_dir = (Elem **)calloc(new_cap, sizeof(Elem *));
      if (blocks != NULL) {
        memcpy(new_dir, blocks, num_blocks * sizeof(Elem *));
        free(blocks);
      }
      blocks = new_dir;
      dir_capacity = new_cap;
    }
    // Allocate new blocks (zero-initialized)
    for (int b = num_blocks; b < needed_blocks; b++) {
      blocks[b] = (Elem *)calloc(BLOCK_SIZE, sizeof(Elem));
    }
    num_blocks = needed_blocks;
  }

  // Return current capacity (total elements allocated)
  int capacity() const { return num_blocks * BLOCK_SIZE; }

  // Insert val at position pos, shifting existing elements up.
  // Automatically grows capacity and increments count.
  void insert(int pos, const Elem& val) {
    ensure_capacity(count + 1);
    for (int i = count; i > pos; i--)
      (*this)[i] = (*this)[i-1];
    (*this)[pos] = val;
    count++;
  }

  // Append val at the end. Grows capacity and increments count.
  void push_back(const Elem& val) {
    ensure_capacity(count + 1);
    (*this)[count] = val;
    count++;
  }

  // Erase element at position pos, shifting elements down. Decrements count.
  void erase(int pos) {
    for (int i = pos; i < count - 1; i++)
      (*this)[i] = (*this)[i+1];
    count--;
  }

  // Free all blocks and the directory, reset count
  void free_all() {
    for (int b = 0; b < num_blocks; b++) {
      free(blocks[b]);
    }
    if (blocks != NULL) free(blocks);
    blocks = NULL;
    num_blocks = 0;
    dir_capacity = 0;
    count = 0;
  }
};

typedef class sTObject TObject;
typedef class sLObject LObject;

typedef struct denominator_list_s denominator_list_s;
typedef denominator_list_s *denominator_list;

struct denominator_list_s{number n; denominator_list next;};
EXTERN_VAR denominator_list DENOMINATOR_LIST;

class sTObject
{
public:
  unsigned long sevSig;
  poly sig;   // the signature of the element
  poly p;       // Lm(p) \in currRing Tail(p) \in tailRing
  poly t_p;     // t_p \in tailRing: as monomials Lm(t_p) == Lm(p)
  poly max_exp;     // p_GetMaxExpP(pNext(p))
  ring tailRing;
  long FDeg;    // pFDeg(p)
  int ecart,
    length,     // as of pLDeg
    pLength,    // either == 0, or == pLength(p)
    i_r;        // index of TObject in R set, or -1 if not in T

#ifdef HAVE_SHIFTBBA
  int shift;
#endif

  /*BOOLEAN*/ char is_normalized; // true, if pNorm was called on p, false otherwise
  // used in incremental sba() with F5C:
  // we know some of the redundant elements in
  // strat->T beforehand, so we can just discard
  // them and do not need to consider them in the
  // interreduction process
  /*BOOLEAN*/ char is_redundant;
  // used in sba's sig-safe reduction:
  // sometimes we already know that a reducer
  // is sig-safe, so no need for a real
  // sig-safeness check
  /*BOOLEAN*/ char is_sigsafe;


#ifdef HAVE_PLURAL
  /*BOOLEAN*/ char is_special; // true, it is a new special S-poly (e.g. for SCA)
#endif

  // initialization
  KINLINE void Init(ring r = currRing);
  KINLINE sTObject(ring tailRing = currRing);
  KINLINE sTObject(poly p, ring tailRing = currRing);
  KINLINE sTObject(poly p, ring c_r, ring tailRing);
  KINLINE sTObject(sTObject* T, int copy);

  KINLINE void Set(ring r=currRing);
  KINLINE void Set(poly p_in, ring r=currRing);
  KINLINE void Set(poly p_in, ring c_r, ring t_r);

  // Frees the polys of T
  KINLINE void Delete();
  // Sets polys to NULL
  KINLINE void Clear();
  // makes a copy of the poly of T
  KINLINE void Copy();

  // ring-dependent Lm access: these might result in allocation of monomials
  KINLINE poly GetLmCurrRing();
  KINLINE poly GetLmTailRing();
  KINLINE poly GetLm(ring r);
  // this returns Lm and ring r (preferably from tailRing), but does not
  // allocate a new poly
  KINLINE void GetLm(poly &p, ring &r) const;

#ifdef OLIVER_PRIVAT_LT
  // routines for calc. with rings
  KINLINE poly GetLtCurrRing();
  KINLINE poly GetLtTailRing();
  KINLINE poly GetLt(ring r);
  KINLINE void GetLt(poly &p, ring &r) const;
#endif

  KINLINE BOOLEAN IsNull() const;

  KINLINE int GetpLength();

  // makes sure that T.p exists
  KINLINE void SetLmCurrRing();

  // Iterations
  // simply get the next monomial
  KINLINE poly Next();
  KINLINE void LmDeleteAndIter();

  // deg stuff
  // compute pTotalDegree
  KINLINE long pTotalDeg() const;
  // computes pFDeg
  KINLINE long pFDeg() const;
  // computes and sets FDeg
  KINLINE long SetpFDeg();
  // gets stored FDeg
  KINLINE long GetpFDeg() const;

  // computes pLDeg
  KINLINE long pLDeg();
  // sets length, FDeg, returns LDeg
  KINLINE long SetDegStuffReturnLDeg();

  // arithmetic
  KINLINE void Mult_nn(number n);
  KINLINE void ShallowCopyDelete(ring new_tailRing, omBin new_tailBin,
                                 pShallowCopyDeleteProc p_shallow_copy_delete,
                                 BOOLEAN set_max = TRUE);
  // manipulations
  KINLINE void pNorm();
  KINLINE void pCleardenom();
  KINLINE void pContent();

#ifdef KDEBUG
  void wrp();
#endif
};

EXTERN_VAR int strat_nr;

class sLObject : public sTObject
{

public:
  unsigned long sev;
  unsigned long sev_lcm;  // short exp vector of lcm, for fast divisibility pre-check
  poly  p1,p2; /*- the pair p comes from,
                 lm(pi) in currRing, tail(pi) in tailring -*/

  poly  lcm;   /*- the lcm of p1,p2 -*/
  kBucket_pt bucket;
  int   i_r1, i_r2;
  unsigned seq;       // the sequence number of the LSet when this LObject was inserted
                      // used to determine LSet ordering for equal LObjects
  size_t flat_index;      // index in writable_set's flat array for unordered iteration
  unsigned checked; // this is the index of S up to which
                      // the corresponding LObject was already checked in
                      // critical pair creation => when entering the
                      // reduction process it is enough to start a second
                      // rewritten criterion check from checked+1 onwards
  BOOLEAN prod_crit;
                      // NOTE: If prod_crit = TRUE then the corresponding pair is
                      // detected by Buchberger's Product Criterion and can be
                      // deleted

  // initialization
  KINLINE void Init(ring tailRing = currRing);
  KINLINE sLObject(ring tailRing = currRing);
  KINLINE sLObject(poly p, ring tailRing = currRing);
  KINLINE sLObject(poly p, ring c_r, ring tailRing);

  // Frees the polys of L
  KINLINE void Delete();
  KINLINE void Clear();

  // Iterations
  KINLINE void LmDeleteAndIter();
  KINLINE poly LmExtractAndIter();

  // spoly related things
  // preparation for reduction if not spoly
  KINLINE void PrepareRed(BOOLEAN use_bucket);
  KINLINE void SetLmTail(poly lm, poly new_p, int length,
                         int use_bucket, ring r);
  KINLINE void Tail_Minus_mm_Mult_qq(poly m, poly qq, int lq, poly spNoether);
  KINLINE void Tail_Mult_nn(number n);
  // deletes bucket, makes sure that p and t_p exists
  KINLINE poly GetP(omBin lmBin = (omBin)NULL);
  // similar, except that only t_p exists
  KINLINE poly GetTP();

  // does not delete bucket, just canonicalizes it
  // returned poly is such that Lm(p) \in currRing, Tail(p) \in tailRing
  KINLINE void CanonicalizeP();

  // makes a copy of the poly of L
  KINLINE void Copy();

  KINLINE int GetpLength();
  KINLINE long pLDeg(BOOLEAN use_last);
  KINLINE long pLDeg();
  KINLINE int SetLength(BOOLEAN lengt_pLength = FALSE);
  KINLINE long SetDegStuffReturnLDeg();
  KINLINE long SetDegStuffReturnLDeg(BOOLEAN use_last);

  // returns minimal component of p
  KINLINE long MinComp();
  // returns component of p
  KINLINE long Comp();

  KINLINE void ShallowCopyDelete(ring new_tailRing,
                                 pShallowCopyDeleteProc p_shallow_copy_delete);

  // sets sev
  KINLINE void SetShortExpVector();

  // enable assignment from TObject
  KINLINE sLObject& operator=(const sTObject&);

  // get T's corresponding to p1, p2: they might return NULL
  KINLINE TObject* T_1(const skStrategy* strat);
  KINLINE TObject* T_2(const skStrategy* strat);
  KINLINE void     T_1_2(const skStrategy* strat,
                         TObject* &T_1, TObject* &T_2);

  // simplify coefficients
  KINLINE void Normalize();
  KINLINE void HeadNormalize();
};

EXTERN_VAR int HCord;

/** @class LSet
 *
 * "L" is the sorted set of critical pairs, and we wish to regularly
 * pop the top item from the queue.  However, we also wish to iterate
 * over the entire set in order, which precludes organizing it as a
 * heap.  Also, our legacy C code modifies the objects once they're in
 * the queue, which precludes using std::map or std::multimap, since
 * their objects are immutable.  We use a custom class writable_set,
 * which is derived from std::multimap, but stores pointers in the
 * tree and therefore allows the objects themselves to be modified.
 */

class CompareLObject {
public:
  skStrategy * strat = NULL;
  int (*compareL) (const LObject &lhs, const LObject &rhs, const kStrategy strat) = NULL;
  KINLINE bool operator() (const LObject &lhs, const LObject &rhs) const;
};

// Hash function for (poly, poly) pairs used in LSet pair index
struct PolyPairHash {
  std::size_t operator()(const std::pair<poly, poly>& p) const {
    // Combine the two pointer hashes
    std::size_t h1 = std::hash<poly>{}(p.first);
    std::size_t h2 = std::hash<poly>{}(p.second);
    // Simple combination: XOR and shift
    return h1 ^ (h2 << 1);
  }
};

class LSet : public writable_set<LObject, CompareLObject> {
private:
  unsigned seq = 0;   // increments by one on every insertion; used to determine ordering
  // Parallel flat array of sev_lcm values for cache-friendly scanning.
  // Same indexing as the flat_ vector in writable_set: sev_flat_[i] holds
  // the sev_lcm of the element at flat_[i].  Erased entries are set to 0
  // (sentinel: causes !(sev_p & ~0UL) == false, so erased entries always
  // fail the pre-filter).
  std::vector<unsigned long> sev_flat_;
  // Parallel flat array of sevSig values for cache-friendly signature scanning.
  // Same indexing as sev_flat_: sevSig_flat_[i] holds the sevSig of the
  // element at flat_[i].  Erased entries are set to 0 (sentinel).
  // In bba mode (non-signature), sevSig is uninitialized — that's fine,
  // the array contains junk that is never queried.
  std::vector<unsigned long> sevSig_flat_;

public:
  std::unordered_map<std::pair<poly, poly>, iterator, PolyPairHash> pair_index;

  // Helper to canonicalize a (p1,p2) pair by pointer value
  static std::pair<poly, poly> canonicalize_pair(poly p1, poly p2) {
    return (p1 <= p2) ? std::make_pair(p1, p2) : std::make_pair(p2, p1);
  }

  using writable_set<LObject, CompareLObject>::iterator;
  using writable_set<LObject, CompareLObject>::begin;
  using writable_set<LObject, CompareLObject>::end;
  using writable_set<LObject, CompareLObject>::rbegin;
  using writable_set<LObject, CompareLObject>::rend;
  using writable_set<LObject, CompareLObject>::key_comp;
  using writable_set<LObject, CompareLObject>::empty;
  using writable_set<LObject, CompareLObject>::size;
  using writable_set<LObject, CompareLObject>::size_type;
  using writable_set<LObject, CompareLObject>::unordered_iterator;
  using writable_set<LObject, CompareLObject>::ubegin;
  using writable_set<LObject, CompareLObject>::uend;

  // Filtered unordered iterator: scans a contiguous sev array
  // for cache-friendly pre-filtering, only visiting elements whose sev
  // passes the filter.  Two filter modes:
  //   sev2_=0: divisibility — skip where (sev1_ & ~sev_array[i]) != 0
  //   sev2_!=0: incomparability — skip where both
  //     (sev_array[i] & ~sev1_) and (sev2_ & ~sev_array[i]) are nonzero
  // The same iterator type is used for both sev_flat_ (lcm) and
  // sevSig_flat_ (signature) scans — only the array pointer differs.
  class filtered_iterator {
    LSet* owner_;
    size_t pos_;
    unsigned long sev1_;
    unsigned long sev2_;
    const std::vector<unsigned long>* sev_array_;
    friend class LSet;

    void advance() {
      // Always use sev_flat_ to detect deleted entries (sev_lcm sentinel 0).
      // Use sev_array_ (which may be sev_flat_ or sevSig_flat_) for the
      // actual filter check. This avoids the problem that sevSig can
      // legitimately be 0, which would be confused with the deleted sentinel.
      const unsigned long* del = owner_->sev_flat_.data();
      const unsigned long* sev = sev_array_->data();
      const size_t sz = sev_array_->size();
      while (pos_ < sz) {
        if (del[pos_] == 0) { ++pos_; continue; }     // deleted sentinel
        unsigned long s = sev[pos_];
        if (sev2_ == 0) {
          if (sev1_ & ~s) { ++pos_; continue; }       // divisibility: skip
        } else {
          if ((s & ~sev1_) && (sev2_ & ~s)) { ++pos_; continue; }  // incomp: skip
        }
        break;
      }
    }

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = LObject;
    using difference_type = std::ptrdiff_t;
    using pointer = LObject*;
    using reference = LObject&;

    filtered_iterator() : owner_(nullptr), pos_(0), sev1_(0), sev2_(0), sev_array_(nullptr) {}
    filtered_iterator(LSet* owner, size_t pos, unsigned long sev1, unsigned long sev2,
                      const std::vector<unsigned long>* sev_array)
      : owner_(owner), pos_(pos), sev1_(sev1), sev2_(sev2), sev_array_(sev_array) {
      advance();
    }

    reference operator*() const { return *owner_->flat_ptr(pos_); }
    pointer operator->() const { return owner_->flat_ptr(pos_); }

    filtered_iterator& operator++() { ++pos_; advance(); return *this; }
    filtered_iterator operator++(int) {
      filtered_iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const filtered_iterator& other) const { return pos_ == other.pos_; }
    bool operator!=(const filtered_iterator& other) const { return pos_ != other.pos_; }
  };

  // Filtered unordered iteration over sev_flat_ (sev_lcm): divisibility filter
  filtered_iterator ufbegin_lcm(unsigned long sev) {
    return filtered_iterator(this, 0, sev, 0, &sev_flat_);
  }
  // Filtered unordered iteration over sev_flat_ (sev_lcm): incomparability filter
  filtered_iterator ufbegin_lcm(unsigned long sev1, unsigned long sev2) {
    return filtered_iterator(this, 0, sev1, sev2, &sev_flat_);
  }
  // Sentinel for filtered lcm iteration (compares by pos_)
  filtered_iterator ufend_lcm() {
    return filtered_iterator(this, sev_flat_.size(), 0, 0, &sev_flat_);
  }
  // Filtered unordered iteration over sevSig_flat_ (signature): divisibility filter
  filtered_iterator ufbegin_sig(unsigned long sev) {
    return filtered_iterator(this, 0, sev, 0, &sevSig_flat_);
  }
  // Sentinel for filtered sig iteration (compares by pos_)
  filtered_iterator ufend_sig() {
    return filtered_iterator(this, sevSig_flat_.size(), 0, 0, &sevSig_flat_);
  }
  // Erase via filtered_iterator; returns next valid filtered position
  filtered_iterator erase(filtered_iterator it);

  // Direct access to sev_flat_ for index-based pair iteration (B dedup)
  const unsigned long* sev_flat_data() const { return sev_flat_.data(); }
  size_t sev_flat_size() const { return sev_flat_.size(); }

  // Default constructor
  LSet() = default;

  // Copy constructor: base class copy rebuilds flat_ with new indices,
  // so we must rebuild sev_flat_, sevSig_flat_, and pair_index to match.
  LSet(const LSet& other)
    : writable_set<LObject, CompareLObject>(other), seq(other.seq) {
    rebuild_sev_flat();
    rebuild_sevSig_flat();
    rebuild_pair_index();
  }

  // Move constructor: base class move rebuilds flat_ with new indices.
  LSet(LSet&& other) noexcept
    : writable_set<LObject, CompareLObject>(std::move(other)), seq(other.seq) {
    rebuild_sev_flat();
    rebuild_sevSig_flat();
    rebuild_pair_index();
  }

  // Copy assignment: same issue — base class rebuilds flat_ from scratch.
  LSet& operator=(const LSet& other) {
    if (this != &other) {
      writable_set<LObject, CompareLObject>::operator=(other);
      seq = other.seq;
      rebuild_sev_flat();
      rebuild_sevSig_flat();
      rebuild_pair_index();
    }
    return *this;
  }

  // Move assignment
  LSet& operator=(LSet&& other) noexcept {
    if (this != &other) {
      writable_set<LObject, CompareLObject>::operator=(std::move(other));
      seq = other.seq;
      rebuild_sev_flat();
      rebuild_sevSig_flat();
      rebuild_pair_index();
    }
    return *this;
  }

  // Override insert to maintain pair_index, sev_flat_, and sevSig_flat_
  iterator insert(const LObject& lobject) {
    iterator it = writable_set<LObject, CompareLObject>::insert(lobject);
    sev_flat_.push_back(it->sev_lcm);
    sevSig_flat_.push_back(it->sevSig);
    if (it->p1 != NULL && it->p2 != NULL) {
      auto key = canonicalize_pair(it->p1, it->p2);
      pair_index.emplace(key, it);
    }
    return it;
  }

  // Override clear to also clear pair_index, sev_flat_, and sevSig_flat_
  void clear() {
    pair_index.clear();
    sev_flat_.clear();
    sevSig_flat_.clear();
    writable_set<LObject, CompareLObject>::clear();
  }

  // Override reorder to rebuild pair_index, sev_flat_, and sevSig_flat_
  void reorder() {
    writable_set<LObject, CompareLObject>::reorder();
    // After reorder, flat_ has no gaps and flat_index is reassigned 0..size()-1
    rebuild_sev_flat();
    rebuild_sevSig_flat();
    pair_index.clear();
    for (iterator it = begin(); it != end(); ++it) {
      if (it->p1 != NULL && it->p2 != NULL) {
        auto key = canonicalize_pair(it->p1, it->p2);
        pair_index.emplace(key, it);
      }
    }
  }

  // Rebuild sev_flat_ from the current flat_ array.
  // Deleted entries get sentinel 0, valid entries get their sev_lcm.
  void rebuild_sev_flat() {
    const size_t n = this->flat_size();
    sev_flat_.resize(n);
    for (size_t i = 0; i < n; i++) {
      LObject* p = this->flat_ptr(i);
      sev_flat_[i] = (p != NULL) ? p->sev_lcm : 0;
    }
  }
  // Rebuild sevSig_flat_ from the current flat_ array.
  // Deleted entries get sentinel 0, valid entries get their sevSig.
  void rebuild_sevSig_flat() {
    const size_t n = this->flat_size();
    sevSig_flat_.resize(n);
    for (size_t i = 0; i < n; i++) {
      LObject* p = this->flat_ptr(i);
      sevSig_flat_[i] = (p != NULL) ? p->sevSig : 0;
    }
  }
  // Rebuild pair_index from all live elements.
  void rebuild_pair_index() {
    pair_index.clear();
    for (iterator it = begin(); it != end(); ++it) {
      if (it->p1 != NULL && it->p2 != NULL) {
        auto key = canonicalize_pair(it->p1, it->p2);
        pair_index.emplace(key, it);
      }
    }
  }

  KINLINE iterator push(LObject& lobject);
  KINLINE bool would_be_top(LObject& lobject);
  KINLINE void pop(void);
  KINLINE void pop_and_erase(void);
  KINLINE const LObject& top(void);
  iterator erase(iterator it);
  unordered_iterator erase(unordered_iterator it);
};

// SElement: one element of the S set (standard basis).
struct SElement {
  poly p;              // the polynomial
  int ecart;           // ecart
  unsigned long sev;   // short exponent vector
  int s_2_r;           // index into R array
  int length;          // number of terms (replaces lenS)
  wlen_type wlength;   // weighted length (replaces lenSw)
  int fromQ;           // from quotient ideal
  poly sig;            // signature (sba only)
  unsigned long sevSig;// short exponent vector of signature (sba only)
  bool deleted;              // lazy-delete flag for parallel phase
  mutable bool pairtest;     // transient: true if spoly(this, h)==0 during enterOnePair

  SElement() : p(NULL), ecart(0), sev(0), s_2_r(0), length(0), wlength(0),
               fromQ(0), sig(NULL), sevSig(0), deleted(false), pairtest(false) {}
};

/**
 * @class sBasisSet
 * @brief The S-set (standard basis) for Groebner basis computation.
 *
 * Inherits privately from BlockArray<SElement> — all element access is
 * through iterators, not integer indices.
 *
 * Ordering modes control insert() and erase() behavior:
 *   SORDER_STANDARD:  sorted by leading monomial + ecart (posInS).
 *                     insert() finds position via binary search, shifts.
 *                     erase() shifts elements down.
 *   SORDER_MONFIRST:  monomials first, then by degree (posInSMonFirst).
 *                     Used by sba() over rings.
 *   SORDER_APPEND:    append at end, no sorting (parallel phase).
 *                     erase() sets deleted flag (no shifting).
 */
enum SOrderMode { SORDER_STANDARD, SORDER_MONFIRST, SORDER_APPEND };

class sBasisSet : private BlockArray<SElement> {
public:
  // --- Iterator (skips deleted entries) ---
  class iterator {
    sBasisSet* set_;
    int pos_;
    friend class sBasisSet;

    void skip_deleted_forward() {
      while (pos_ < set_->count && set_->elem(pos_).deleted) pos_++;
    }
    void skip_deleted_backward() {
      while (pos_ > 0 && set_->elem(pos_).deleted) pos_--;
    }

    iterator(sBasisSet* s, int pos) : set_(s), pos_(pos) { }

  public:
    iterator() : set_(NULL), pos_(0) {}

    SElement& operator*() const { return set_->elem(pos_); }
    SElement* operator->() const { return &set_->elem(pos_); }

    // Raw index — available for internal use (e.g., pairtest migration).
    // Prefer iterator-based APIs in new code.
    int index() const { return pos_; }

    iterator& operator++() {
      pos_++;
      skip_deleted_forward();
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator& operator--() {
      pos_--;
      skip_deleted_backward();
      return *this;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      --(*this);
      return tmp;
    }

    bool operator==(const iterator& other) const { return pos_ == other.pos_; }
    bool operator!=(const iterator& other) const { return pos_ != other.pos_; }
  };

  class const_iterator {
    const sBasisSet* set_;
    int pos_;
    friend class sBasisSet;

    void skip_deleted_forward() {
      while (pos_ < set_->count && set_->elem(pos_).deleted) pos_++;
    }

    const_iterator(const sBasisSet* s, int pos) : set_(s), pos_(pos) { }

  public:
    const_iterator() : set_(NULL), pos_(0) {}
    const_iterator(const iterator& it) : set_(it.set_), pos_(it.pos_) {}

    const SElement& operator*() const { return set_->elem(pos_); }
    const SElement* operator->() const { return &set_->elem(pos_); }
    int index() const { return pos_; }

    const_iterator& operator++() {
      pos_++;
      skip_deleted_forward();
      return *this;
    }
    const_iterator operator++(int) {
      const_iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const const_iterator& other) const { return pos_ == other.pos_; }
    bool operator!=(const const_iterator& other) const { return pos_ != other.pos_; }
  };

  // --- Construction / mode ---
  sBasisSet() : order_(SORDER_STANDARD), live_count_(0), pairtest_any_(false) {}

  SOrderMode order() const { return order_; }
  void set_order(SOrderMode m) { order_ = m; }
  bool is_append_mode() const { return order_ == SORDER_APPEND; }

  // --- Iterators ---
  iterator begin() {
    iterator it(this, 0);
    it.skip_deleted_forward();
    return it;
  }
  iterator end() { return iterator(this, count); }

  const_iterator begin() const {
    const_iterator it(this, 0);
    it.skip_deleted_forward();
    return it;
  }
  const_iterator end() const { return const_iterator(this, count); }

  // --- Size ---
  int size() const { return live_count_; }
  bool empty() const { return live_count_ == 0; }

  // --- Insert ---
  // Inserts val at the position determined by the current ordering mode.
  // SORDER_STANDARD/SORDER_MONFIRST: sorted insert (needs strat for comparison).
  // SORDER_APPEND: append at end.
  // Implemented in kutil.cc.
  iterator insert(const SElement& val, kStrategy strat);

  // Append at end regardless of mode (used during initialization).
  iterator push_back(const SElement& val) {
    BlockArray<SElement>::push_back(val);
    live_count_++;
    return iterator(this, count - 1);
  }

  // --- Erase ---
  // SORDER_APPEND: set deleted flag (no shifting).
  // Other modes: shift elements down.
  void erase(iterator it) {
    if (order_ == SORDER_APPEND) {
      elem(it.pos_).deleted = true;
    } else {
      BlockArray<SElement>::erase(it.pos_);
    }
    live_count_--;
  }

  // --- Capacity / storage ---
  using BlockArray<SElement>::ensure_capacity;
  using BlockArray<SElement>::capacity;
  using BlockArray<SElement>::free_all;

  // Set element count directly (used during initialization).
  void setsize(int n) {
    BlockArray<SElement>::setsize(n);
    live_count_ = n;
  }

  // --- Pairtest ---
  // Clear all pairtest flags and the sentinel.
  void clear_pairtest() {
    for (int i = 0; i < count; i++) elem(i).pairtest = false;
    pairtest_any_ = false;
  }
  // Set the sentinel (some zero spoly was found).
  void set_pairtest_any() { pairtest_any_ = true; }
  bool has_pairtest() const { return pairtest_any_; }

  // Compact: remove deleted entries, pack remaining entries contiguously.
  // Only meaningful after lazy mode. Resets to non-lazy mode.
  void compact() {
    int dst = 0;
    for (int src = 0; src < count; src++) {
      if (!elem(src).deleted) {
        if (dst != src)
          elem(dst) = elem(src);
        dst++;
      }
    }
    count = dst;
    live_count_ = dst;
    order_ = SORDER_STANDARD;
  }

  // --- Stable pointer access (valid as long as element exists) ---
  SElement* addr(iterator it) {
    return BlockArray<SElement>::addr(it.pos_);
  }

  // --- Member methods (use private elem() for random access) ---
  // Implementations in kutil.cc since they need skStrategy, which is
  // defined after this class.

  // Insert a new basis element from an LObject (replaces enterSBba).
  // Builds an SElement, inserts at the position determined by ordering mode.
  void enter_bba(LObject &p, kStrategy strat, int atR = -1, int atS = -1);

  // Insert for signature-based algorithms (replaces enterSSba).
  // Also copies sig and sevSig fields.
  void enter_sba(LObject &p, kStrategy strat, int atR = -1, int atS = -1);

  // Delete element at iterator and update related structures.
  // In non-lazy mode, shifts elements down (like old deleteInS).
  // Returns next valid iterator after the erased position.
  iterator erase_and_next(iterator it);

  // Reorder S after inter-reduction may have changed leading monomials.
  // (replaces reorderS)
  void reorder(int *suc, kStrategy strat);

  // S-to-T lookups (replace skStrategy::S_2_T / s_2_t).
  // Return the TObject corresponding to the S element at the iterator.
  TObject* S_2_T(const_iterator it, kStrategy strat);
  TObject* s_2_t(const_iterator it, kStrategy strat);

  // clearS: if p divides *at, delete *at and adjust the iterator.
  // (replaces the clearS inline in kInline.h)
  void clear_if_divisible(poly p, unsigned long p_sev,
                          iterator &at, kStrategy strat);

  // For tgb.cc: binary search variant with different comparison.
  iterator simple_find_pos(kStrategy strat, poly p);

  // For tgb.cc: move an element from one position to another,
  // shifting intervening elements. old_it and new_it provide the
  // element positions.
  void move_elem(iterator old_it, iterator new_it);

  // Construct an iterator at a raw index (no skip).
  // Transitional: used by code migrating from int-based to iterator-based API.
  // In lazy mode, the element at index i may be deleted — caller must check.
  iterator iterator_at(int i) { return iterator(this, i); }
  const_iterator const_iterator_at(int i) const { return const_iterator(this, i); }

private:
  SOrderMode order_;
  int live_count_;
  bool pairtest_any_;  // sentinel: true if any SElement.pairtest was set

  // Internal: insert at a specific position (for reorder, etc.)
  iterator insert_at(int pos, const SElement& val) {
    BlockArray<SElement>::insert(pos, val);
    live_count_++;
    return iterator(this, pos);
  }

  // Internal: binary search for sorted insertion position (posInS).
  iterator find_pos(kStrategy strat, const poly p, int ecart_p);

  // Internal: binary search for monfirst insertion position (posInSMonFirst).
  iterator find_pos_monfirst(kStrategy strat, const poly p);

  // Internal raw element access (used by iterators and member functions)
  SElement& elem(int i) { return (*static_cast<BlockArray<SElement>*>(this))[i]; }
  const SElement& elem(int i) const { return (*static_cast<const BlockArray<SElement>*>(this))[i]; }

  friend class skStrategy;
};

class skStrategy
#ifdef HAVE_OMALLOC
                 : public omallocClass
#endif
{
public:
  kStrategy next = NULL;
  int (*red)(LObject * L,kStrategy strat) = NULL;
  int (*red2)(LObject * L,kStrategy strat) = NULL;
  void (*initEcart)(TObject * L) = NULL;
  int (*posInT)(const BlockArray<TObject> &T,const int tl,LObject &h) = NULL;
  int (*compareL) (const LObject &lhs, const LObject &rhs, const kStrategy strat) = NULL;
  int (*compareLOld) (const LObject &lhs, const LObject &rhs, const kStrategy strat) = NULL;
  void (*enterS)(LObject &h, kStrategy strat, int atR /*= -1*/, int atS /*= -1*/) = NULL;
  void (*initEcartPair)(LObject * h, poly f, poly g, int ecartF, int ecartG) = NULL;
  void (*enterOnePair) (const SElement &si,poly p,int ecart, int isFromQ,kStrategy strat, int atR /*= -1*/) = NULL;
  void (*chainCrit) (poly p,int ecart,kStrategy strat) = NULL;
  BOOLEAN (*syzCrit) (poly sig, unsigned long not_sevSig, kStrategy strat) = NULL;
  BOOLEAN (*rewCrit1) (poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start /*= 0*/) = NULL;
  BOOLEAN (*rewCrit2) (poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start /*= 0*/) = NULL;
  BOOLEAN (*rewCrit3) (poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start /*= 0*/) = NULL;
  pFDegProc pOrigFDeg = NULL;
  pLDegProc pOrigLDeg = NULL;
  pFDegProc pOrigFDeg_TailRing = NULL;
  pLDegProc pOrigLDeg_TailRing = NULL;
  s_poly_proc_t s_poly = NULL;

  LObject P;
  ideal D = NULL; /*V(S) is in D(D)*/
  ideal M = NULL; /*set of minimal generators*/
  sBasisSet S;
  int Srank = 1;       // rank for getShdl() result ideal
  polyset syz = NULL;
  intset syzIdx = NULL;// index in the syz array at which the first
                       // syzygy of component i comes up
                       // important for signature-based algorithms
  unsigned sbaOrder = 0;
  int currIdx = 0;
  int max_lower_index = 0;
  unsigned long* sevSyz = NULL;

  // Build an ideal from S on demand. The caller owns the returned ideal.
  // The ideal's m[] entries share poly pointers with S[i].p — they are
  // NOT copies. The caller must not free them while S is still alive,
  // or must take ownership (set S[i].p = NULL afterward).
  ideal getShdl();
  BlockArray<unsigned long> sevT;
  BlockArray<TObject> T;
  LSet L;
  LSet    B;
  poly    kNoether = NULL;
  poly    t_kNoether = NULL; // same polys in tailring
  KINLINE poly    kNoetherTail();
  BOOLEAN * NotUsedAxis = NULL;
  // pairtest is now per-SElement (SElement::pairtest) + sBasisSet::pairtest_any_
  poly tail = NULL;
  intvec * kModW = NULL;
  intvec * kHomW = NULL;
  // procedure for ShalloCopy from tailRing  to currRing
  pShallowCopyDeleteProc p_shallow_copy_delete = NULL;
  // pointers to Tobjects R[i] is ith Tobject which is generated
  BlockArray<TObject*>  R;
  ring tailRing = NULL;
  omBin lmBin = NULL;
  omBin tailBin = NULL;
  int nr = 0;
  int cp = 0,c3 = 0;
  int mu = 0;
  int syzl = 0,syzmax = 0,syzidxmax = 0;
  int ak = 0,LazyDegree = 0,LazyPass = 0;
  int syzComp = 0;
  int lastAxis = 0;
  int newIdeal = 0;
  int minim = 0;
  bool sigdrop = false; //This is used to check sigdrop in sba over Z
  int nrsyzcrit = 0; // counts how many pairs are deleted by SyzCrit
  int nrrewcrit = 0; // counts how many pairs are deleted by FaugereRewCrit
  int sbaEnterS = 0; // sba over Z strategy: if sigdrop element has _*gen(sbaEnterS+1), then
                     // add directly sbaEnterS elements into S
  int blockred = 0;  // counter for blocked reductions in redSig
  int blockredmax = 0;
  #ifdef HAVE_SHIFTBBA
  int cv = 0; // in shift bases: counting V criterion
  /*BOOLEAN*/ char rightGB = FALSE;
  #endif
  /*BOOLEAN*/ char interpt = FALSE;
  /*BOOLEAN*/ char homog = FALSE;
#ifdef HAVE_PLURAL
  /*BOOLEAN*/ char z2homog = FALSE; // Z_2 - homogeneous input allows product criterion in commutative and SCA cases!
#endif
  /*BOOLEAN*/ char kAllAxis = FALSE; // all axis are used -> (re)compute noether
  /*BOOLEAN*/ char honey = FALSE,sugarCrit = FALSE;
  /*BOOLEAN*/ char Gebauer = FALSE,noTailReduction = FALSE;
  /*BOOLEAN*/ char fromT = FALSE;
  /*BOOLEAN*/ char noetherSet = FALSE;
  /*BOOLEAN*/ char update = FALSE;
  /*BOOLEAN*/ char compareLOldFlag = FALSE;
  /*BOOLEAN*/ char use_buckets = FALSE;
  // if set, pLDeg(p, l) == (pFDeg(pLast(p), pLength)
  /*BOOLEAN*/ char LDegLast = FALSE;
  // if set, then L.length == L.pLength
  /*BOOLEAN*/ char length_pLength = FALSE;
  // if set, then compareInL does not depend on L.length
  /*BOOLEAN*/ char compareLDependsOnLength = FALSE;
#ifdef HAVE_PLURAL
  // set this flag to 1 to stop the product criteria
  // use ALLOW_PROD_CRIT(strat) to test
  /*BOOLEAN*/ char no_prod_crit = FALSE;
#define ALLOW_PROD_CRIT(A) (!(A)->no_prod_crit)
#else
#define ALLOW_PROD_CRIT(A) (1)
#endif
  char    redTailChange = FALSE;
  char    news = FALSE;
  char    newt = FALSE;/*used for messageSets*/
  char    noClearS = FALSE;
  char    completeReduce_retry = FALSE;
  char    overflow = FALSE;
  // Flags indicating whether optional S-set fields are in use.
  // With SElement, the fields always exist but may not be meaningful.
  /*BOOLEAN*/ char use_lenS = FALSE;   // replaces lenS != NULL check
  /*BOOLEAN*/ char use_lenSw = FALSE;  // replaces lenSw != NULL check
  /*BOOLEAN*/ char hasFromQ = FALSE;   // replaces fromQ != NULL check

  skStrategy();
  ~skStrategy();

  // return TObject corresponding to S[i]: assume that it exists
  // i.e. no error checking is done
  KINLINE TObject* S_2_T(int i);
  // like S_2_T, except that NULL is returned if it can not be found
  KINLINE TObject* s_2_t(int i);
};

int compareL0 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL0Ring (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareLSig (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareLSigRing (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL10 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL11 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL11Ring (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareLF5C (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareLF5CRing (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL11Ringls (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL110 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL110Ring (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL13 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL15 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL15Ring (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL17 (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL17Ring (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL17_c (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareL17_cRing (const LObject &lhs, const LObject &rhs, const kStrategy strat);
int compareLSpecial (const LObject &lhs, const LObject &rhs, const kStrategy strat);

void deleteHC(poly *p, int *e, int *l, kStrategy strat);
void deleteHC(LObject* L, kStrategy strat, BOOLEAN fromNext = FALSE);
void deleteInS (int i,kStrategy strat);
void cleanT (kStrategy strat);
void enterSBba (LObject &p, kStrategy strat, int atR = -1, int atS = -1);
void enterSBbaShift (LObject &p, kStrategy strat, int atR = -1, int atS = -1);
void enterSSba (LObject &p, kStrategy strat, int atR = -1, int atS = -1);
void initEcartPairBba (LObject* Lp,poly f,poly g,int ecartF,int ecartG);
void initEcartPairMora (LObject* Lp,poly f,poly g,int ecartF,int ecartG);
int posInS (const kStrategy strat, const int length, const poly p,
            const int ecart_p);
int posInSMonFirst (const kStrategy strat, const int length, const poly p);
int posInIdealMonFirst (const ideal F, const poly p,int start = 0,int end = -1);
int posInT0 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT1 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT2 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT11 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInTSig (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT110 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT13 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT15 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT17 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT17_c (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT19 (const BlockArray<TObject> &set,const int length,LObject &p);
int posInT_EcartpLength(const BlockArray<TObject> &set,const int length,LObject &p);
int posInT_EcartFDegpLength(const BlockArray<TObject> &set,const int length,LObject &p);
int posInT_FDegpLength(const BlockArray<TObject> &set,const int length,LObject &p);
int posInT_pLength(const BlockArray<TObject> &set,const int length,LObject &p);

#ifdef HAVE_MORE_POS_IN_T
int posInT_EcartFDegpLength(const BlockArray<TObject> &set,const int length,LObject &p);
int posInT_FDegpLength(const BlockArray<TObject> &set,const int length,LObject &p);
int posInT_pLength(const BlockArray<TObject> &set,const int length,LObject &p);
#endif


void reorderS (int* suc,kStrategy strat);
int posInSyz (const kStrategy strat, const poly sig);
KINLINE poly redtailBba (poly p,int end_pos,kStrategy strat,BOOLEAN normalize=FALSE);
KINLINE poly redtailBbaBound (poly p,int end_pos,kStrategy strat,int bound,BOOLEAN normalize=FALSE);
KINLINE poly redtailBba_Ring (poly p,int end_pos,kStrategy strat);
KINLINE poly redtailBba_Z (poly p,int end_pos,kStrategy strat);
poly redtailBba_NF (poly p, kStrategy strat );
poly redtailBba_Ring (LObject* L, int end_pos, kStrategy strat );
poly redtailBba_Z (LObject* L, int end_pos, kStrategy strat );
void redtailBbaAlsoLC_Z (LObject* L, int end_pos, kStrategy strat );
poly redtailBba (LObject *L, int end_pos,kStrategy strat,
                 BOOLEAN withT = FALSE,BOOLEAN normalize=FALSE);
poly redtailBbaBound (LObject *L, int end_pos,kStrategy strat,int bound,
                 BOOLEAN withT = FALSE,BOOLEAN normalize=FALSE);
poly redtailSba (LObject *L, int end_pos,kStrategy strat,
                 BOOLEAN withT = FALSE,BOOLEAN normalize=FALSE);
poly redtailBba (TObject *T, int end_pos,kStrategy strat);
poly redtail (poly p,int end_pos,kStrategy strat);
poly redtail (LObject *L,int end_pos,kStrategy strat);
poly redNF (poly h,int & max_ind,int nonorm,kStrategy strat);
int redNF0 (LObject *P,kStrategy strat);
poly redNFTail (poly h,const int sl,kStrategy strat);
int redHoney (LObject* h, kStrategy strat);
int redHoneyM (LObject* h, kStrategy strat);
int redLiftstd (LObject* h, kStrategy strat);
int redRing (LObject* h,kStrategy strat);
int redRing_Z (LObject* h,kStrategy strat);
int redRiloc (LObject* h,kStrategy strat);
void enterExtendedSpoly(poly h,kStrategy strat);
void enterExtendedSpolySig(poly h,poly hSig,kStrategy strat);
void superenterpairs (poly h,int k,int ecart,int pos,kStrategy strat, int atR = -1);
void superenterpairsSig (poly h,poly hSig,int hFrom,int k,int ecart,int pos,kStrategy strat, int atR = -1);
int redLazy (LObject* h,kStrategy strat);
int redHomog (LObject* h,kStrategy strat);
int redSig (LObject* h,kStrategy strat);
int redSigRing (LObject* h,kStrategy strat);
//adds hSig to be able to check with F5's criteria when entering pairs!
void enterpairsSig (poly h, poly hSig, int from, int k, int ec, int pos,kStrategy strat, int atR = -1);
void enterpairs (poly h, int k, int ec, int pos,kStrategy strat, int atR = -1);
void entersets (LObject h);
void pairs ();
BOOLEAN sbaCheckGcdPair (LObject* h,kStrategy strat);
void message (int i,int* olddeg,LSet::size_type* reduc,kStrategy strat,int red_result);
void messageStat (int hilbcount,kStrategy strat);
void messageStatSBA (int hilbcount,kStrategy strat);
#ifdef KDEBUG
void messageSets (kStrategy strat);
#else
#define messageSets(s)  do {} while (0)
#endif

void initEcartNormal (TObject* h);
void initEcartBBA (TObject* h);
void initS (ideal F, ideal Q,kStrategy strat);
void initSL (ideal F, ideal Q,kStrategy strat);
void initSLSba (ideal F, ideal Q,kStrategy strat);
/*************************************************
 * when initializing a new bunch of principal
 * syzygies at the beginning of a new iteration
 * step in a signature-based algorithm we
 * compute ONLY the leading elements of those
 * syzygies, NOT the whole syzygy
 * NOTE: this needs to be adjusted for a more
 * general approach on signature-based algorithms
 ***********************************************/
void initSyzRules (kStrategy strat);
void updateS(BOOLEAN toT,kStrategy strat);
void enterSyz (LObject &p,kStrategy strat, int atT);
void enterT (LObject &p,kStrategy strat, int atT = -1);
// enlargeT removed — BlockArray grows automatically via ensure_capacity
void replaceInLAndSAndT(LObject &p, int tj, kStrategy strat);
void enterT_strong (LObject &p,kStrategy strat, int atT = -1);
void cancelunit (LObject* p,BOOLEAN inNF=FALSE);
void HEckeTest (poly pp,kStrategy strat);
void initBuchMoraCrit(kStrategy strat);
void initSbaCrit(kStrategy strat);
void initHilbCrit(ideal F, ideal Q, bigintmat **hilb,kStrategy strat);
void initBuchMoraPos(kStrategy strat);
void initBuchMoraPosRing(kStrategy strat);
void initSbaPos(kStrategy strat);
void initBuchMora (ideal F, ideal Q,kStrategy strat);
void initSbaBuchMora (ideal F, ideal Q,kStrategy strat);
void exitBuchMora (kStrategy strat);
void exitSba (kStrategy strat);
void updateResult(ideal Q,kStrategy strat);
void completeReduce (kStrategy strat, BOOLEAN withT=FALSE);
void kFreeStrat(kStrategy strat);
void enterOnePairNormal (const SElement &si,poly p,int ecart, int isFromQ,kStrategy strat, int atR);
void chainCritNormal (poly p,int ecart,kStrategy strat);
void chainCritOpt_1 (poly,int,kStrategy strat);
void chainCritSig (poly p,int ecart,kStrategy strat);
BOOLEAN homogTest(polyset F, int Fmax);
BOOLEAN newHEdge(kStrategy strat);
BOOLEAN syzCriterion(poly sig, unsigned long not_sevSig, kStrategy strat);
BOOLEAN syzCriterionInc(poly sig, unsigned long not_sevSig, kStrategy strat);
KINLINE BOOLEAN arriRewDummy(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start);
BOOLEAN arriRewCriterion(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start);
BOOLEAN arriRewCriterionPre(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start);
BOOLEAN faugereRewCriterion(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, int start);
BOOLEAN findMinLMPair(poly sig, unsigned long not_sevSig, kStrategy strat, int start);

/// returns index of p in TSet, or -1 if not found
int kFindInT(poly p, const BlockArray<TObject> &T, int tlength);
#ifdef HAVE_SHIFTBBA
int kFindInTShift(poly p, const BlockArray<TObject> &T, int tlength);
#endif

/// return -1 if no divisor is found
///        number of first divisor in T, otherwise
int kFindDivisibleByInT(const kStrategy strat, const LObject* L, const int start=0);
int kFindDivisibleByInT_ecart(const kStrategy strat, const LObject* L, const int ecart);
int kFindDivisibleByInT_Z(const kStrategy strat, const LObject* L, const int start=0);
int kFindSameLMInT_Z(const kStrategy strat, const LObject* L, const int start=0);

/// tests if T[0] divides the leading monomial of L, returns -1 if not
int kTestDivisibleByT0_Z(const kStrategy strat, const LObject* L);
/// return -1 if no divisor is found
///        number of first divisor in S, otherwise
int kFindDivisibleByInS(const kStrategy strat, int *max_ind, LObject* L);

int kFindNextDivisibleByInS(const kStrategy strat, int start,int max_ind, LObject* L);
TObject* kFindDivisibleByInS_T(kStrategy strat, int end_pos, LObject* L, TObject *T, long ecart = LONG_MAX);

/***************************************************************
 *
 * stuff to be inlined
 *
 ***************************************************************/

KINLINE void initT (BlockArray<TObject> &T);
KINLINE void initR (BlockArray<TObject*> &R);
KINLINE void initsevT (BlockArray<unsigned long> &sevT);
KINLINE poly k_LmInit_currRing_2_tailRing(poly p, ring tailRing, omBin bin);
KINLINE poly k_LmInit_tailRing_2_currRing(poly p, ring tailRing, omBin bin);
KINLINE poly k_LmShallowCopyDelete_currRing_2_tailRing(poly p, ring tailRing, omBin bin);
KINLINE poly k_LmShallowCopyDelete_tailRing_2_currRing(poly p, ring tailRing,  omBin bin);

KINLINE poly k_LmInit_currRing_2_tailRing(poly p, ring tailRing);
KINLINE poly k_LmInit_tailRing_2_currRing(poly p, ring tailRing);
KINLINE poly k_LmShallowCopyDelete_currRing_2_tailRing(poly p, ring tailRing);
KINLINE poly k_LmShallowCopyDelete_tailRing_2_currRing(poly p, ring tailRing);

// if exp bound is not violated, return TRUE and
//                               get m1 = LCM(LM(p1), LM(p2))/LM(p1)
//                                   m2 = LCM(LM(p1), LM(p2))/LM(p2)
// return FALSE and m1 == NULL, m2 == NULL     , otherwise
KINLINE BOOLEAN k_GetLeadTerms(const poly p1, const poly p2, const ring p_r,
                               poly &m1, poly &m2, const ring m_r);
KINLINE void k_GetStrongLeadTerms(const poly p1, const poly p2, const ring leadRing,
                               poly &m1, poly &m2, poly &lcm, const ring taiRing);
#ifdef KDEBUG
// test strat
BOOLEAN kTest(kStrategy strat);
// test strat, and test that S is contained in T
BOOLEAN kTest_TS(kStrategy strat);
// test LObject
BOOLEAN kTest_L(LObject* L, kStrategy strat,
                 BOOLEAN testp = FALSE, int lpos = -1,
                 BlockArray<TObject> *T = NULL, int tlength = -1);
// test TObject
BOOLEAN kTest_T(TObject* T, kStrategy strat, int tpos = -1, char TN = '?');
// test set strat->SevS
BOOLEAN kTest_S(kStrategy strat);
#else
#define kTest(A)        (TRUE)
#define kTest_TS(A)     (TRUE)
#define kTest_T(T,S)    (TRUE)
#define kTest_S(T)      (TRUE)
#define kTest_L(T,R)    (TRUE)
#endif


/***************************************************************
 *
 * From kstd1.cc
 *
 ***************************************************************/
int redFirst (LObject* h,kStrategy strat);
int redEcart (LObject* h,kStrategy strat);
void enterSMora (LObject &p, kStrategy strat, int atR = -1, int atS = -1);
void enterSMoraNF (LObject &p, kStrategy strat, int atR = -1, int atS = -1);


/***************************************************************
 *
 * From kstd2.cc
 *
 ***************************************************************/
poly kFindZeroPoly(poly input_p, ring leadRing, ring tailRing);
ideal bba (ideal F, ideal Q,intvec *w,bigintmat *hilb,kStrategy strat);
ideal sba (ideal F, ideal Q,intvec *w,bigintmat *hilb,kStrategy strat);
poly kNF2 (ideal F, ideal Q, poly q, kStrategy strat, int lazyReduce);
ideal kNF2 (ideal F,ideal Q,ideal q, kStrategy strat, int lazyReduce);
poly kNF2Bound (ideal F, ideal Q, poly q,int bound, kStrategy strat, int lazyReduce);
ideal kNF2Bound (ideal F,ideal Q,ideal q,int bound, kStrategy strat, int lazyReduce);
void initBba(kStrategy strat);
void initSba(ideal F,kStrategy strat);
void f5c (kStrategy strat, int& olddeg, int& minimcnt, int& hilbeledeg,
          int& hilbcount, int& srmax, LSet::size_type& reduc, ideal Q,
          intvec *w,bigintmat *hilb );

/***************************************************************
 *
 * From kspoly.cc
 *
 ***************************************************************/
// Reduces PR with PW
// Assumes PR != NULL, PW != NULL, Lm(PW) divides Lm(PR)
// Changes: PR
// Const:   PW
// If coef != NULL, then *coef is a/gcd(a,b), where a = LC(PR), b = LC(PW)
// If strat != NULL, tailRing is changed if reduction would violate exp bound
// of tailRing
// Returns: 0 everything ok, no tailRing change
//          1 tailRing has successfully changed (strat != NULL)
//          2 no reduction performed, tailRing needs to be changed first
//            (strat == NULL)
//         -1 tailRing change could not be performed due to exceeding exp
//            bound of currRing
//  reduce should be set inside "kNF" (only relevant for rings)
int ksReducePoly(LObject* PR,
                 TObject* PW,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 poly *mon =NULL,
                 kStrategy strat = NULL,
                 BOOLEAN redtail = FALSE);

/* like ksReducePoly, but if the reducer has only 1 term we still
 * compute a possible coefficient multiplier for PR. this comes from
 * a special situation in redRing_Z and it is used only there. */
int ksReducePolyZ(LObject* PR,
                 TObject* PW,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 kStrategy strat = NULL);

int ksReducePolyLC(LObject* PR,
                 TObject* PW,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 kStrategy strat = NULL);


int ksReducePolyGCD(LObject* PR,
                 TObject* PW,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 kStrategy strat = NULL);

int ksReducePolyBound(LObject* PR,
                 TObject* PW,
                 int bound,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 kStrategy strat = NULL);

// Reduces PR with PW
// Assumes PR != NULL, PW != NULL, Lm(PW) divides Lm(PR)
// Changes: PR
// Const:   PW
// If coef != NULL, then *coef is a/gcd(a,b), where a = LC(PR), b = LC(PW)
// If strat != NULL, tailRing is changed if reduction would violate exp bound
// of tailRing
// Returns: 0 everything ok, no tailRing change
//          1 tailRing has successfully changed (strat != NULL)
//          2 no reduction performed, tailRing needs to be changed first
//            (strat == NULL)
//          3 no reduction performed, not sig-safe!!!
//         -1 tailRing change could not be performed due to exceeding exp
//            bound of currRing
int ksReducePolySig(LObject* PR,
                 TObject* PW,
                 long idx,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 kStrategy strat = NULL);

int ksReducePolySigRing(LObject* PR,
                 TObject* PW,
                 long idx,
                 poly spNoether = NULL,
                 number *coef = NULL,
                 kStrategy strat = NULL);

// Reduces PR at Current->next with PW
// Assumes PR != NULL, Current contained in PR
//         Current->next != NULL, LM(PW) divides LM(Current->next)
// Changes: PR
// Const:   PW
// Return: see ksReducePoly
int ksReducePolyTail(LObject* PR,
                     TObject* PW,
                     poly Current,
                     poly spNoether = NULL);

KINLINE int ksReducePolyTail(LObject* PR, TObject* PW, LObject* Red);

// Creates S-Poly of Pair
// Const:   Pair->p1, Pair->p2
// Changes: Pair->p == S-Poly of p1, p2
// Assume:  Pair->p1 != NULL && Pair->p2
void ksCreateSpoly(LObject* Pair, poly spNoether = NULL,
                   int use_buckets=0, ring tailRing=currRing,
                   poly m1 = NULL, poly m2 = NULL, BlockArray<TObject*>* R = NULL);

/*2
* creates the leading term of the S-polynomial of p1 and p2
* do not destroy p1 and p2
* remarks:
*   1. the coefficient is 0 (nNew)
*   2. pNext is undefined
*/
poly ksCreateShortSpoly(poly p1, poly p2, ring tailRing);


// old stuff
KINLINE poly ksOldSpolyRed(poly p1, poly p2, poly spNoether = NULL);
KINLINE poly ksOldSpolyRedNew(poly p1, poly p2, poly spNoether = NULL);
KINLINE poly ksOldCreateSpoly(poly p1, poly p2, poly spNoether = NULL, ring r = currRing);
KINLINE void ksOldSpolyTail(poly p1, poly q, poly q2, poly spNoether, ring r = currRing);

/***************************************************************
 *
 * Routines related for ring changes during std computations
 *
 ***************************************************************/
// return TRUE and set m1, m2 to k_GetLcmTerms,
//             if spoly creation of strat->P does not violate
//             exponent bound of strat->tailRing
//      FALSE, otherwise
BOOLEAN kCheckSpolyCreation(LObject* L, kStrategy strat, poly &m1, poly &m2);
// return TRUE if gcdpoly creation of R[atR] and S[atS] does not violate
//             exponent bound of strat->tailRing
//      FALSE, otherwise
BOOLEAN kCheckStrongCreation(int atR, poly m1, int atS, poly m2, kStrategy strat);
poly preIntegerCheck(ideal F, ideal Q);
void postReduceByMon(LObject* h, kStrategy strat);
void postReduceByMonSig(LObject* h, kStrategy strat);
void finalReduceByMon(kStrategy strat);
// change strat->tailRing and adjust all data in strat, L, and T:
// new tailRing has larger exponent bound
// do nothing and return FALSE if exponent bound increase would result in
// larger exponent bound that that of currRing
BOOLEAN kStratChangeTailRing(kStrategy strat,
                             LObject* L = NULL, TObject* T = NULL,
                             // take this as new_expbound: if 0
                             // new expbound is 2*expbound of tailRing
                             unsigned long new_expbound = 0);
// initiate a change of the tailRing of strat -- should be called
// right before main loop in bba
void kStratInitChangeTailRing(kStrategy strat);

/// Output some debug info about a given strategy
void kDebugPrint(kStrategy strat);

// getting sb order for sba computations
ring sbaRing(kStrategy strat, const ring r=currRing, BOOLEAN complete=TRUE, int sgn=1);

KINLINE void clearS (poly p, unsigned long p_sev, int* at, int* k,
  kStrategy strat);

#include "kernel/GBEngine/kInline.h"

/* shiftgb stuff */
#include "kernel/GBEngine/shiftgb.h"

#ifdef HAVE_SHIFTBBA
static inline BOOLEAN kExistsInL1(const poly p, const kStrategy strat)
{
  for(auto it = strat->L.begin(); it != strat->L.end(); ++it)
  {
    if (p == it->p1) return TRUE;
  }
  return FALSE;
}

void enterTShift(LObject p, kStrategy strat, int atT = -1);

BOOLEAN enterOnePairShift (poly q, poly p, int ecart, int isFromQ, kStrategy strat, int atR, int ecartq, int qisFromQ, int shiftcount, int ifromS);

void enterpairsShift (poly h,int k,int ecart,int pos,kStrategy strat, int atR);

void superenterpairsShift (poly h,int k,int ecart,int pos,kStrategy strat, int atR);

poly redtailBbaShift (LObject* L, int pos, kStrategy strat, BOOLEAN withT, BOOLEAN normalize);

int redFirstShift (LObject* h,kStrategy strat); // ok

ideal bbaShift(ideal F, ideal Q,intvec *w,bigintmat *hilb,kStrategy strat);
#endif

static inline void kDeleteLcm(LObject *P)
{
 if (P->lcm!=NULL)
 {
   if (rField_is_Ring(currRing))
     pLmDelete(P->lcm);
   else
     pLmFree(P->lcm);
   P->lcm=NULL;
 }
}

void initenterpairs (poly h,int k,int ecart,int isFromQ,kStrategy strat, int atR = -1);
#endif
