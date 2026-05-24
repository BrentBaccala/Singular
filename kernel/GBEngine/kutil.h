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
#include <pthread.h>
#include <atomic>
#include <cstdint>

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

// Forward declarations for bench-only elision toggles (task 363; run 587)
// — needed in BlockArray::operator[] and sBasisSet's lock wrappers below.
// Full declarations / rationale near line ~330 (alongside kt_current_ctx).
struct SweepContext;
extern SweepContext* kt_current_ctx;
extern bool g_bench_elide_blockarray_atomic;
extern bool g_bench_elide_sbasis_rwlock;
extern bool g_bench_elide_lset_wrapper;
// Serial compact-on-pop knob (task 368): 0=off 1=threshold 2=everypop.
// Read once at startup from SINGULAR_BENCH_SERIAL_COMPACT; used only
// at serial-mode (kt_current_ctx == NULL) LSet::pop / pop_and_erase.
extern int g_bench_serial_compact;

// Serial L-erase contract knob (task 370): 0=tombstone 1=physical.
// Read once at startup from SINGULAR_BENCH_SERIAL_LSET_ERASE; used
// only at serial-mode (kt_current_ctx == NULL) LSet::erase via the
// filtered_iterator overload (the chainCritNormal chain-crit deletion
// path).  0=tombstone is the current behaviour (true no-op; non-bench
// / default builds entirely unaffected); 1=physical routes those
// serial deletions through LSet::physical_erase (immediate tree-node
// removal + poly-free at deletion time, the next-opt erase contract /
// the pre-task-360 named primitive LSetChunk::physical_erase).
// Parallel mode never takes the physical branch (gated kt_current_ctx
// == NULL), so T>1 is byte- and TSan-identical regardless of value.
extern int g_bench_serial_lset_erase;
enum { KT_SERIAL_LSET_ERASE_TOMBSTONE = 0, KT_SERIAL_LSET_ERASE_PHYSICAL = 1 };

// Serial plain-load scan knob (prototype): when true, the hot
// chainCritNormal chain-criterion scan (filtered_iterator::advance)
// reads the sev_flat_ deletion sentinel with a PLAIN non-atomic load
// instead of __atomic_load_n(ACQUIRE).  The atomic acquire emits the
// same x86 mov as a plain load, BUT it is a compiler optimization
// barrier: GCC will not unroll the scan loop across an atomic load.
// next-opt's advance (plain load) is 4x-unrolled; parallel-bba's
// (atomic) is a scalar non-unrolled loop, which is the bulk of the
// advance 16.5%->23.0% self-time gap.  SERIAL-MODE ONLY: a plain load
// is racy against concurrent worker writers, so this must only be set
// at SINGULAR_THREADS=1.  Default false -> atomic path -> byte- and
// TSan-identical to the shipped behaviour.  Read once at startup from
// SINGULAR_BENCH_SERIAL_PLAIN_SCAN.
extern bool g_bench_serial_plain_scan;

// Block-allocated array: elements are stored in fixed-size blocks
// indexed via a two-level directory.  blocks[b] is a calloc'd block
// of BLOCK_SIZE elements; operator[](i) returns blocks[i>>SHIFT][i&MASK].
// Existing blocks never move when new blocks are appended, so pointers
// to elements (&arr[i]) remain stable across growth.
//
// Concurrency: parallel-bba readers (sweep_one_tile, reduce_slot_from_sweep)
// access strat->T / strat->sevT / strat->R without taking the strat's
// rwlock; only the strat's `count` atomic, the per-T-entry published
// release/acquire flag, and the BlockArray directory pointer's atomicity
// (this class) protect the read path.  Writers (enterT) hold the strat
// wrlock.
template<typename Elem, int BLOCK_SHIFT = 10>
class BlockArray {
  static const int BLOCK_SIZE = 1 << BLOCK_SHIFT;
  static const int BLOCK_MASK = BLOCK_SIZE - 1;
  // Directory pointer.  Atomic because the lockless reader path
  // (operator[], addr()) loads it on every access while ensure_capacity
  // can replace it from a concurrent writer thread.  See ensure_capacity
  // for the full memory-ordering story.
  std::atomic<Elem**> blocks;
  int num_blocks;       // current number of allocated blocks (writer-only)
  int dir_capacity;     // allocated directory slots (writer-only)
  // Old directory pointers retained until free_all().  When the
  // directory grows we cannot free the old one immediately because
  // lockless readers may still be dereferencing it (see ensure_capacity
  // comment).  retained_dirs holds them so the memory stays valid for
  // the lifetime of the BlockArray; free_all() drops them all at
  // teardown.  Per-grow size is geometric (4, 8, 16, ...), so total
  // retained bytes are bounded by ~2× the final directory size — a few
  // KB for realistic strat->T sizes (200K entries → <2 KB).
  std::vector<Elem**> retained_dirs;
protected:
  // Atomic so parallel-bba readers (pop_and_prepare under L_lock,
  // sweep workers without an S-lock) can load the size without
  // racing the enterT writer (which holds S-exclusive).  release on
  // setsize / acquire on size pairs with the existing tobject_publish
  // release-store on the per-entry published flag, giving readers a
  // stable size→entry-content ordering on weakly-ordered hardware.
  std::atomic<int> count;
public:
  BlockArray() : blocks(nullptr), num_blocks(0), dir_capacity(0), count(0) {}

  // Bench toggle (task 363; run 587): when g_bench_elide_blockarray_atomic
  // is set AND we're in serial mode (kt_current_ctx == NULL), use a
  // relaxed load instead of acquire on the directory pointer.  No-op
  // in parallel mode — parallel readers MUST stay acquire to pair
  // with ensure_capacity's release-store.  See kutil.h header comment
  // for the contract.
  //
  // The toggle's branch is __builtin_expect-marked false to keep the
  // common-case fast path branch-predicted correctly: when no toggle
  // is set (production behaviour), the branch falls through and the
  // hot-path operator[] is identical to its pre-task-363 form.
  // Note that on x86-64, acquire-load on a naturally-aligned word is
  // implemented as a plain MOV — no fence, no LOCK prefix — so even
  // without the toggle the cost is just the load itself.
  Elem& operator[](int i) {
    Elem** dir;
    if (__builtin_expect(g_bench_elide_blockarray_atomic && kt_current_ctx == NULL, 0)) {
      dir = blocks.load(std::memory_order_relaxed);
    } else {
      dir = blocks.load(std::memory_order_acquire);
    }
    return dir[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }
  const Elem& operator[](int i) const {
    Elem** dir;
    if (__builtin_expect(g_bench_elide_blockarray_atomic && kt_current_ctx == NULL, 0)) {
      dir = blocks.load(std::memory_order_relaxed);
    } else {
      dir = blocks.load(std::memory_order_acquire);
    }
    return dir[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }

  // Return pointer to element i (stable across growth — see class comment)
  Elem* addr(int i) {
    Elem** dir;
    if (__builtin_expect(g_bench_elide_blockarray_atomic && kt_current_ctx == NULL, 0)) {
      dir = blocks.load(std::memory_order_relaxed);
    } else {
      dir = blocks.load(std::memory_order_acquire);
    }
    return &dir[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }

  // const overload of addr(): same stable-pointer / lockless-reader
  // contract, used by the block-wise SIMD sev scan which holds a
  // `const BlockArray<unsigned long>&`.  Must take the SAME atomic
  // directory-load path as operator[]/addr so the acquire-load still
  // pairs with ensure_capacity's release-store (see ensure_capacity's
  // memory-ordering comment); do NOT cache the directory across calls.
  const Elem* addr(int i) const {
    Elem** dir;
    if (__builtin_expect(g_bench_elide_blockarray_atomic && kt_current_ctx == NULL, 0)) {
      dir = blocks.load(std::memory_order_relaxed);
    } else {
      dir = blocks.load(std::memory_order_acquire);
    }
    return &dir[i >> BLOCK_SHIFT][i & BLOCK_MASK];
  }

  // Block geometry, exposed for the block-wise SIMD sev scan.  Elements
  // are contiguous within a block: indices i .. (i | BLOCK_MASK) of one
  // block share the same inner allocation, so a SIMD kernel may run on
  // that contiguous span.  block_size() is the run length; the in-block
  // remaining run from index i is contiguous_run(i).
  static constexpr int block_size() { return BLOCK_SIZE; }
  static int contiguous_run(int i) { return BLOCK_SIZE - (i & BLOCK_MASK); }

  // Number of elements in use
  int size() const { return count.load(std::memory_order_acquire); }
  bool empty() const { return size() == 0; }

  // Set the count directly (for migration from external tl/sl counters)
  void setsize(int n) { count.store(n, std::memory_order_release); }

  // Ensure at least n elements are allocated (indices 0..n-1).
  //
  // Lock discipline.  Caller must hold the BlockArray's wrlock (in
  // parallel-bba: skStrategy's S-exclusive rwlock, taken by enterT
  // and friends).  Concurrent readers may dereference the directory
  // through operator[] / addr() without holding any lock; the protocol
  // below makes that safe even across a directory regrow.
  //
  // Two-level layout.  `blocks` is a directory of length `dir_capacity`
  // pointing at calloc'd blocks of BLOCK_SIZE elements each.  Element i
  // lives at blocks[i >> BLOCK_SHIFT][i & BLOCK_MASK].  Per-block
  // pointers (the inner allocations) never move once allocated:
  // ensure_capacity only ever fills new directory slots (b in
  // [num_blocks, needed_blocks)) with fresh calloc'd blocks, never
  // touches existing entries 0..num_blocks-1.  Hence Elem* pointers
  // returned by addr(i) stay valid forever.
  //
  // Directory growth is the load-bearing piece for thread safety.
  // When needed_blocks exceeds dir_capacity we calloc a larger
  // directory, memcpy the existing block pointers into it, and
  // publish it via a release-store on `blocks`.  Concurrent readers
  // either see the old directory (whose contents are unchanged for
  // their indices — same per-block pointers) or the new one
  // (everything they need is at the same indices).  Either way the
  // dereference of dir[i >> BLOCK_SHIFT] is safe.
  //
  // We MUST NOT free the old directory here.  A reader may have
  // already loaded the old `blocks` pointer and be about to
  // dereference dir[i >> BLOCK_SHIFT] when the writer reaches this
  // line; calling free() would yank the directory out from under
  // it.  Instead the old directory is retained in retained_dirs and
  // freed at free_all() time, when the BlockArray is single-threaded.
  // Per-grow leak is geometric, so retained_dirs bytes are bounded
  // by the largest directory size — negligible.
  //
  // Memory ordering.  The release-store on `blocks` synchronises
  // with the acquire-load in operator[] / addr(): a reader that sees
  // the new directory pointer is guaranteed to see the memcpy'd
  // contents (the old block pointers in the new directory) and the
  // newly-allocated block pointers in slots [num_blocks, needed_blocks).
  // Readers using indices below the writer's current `count` will
  // pair with `count`'s release/acquire — they only ever load `blocks`
  // after observing a `count` value that already implies the directory
  // was sized large enough.
  //
  // The relaxed-load of `blocks` at the top is fine because we hold
  // the wrlock — no other writer can race us — and we only need the
  // value for our own write.
  void ensure_capacity(int n) {
    int needed_blocks = (n + BLOCK_SIZE - 1) >> BLOCK_SHIFT;
    if (needed_blocks <= num_blocks) return;
    Elem** cur = blocks.load(std::memory_order_relaxed);
    if (needed_blocks > dir_capacity) {
      int new_cap = dir_capacity == 0 ? 4 : dir_capacity;
      while (new_cap < needed_blocks) new_cap *= 2;
      Elem **new_dir = (Elem **)calloc(new_cap, sizeof(Elem *));
      if (cur != NULL) {
        memcpy(new_dir, cur, num_blocks * sizeof(Elem *));
        retained_dirs.push_back(cur);  // can't free — readers may still hold it
      }
      blocks.store(new_dir, std::memory_order_release);
      cur = new_dir;
      dir_capacity = new_cap;
    }
    // Allocate new blocks (zero-initialized).  Existing entries
    // cur[0..num_blocks-1] are untouched, so concurrent readers below
    // the published count see no change.
    for (int b = num_blocks; b < needed_blocks; b++) {
      cur[b] = (Elem *)calloc(BLOCK_SIZE, sizeof(Elem));
    }
    num_blocks = needed_blocks;
  }

  // Return current capacity (total elements allocated)
  int capacity() const { return num_blocks * BLOCK_SIZE; }

  // Insert val at position pos, shifting existing elements up.
  // Automatically grows capacity and increments count.  Caller is
  // expected to hold an exclusive lock on the array (BlockArray's
  // own mutators are not safe against concurrent readers regardless
  // of the count atomicity — the element-shift loop is non-atomic).
  void insert(int pos, const Elem& val) {
    int c = count.load(std::memory_order_relaxed);
    ensure_capacity(c + 1);
    for (int i = c; i > pos; i--)
      (*this)[i] = (*this)[i-1];
    (*this)[pos] = val;
    count.store(c + 1, std::memory_order_release);
  }

  // Append val at the end. Grows capacity and increments count.
  void push_back(const Elem& val) {
    int c = count.load(std::memory_order_relaxed);
    ensure_capacity(c + 1);
    (*this)[c] = val;
    count.store(c + 1, std::memory_order_release);
  }

  // Erase element at position pos, shifting elements down. Decrements count.
  void erase(int pos) {
    int c = count.load(std::memory_order_relaxed);
    for (int i = pos; i < c - 1; i++)
      (*this)[i] = (*this)[i+1];
    count.store(c - 1, std::memory_order_release);
  }

  // Free all blocks, the current directory, and any retained old
  // directories from prior growths.  Resets count.  Single-threaded
  // by construction (called from skStrategy teardown after workers
  // are joined).
  void free_all() {
    Elem** cur = blocks.load(std::memory_order_relaxed);
    for (int b = 0; b < num_blocks; b++) {
      free(cur[b]);
    }
    if (cur != NULL) free(cur);
    for (Elem** old : retained_dirs) free(old);
    retained_dirs.clear();
    blocks.store(nullptr, std::memory_order_relaxed);
    num_blocks = 0;
    dir_capacity = 0;
    count.store(0, std::memory_order_relaxed);
  }
};

typedef class sTObject TObject;
typedef class sLObject LObject;

typedef struct denominator_list_s denominator_list_s;
typedef denominator_list_s *denominator_list;

struct denominator_list_s{number n; denominator_list next;};
EXTERN_VAR denominator_list DENOMINATOR_LIST;

// SElement: one element of the S set (standard basis).
//
// The `deleted` field is a plain bool but is accessed via atomic ops
// from the parallel phase-1 drain path (task 299; run 506).  Writers in phase 1
// hold only a shared lock on the enclosing sBasisSet, so concurrent
// writers must CAS to claim "tombstone ownership" (the winning thread
// is responsible for any owner-side cleanup such as freeing polys).
// Readers under shared lock use a relaxed load; a stale `false` read
// is benign (the walk just considers the entry live a little longer).
//
// Helper functions provide the atomic ops while keeping the struct
// copyable (std::atomic<bool> is not copy-assignable, which would
// break existing code that assigns SElements by value during
// reorder/compact/shift).  See deleted_atomic_load / cas helpers
// below the struct.
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
  // arrival_id: monotonic counter captured at enterS() time, used by
  // parallel phase-1 drainers to iterate "entries that arrived before
  // me" — enterS places at a sorted position, so insertion index does
  // not reflect arrival order under SORDER_STANDARD.  Serial code also
  // bumps the counter to keep arrival_id monotonic over runs.  0 is a
  // valid id for the very first S entry; set to UINT64_MAX in the
  // default constructor so that pre-parallel-drain entries (entered
  // before the counter was started) trivially satisfy "arrived before
  // me" checks regardless of my_arrival.
  uint64_t arrival_id;
  bool deleted;              // lazy-delete flag for parallel phase (atomic access)
  mutable bool pairtest;     // transient: true if spoly(this, h)==0 during enterOnePair

  SElement() : p(NULL), ecart(0), sev(0), s_2_r(0), length(0), wlength(0),
               fromQ(0), sig(NULL), sevSig(0), arrival_id(0),
               deleted(false), pairtest(false) {}
};

// Thread-local override for strat->B (task 301; run 508).
// When non-NULL, the parallel phase-1 drain sets this before calling into
// enterpairs / initenterpairs; the enterOnePair family and chainCrit family
// write into *t_local_B_override instead of strat->B.  Serial code leaves
// this NULL throughout, so the helper strat_B(strat) is a single-pointer
// compare + branch and falls through to strat->B with no behavioural
// change.  Using a thread-local pointer (rather than threading an LSetChunk*
// parameter through ~10 call signatures) keeps the diff surgical while
// giving phase-1 drainers their own private B.
//
// Safety: the pointed-to LSetChunk must outlive every call into enterpairs made
// while the pointer is set.  The drain allocates the LSetChunk on its own stack
// in process_survivor_lobject and restores the pointer (to NULL or the
// previous value) on exit.
class LSetChunk;  // forward declare (defined later in this header)
extern __thread LSetChunk* t_local_B_override;

// kt_current_ctx is forward-declared near the top of this header for
// early use in BlockArray's bench-elision gates.  Defined in kthread.cc.
// Used by inline LSet::erase wrappers to dispatch between tombstone-only
// erase (when running inside bba_parallel_loop, where worker rdlock
// holders would UAF on physical removal) and physical erase + poly-free
// (when running in serial mode).  bba_parallel_loop sets kt_current_ctx
// at entry (kthread.cc:2070) and clears it at exit (kthread.cc:2373);
// the dispatch is one atomic load + branch per LSet::erase call.

// --- Bench-only elision toggles (task 363; run 587) ---
//
// Three runtime knobs to attribute the parallel-bba T=1 structural gap
// (post-f392b3c33) among three suspected serial-mode overheads:
//
//   SINGULAR_BENCH_ELIDE_SBASIS_RWLOCK=1   skip rwlock/sBasisSet plumbing
//   SINGULAR_BENCH_ELIDE_BLOCKARRAY_ATOMIC=1 use relaxed load on
//                                          BlockArray<Elem>::blocks
//                                          (read on every T access)
//   SINGULAR_BENCH_ELIDE_LSET_WRAPPER=1    short-circuit chunked LSet
//                                          to first/only chunk in serial
//                                          mode (skip chunk-walk indirection)
//
// All three are SERIAL-MODE ONLY: each elision site is gated on
// (g_bench_elide_X && kt_current_ctx == NULL).  Parallel-mode behaviour
// is unchanged.  Toggles are read ONCE at process startup (from
// kt_bench_toggles_init() in kutil.cc, invoked at the top of siInit)
// into module-static booleans — getenv on every call would itself
// perturb the measurement.
//
// **Measurement scaffolding, not permanent feature**.  After the
// attribution-report task lands, the user decides whether to keep
// (convert to proper kt_current_ctx dispatch like f392b3c33), revert,
// or leave behind a documented dev knob.
// g_bench_elide_* are forward-declared near the top of this header for
// early use in BlockArray's gates and sBasisSet's lock wrappers.
void kt_bench_toggles_init();

// strat_B(strat) is the LSetChunk enterOnePair/chainCrit should write into.
// Returns *t_local_B_override if set, else strat->B.  Defined after
// skStrategy as an inline function (needs both LSetChunk and skStrategy
// complete).

// Thread-local "my arrival" for the parallel phase-1 drain
// (task 301; run 508).
//
// When set to a value != UINT64_MAX, the enterpairs / initenterpairs /
// chainCrit family filter their S-iteration to entries with
// arrival_id < t_local_my_arrival.  This realises the "iterate entries
// that arrived before h" semantics needed under SORDER_STANDARD —
// because enterS places h at a sorted position, the insertion index
// does not reflect arrival order.
//
// Default UINT64_MAX = no filter (serial path; arrival_id < UINT64_MAX
// is trivially true for every finite arrival_id ever stamped).
extern __thread uint64_t t_local_my_arrival;

// Thread-local pairtest-hit list (task 301; run 508).
//
// In serial mode SElement.pairtest is used: enterOnePair sets it for
// S[i] whose spoly with h is zero, chainCritNormal scans S for pairtest
// entries, clear_pairtest clears them all.  This pattern races across
// concurrent drainers: drainer A setting pairtest[i] (because
// spoly(S[i], h_A) == 0) could be observed by drainer B's
// chainCritNormal and applied to B's (S[*], h_B) pairs — but the
// theoretical justification only holds for h_A.
//
// To avoid the race under parallel phase 1, drainers push the
// SElement* of every pairtest-hit into a thread-local vector instead
// of setting SElement.pairtest.  chainCritNormal iterates this vector
// when t_local_my_arrival != UINT64_MAX; the SElement.pairtest field
// is left alone.  The vector is cleared at the start of each phase 1
// and again at the end (defence in depth).
//
// Serial callers (t_local_my_arrival == UINT64_MAX) keep the existing
// SElement.pairtest behaviour unchanged.
extern __thread std::vector<SElement*>* t_local_pairtest_hits;

// Atomic helpers for SElement.deleted.  Serial code can still read/write
// the field directly; parallel phase-1 drainers must use these.
static inline bool selement_deleted_load(const SElement &e) {
  return __atomic_load_n(&e.deleted, __ATOMIC_RELAXED);
}
// CAS: attempt to flip deleted from false->true.  Returns true if the
// caller won the race (is "owner" of cleanup).  Returns false if another
// thread already tombstoned this element.
static inline bool selement_deleted_cas(SElement &e) {
  bool expected = false;
  return __atomic_compare_exchange_n(&e.deleted, &expected, true,
                                     /*weak=*/false,
                                     __ATOMIC_ACQ_REL, __ATOMIC_RELAXED);
}
static inline void selement_deleted_clear(SElement &e) {
  __atomic_store_n(&e.deleted, false, __ATOMIC_RELAXED);
}

// Atomic helpers for SElement.pairtest.  Multiple parallel phase-1
// drainers may concurrently set pairtest=true on the same SElement
// (from enterOnePair when spoly is zero), and concurrently read it
// inside chainCritNormal, and clear_pairtest clears them all.  All
// three operations use atomic RELAXED: writes are either idempotent
// (both threads write true) or mass-clear (clear_pairtest) and a
// stale read just means a pair survives to the next scan.
static inline bool selement_pairtest_load(const SElement &e) {
  return __atomic_load_n(&e.pairtest, __ATOMIC_RELAXED);
}
static inline void selement_pairtest_set(const SElement &e) {
  // pairtest is mutable; cast away const as the C API expects non-const.
  __atomic_store_n(const_cast<bool*>(&e.pairtest), true, __ATOMIC_RELAXED);
}
static inline void selement_pairtest_clear(SElement &e) {
  __atomic_store_n(&e.pairtest, false, __ATOMIC_RELAXED);
}

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
 *                     erase() sets deleted flag (no shifting).
 *   SORDER_MONFIRST:  monomials first, then by degree (posInSMonFirst).
 *                     Used by sba() over rings.
 *                     erase() shifts elements down.
 *   SORDER_APPEND:    append at end, no sorting (parallel phase).
 *                     erase() sets deleted flag (no shifting).
 *
 * Tombstone-on-erase (SORDER_STANDARD, SORDER_APPEND): physical `count`
 * of the underlying BlockArray is not decremented by erase(); only the
 * `deleted` flag on SElement is flipped and live_count_ is decremented.
 * Iterators automatically skip deleted entries via skip_deleted_forward /
 * skip_deleted_backward. Binary search in find_pos still works because
 * sort order is preserved (polys on tombstoned entries are not freed;
 * the enclosing T set owns them). Compaction via compact() packs live
 * entries contiguously; it is called implicitly at the start of
 * reorder() and cleanT().
 *
 * --- Iterator invalidation contract ------------------------------------
 *
 * Iterators are position-based handles {sBasisSet*, int pos_}. Because
 * the underlying BlockArray is append-only for storage blocks and never
 * reallocates existing blocks, iterators remain valid across BlockArray
 * directory reallocation — their held integer index stays correct in
 * the absence of set-level mutations listed below.
 *
 * Mutation operations have the following effects on outstanding iterators:
 *
 * 1. Append (enter_bba / push_back at end / insert() when order_ ==
 *    SORDER_APPEND): does NOT invalidate any existing iterator. An
 *    iterator previously captured at end() now validly refers to the
 *    first appended element. This intentionally diverges from STL
 *    vector::end() semantics, and is relied upon by the parallel
 *    survivor-drain path (e.g. sLObject::checked, to land in Stage B).
 *
 * 2. Insert-in-middle (enter_bba / insert() in sorted mode at a
 *    find_pos-computed position): invalidates all iterators at
 *    positions >= the insert position. Iterators at positions < insert
 *    remain valid.
 *
 * 3. Non-lazy erase (erase_and_next in SORDER_MONFIRST, or erase() in
 *    SORDER_MONFIRST mode): shifts elements down. Invalidates all
 *    iterators at positions > the erased position. The iterator AT
 *    the erased position is returned by erase_and_next repointed at
 *    whatever element slid into the vacated slot (or end() if the
 *    erased element was last).
 *
 * 4. Lazy erase (erase() in SORDER_STANDARD or SORDER_APPEND mode):
 *    sets the deleted tombstone flag on the element. Invalidates only
 *    the iterator at the tombstoned position — it becomes a "tombstone
 *    iterator" whose operator++ will skip forward past all consecutive
 *    tombstones and land on the next live element.
 *
 * 5. Reorder (reorder()): invalidates ALL outstanding iterators. Callers
 *    must drop every stored iterator before calling reorder(). reorder()
 *    also compacts tombstones as a side effect.
 *
 * 6. Compact (compact()): invalidates ALL outstanding iterators. Callers
 *    must drop every stored iterator before calling compact().
 *
 * 7. Responsibility is on callers. There is no runtime enforcement of
 *    these rules; violations are undefined behaviour (in practice: a
 *    stale iterator refers to a different logical element than the one
 *    the caller captured).
 */
enum SOrderMode { SORDER_STANDARD, SORDER_MONFIRST, SORDER_APPEND };

class sBasisSet : private BlockArray<SElement> {
public:
  // --- Iterator (skips deleted entries) ---
  class iterator {
    sBasisSet* set_;
    int pos_;
    friend class sBasisSet;

    // Uses atomic load on SElement.deleted so phase-1 readers under
    // shared lock see a consistent value even while peer drainers
    // CAS-tombstone entries (task 299; run 506).
    void skip_deleted_forward() {
      while (pos_ < set_->count.load(std::memory_order_relaxed) && selement_deleted_load(set_->elem(pos_))) pos_++;
    }
    void skip_deleted_backward() {
      while (pos_ > 0 && selement_deleted_load(set_->elem(pos_))) pos_--;
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
      // Caller-gated elision: in serial mode sBasisSet::erase
      // physically removes the element, so no tombstones can exist;
      // skipping over them is dead weight.  Function-level dispatch
      // is insufficient because libSingular.so's default visibility
      // forces self-calls through the PLT (see the
      // TObjectIterator::skip_unpublished_forward analysis in
      // docs/profile-parallel-bba-43277bf8d-c200-1.md).
      if (kt_current_ctx != NULL) skip_deleted_forward();
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator& operator--() {
      pos_--;
      if (kt_current_ctx != NULL) skip_deleted_backward();
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
      while (pos_ < set_->count.load(std::memory_order_relaxed) && selement_deleted_load(set_->elem(pos_))) pos_++;
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
      if (kt_current_ctx != NULL) skip_deleted_forward();
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
  sBasisSet() : order_(SORDER_STANDARD), live_count_(0), pairtest_any_(false),
                deleted_count_(0), peak_deleted_count_(0),
                erase_call_count_(0), compact_call_count_(0) {
    pthread_rwlock_init(&rwlock_, NULL);
  }
  ~sBasisSet() {
    pthread_rwlock_destroy(&rwlock_);
  }
  // Non-copyable / non-movable to keep rwlock identity stable.
  sBasisSet(const sBasisSet&) = delete;
  sBasisSet& operator=(const sBasisSet&) = delete;

  // --- Parallel locking (task 299; run 506) ---
  // Replaced the prior pthread_mutex with a pthread_rwlock to support the
  // phase-0/1/2 drain design: phase-0 writers (enterS append) take
  // lock_exclusive(); phase-1 readers (iteration + tombstone writes via
  // atomic CAS on SElement.deleted) take lock_shared().  The serial code
  // paths (THREADS=1, bba(), sba(), slimgb(), etc.) never acquire this
  // lock; they are single-threaded so there is no contention.
  //
  // The legacy lock()/unlock() names map to exclusive mode so callers that
  // were using the mutex prior to task 299 (run 506) continue to get the same
  // serialisation semantics unchanged.
  // See ~/project/docs/parallel-bba-thread-safety-report.md for rationale.
  //
  // Bench toggle (task 363; run 587): when g_bench_elide_sbasis_rwlock is
  // set AND we're in serial mode (kt_current_ctx == NULL), make the
  // pthread call a no-op.  Note: at T=1 serial mode, no caller in the
  // current tree actually invokes these wrappers — all sBasisSet rwlock
  // acquisition happens in kthread.cc parallel paths.  The toggle is
  // wired up for completeness in case future serial-mode code paths
  // acquire it.
  void lock() {
    if (g_bench_elide_sbasis_rwlock && kt_current_ctx == NULL) return;
    pthread_rwlock_wrlock(&rwlock_);
  }
  void unlock() {
    if (g_bench_elide_sbasis_rwlock && kt_current_ctx == NULL) return;
    pthread_rwlock_unlock(&rwlock_);
  }
  void lock_exclusive() {
    if (g_bench_elide_sbasis_rwlock && kt_current_ctx == NULL) return;
    pthread_rwlock_wrlock(&rwlock_);
  }
  void unlock_exclusive() {
    if (g_bench_elide_sbasis_rwlock && kt_current_ctx == NULL) return;
    pthread_rwlock_unlock(&rwlock_);
  }
  void lock_shared() {
    if (g_bench_elide_sbasis_rwlock && kt_current_ctx == NULL) return;
    pthread_rwlock_rdlock(&rwlock_);
  }
  void unlock_shared() {
    if (g_bench_elide_sbasis_rwlock && kt_current_ctx == NULL) return;
    pthread_rwlock_unlock(&rwlock_);
  }
  pthread_rwlock_t* raw_rwlock() { return &rwlock_; }

  SOrderMode order() const { return order_; }
  void set_order(SOrderMode m) { order_ = m; }
  bool is_append_mode() const { return order_ == SORDER_APPEND; }

  // --- Iterators ---
  iterator begin() {
    iterator it(this, 0);
    if (kt_current_ctx != NULL) it.skip_deleted_forward();
    return it;
  }
  iterator end() { return iterator(this, count.load(std::memory_order_relaxed)); }

  const_iterator begin() const {
    const_iterator it(this, 0);
    if (kt_current_ctx != NULL) it.skip_deleted_forward();
    return it;
  }
  const_iterator end() const { return const_iterator(this, count.load(std::memory_order_relaxed)); }

  // cbegin/cend: explicit const_iterator access even when *this is non-const.
  // Useful for capturing a stable "end at snapshot time" iterator
  // (sLObject::checked) from a mutable strat->S.
  const_iterator cbegin() const {
    const_iterator it(this, 0);
    if (kt_current_ctx != NULL) it.skip_deleted_forward();
    return it;
  }
  const_iterator cend() const { return const_iterator(this, count.load(std::memory_order_relaxed)); }

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
    return iterator(this, count.load(std::memory_order_relaxed) - 1);
  }

  // --- Erase ---
  // SORDER_STANDARD, SORDER_APPEND: set deleted flag (no shifting).
  //   The physical `count` of the underlying BlockArray is unchanged;
  //   live_count_ is decremented. Compaction packs the array later via
  //   compact() (called by reorder(), cleanT(), or explicitly by callers).
  // SORDER_MONFIRST: shift elements down (legacy shift semantics).
  //   MONFIRST is only used by sba() over rings and is uncommon; keeping
  //   shift avoids touching the binary-search assumptions of the
  //   monfirst-count logic in find_pos_monfirst.
  //
  // Parallel (task 299; run 506): phase-1 drainers may call erase() under shared
  // S-lock concurrently.  The CAS on SElement.deleted ensures only one
  // thread succeeds; only the winning thread decrements live_count_ and
  // bumps deleted_count_.  live_count_ / deleted_count_ / erase_call_count_
  // are atomically incremented so serial readers stay consistent with
  // the parallel writers.
  void erase(iterator it) {
    // Serial-mode dispatch (mirrors LSet::erase, commit f392b3c33):
    // outside bba_parallel_loop, no peer drainer is reading S under
    // shared lock, so we can physically remove the element.  This
    // avoids tombstone accumulation that would otherwise force every
    // skip_deleted_forward / skip_deleted_backward to walk past
    // dead entries on each iteration.  Physical erase shifts
    // subsequent elements down (BlockArray<SElement>::erase), so
    // iterator stability after erase requires the erase_and_next
    // pattern (already used by the relevant call sites).
    if (order_ == SORDER_MONFIRST || kt_current_ctx == NULL) {
      __atomic_add_fetch(&erase_call_count_, 1, __ATOMIC_RELAXED);
      BlockArray<SElement>::erase(it.pos_);
      live_count_--;
    } else {
      // CAS-based tombstone: only the winning thread mutates counters.
      if (selement_deleted_cas(elem(it.pos_))) {
        __atomic_add_fetch(&erase_call_count_, 1, __ATOMIC_RELAXED);
        __atomic_add_fetch(&deleted_count_, 1, __ATOMIC_RELAXED);
        __atomic_sub_fetch(&live_count_, 1, __ATOMIC_RELAXED);
      }
    }
  }

  // --- Capacity / storage ---
  using BlockArray<SElement>::ensure_capacity;
  using BlockArray<SElement>::capacity;
  using BlockArray<SElement>::free_all;

  // Set element count directly (used during initialization and by
  // paths that rebuild S from scratch, e.g. SCA's kInterRedOld).
  // Clears tombstone flags for the retained n slots so the caller sees
  // a clean array with live_count == physical count.
  void setsize(int n) {
    BlockArray<SElement>::setsize(n);
    live_count_ = n;
    deleted_count_ = 0;
    for (int i = 0; i < n; i++) elem(i).deleted = false;
  }

  // --- Pairtest ---
  // Clear all pairtest flags and the sentinel.  Parallel phase-1 drainers
  // that concurrently read pairtest see atomically-cleared values; a stale
  // `true` read just means one extra (safe) chain-crit pass.
  void clear_pairtest() {
    int n = count.load(std::memory_order_relaxed);
    for (int i = 0; i < n; i++)
      selement_pairtest_clear(elem(i));
    __atomic_store_n(&pairtest_any_, false, __ATOMIC_RELAXED);
  }
  // Set the sentinel (some zero spoly was found).  Atomic because phase-1
  // drainers may race on this flag from enterOnePair.
  void set_pairtest_any() {
    __atomic_store_n(&pairtest_any_, true, __ATOMIC_RELAXED);
  }
  bool has_pairtest() const {
    return __atomic_load_n(&pairtest_any_, __ATOMIC_RELAXED);
  }

  // Compact: remove entries whose .p field is NULL (typically set by
  // updateResult's pDelete pass). Differs from compact() in that this
  // filters on p==NULL rather than on the deleted tombstone flag; some
  // callers pDelete directly without marking the element deleted.
  // Also drops any tombstoned entries (whose p might still be non-NULL)
  // because after compact_null_p the intent is a packed live array.
  void compact_null_p() {
    int n = count.load(std::memory_order_relaxed);
    int dst = 0;
    for (int src = 0; src < n; src++) {
      if (elem(src).p != NULL && !elem(src).deleted) {
        if (dst != src)
          elem(dst) = elem(src);
        dst++;
      }
    }
    count.store(dst, std::memory_order_relaxed);
    live_count_ = dst;
    deleted_count_ = 0;
  }

  // Compact: remove deleted entries, pack remaining entries contiguously.
  // Walks the physical array in place; the live entries retain their
  // relative ordering (so SORDER_STANDARD sort-order is preserved).
  // Resets physical count to live_count_ and clears deleted_count_.
  // Does NOT change order_ — callers that want SORDER_STANDARD should
  // set it explicitly.
  // Invalidates all outstanding iterators.
  void compact() {
    compact_call_count_++;
    if (deleted_count_ > peak_deleted_count_)
      peak_deleted_count_ = deleted_count_;
    int n = count.load(std::memory_order_relaxed);
    if (deleted_count_ == 0 && live_count_ == n) {
      // Fast path: nothing to do.
      return;
    }
    int dst = 0;
    for (int src = 0; src < n; src++) {
      if (!elem(src).deleted) {
        if (dst != src)
          elem(dst) = elem(src);
        dst++;
      }
    }
    count.store(dst, std::memory_order_relaxed);
    live_count_ = dst;
    deleted_count_ = 0;
  }

  // --- Erase / compact statistics ---
  // Cumulative counters for diagnostics; use debug_print_stats() to dump.
  // Not thread-safe with live writes; read from a quiescent state.
  int deleted_count() const { return deleted_count_; }
  int peak_deleted_count() const { return peak_deleted_count_; }
  int erase_call_count() const { return erase_call_count_; }
  int compact_call_count() const { return compact_call_count_; }
  int physical_count() const { return count.load(std::memory_order_relaxed); }
  void debug_print_stats(const char *tag = NULL) const {
    fprintf(stderr,
            "sBasisSet[%s] live=%d physical=%d deleted=%d peak_deleted=%d erase_calls=%d compact_calls=%d\n",
            tag ? tag : "", live_count_, count.load(std::memory_order_relaxed), deleted_count_,
            peak_deleted_count_, erase_call_count_, compact_call_count_);
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
  // Returns an iterator to the just-inserted element.
  // atS: end() means "compute position via find_pos internally" (was int atS = -1).
  // arrival_id is set from strat->arrival_counter.fetch_add(1) so the
  // inserted SElement's arrival_id is unique among all elements ever
  // entered into S.  Parallel drain captures my_arrival =
  // strat->arrival_counter.load() BEFORE calling enter_bba, knowing it
  // holds S-exclusive so this enter_bba will get exactly my_arrival.
  iterator enter_bba(LObject &p, kStrategy strat, int atR, iterator atS);

  // Insert for signature-based algorithms (replaces enterSSba).
  // Also copies sig and sevSig fields.
  // Returns an iterator to the just-inserted element.
  // atS: end() means "compute position via find_pos internally" (was int atS = -1).
  iterator enter_sba(LObject &p, kStrategy strat, int atR, iterator atS);

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

  // For tgb.cc: binary search ordering by length/wlength first, then
  // by leading monomial. Replaces tgb.cc's static simple_posInS.
  // Implementation in kutil.cc (needs complete kStrategy).
  iterator simple_find_pos(poly p, int len, wlen_type wlen, kStrategy strat);

  // For tgb.cc: shift an element within S. Both take raw int positions
  // because tgb's slimgb path manages S by raw-index arithmetic
  // (simple_posInS + adjacent shifts); these methods are the sBasisSet-
  // internal encapsulation of that raw shift, using the private iterator
  // constructor to avoid exposing an int->iterator adapter at the class
  // API level. old_pos > new_pos for move_forward; old_pos < new_pos for
  // move_backward.
  iterator move_forward(int old_pos, int new_pos) {
    assume(old_pos >= new_pos);
    SElement saved = elem(old_pos);
    for (int k = old_pos; k > new_pos; k--) elem(k) = elem(k - 1);
    elem(new_pos) = saved;
    return iterator(this, new_pos);
  }
  iterator move_backward(int old_pos, int new_pos) {
    assume(old_pos <= new_pos);
    SElement saved = elem(old_pos);
    for (int k = old_pos; k < new_pos; k++) elem(k) = elem(k + 1);
    elem(new_pos) = saved;
    return iterator(this, new_pos);
  }

  // --- Reverse iteration ---
  // Thin wrappers using the forward iterator's operator-- . Note: the
  // forward iterator's skip_deleted_backward stops at pos_ > 0, so if
  // slot 0 is tombstoned a reverse walk will observe it. Intended for
  // callers that don't have tombstones in S (e.g. ring-coefficient
  // post-compaction paths).
  class reverse_iterator {
    iterator it_;  // points one past the current element (like std::reverse_iterator)
  public:
    reverse_iterator() {}
    explicit reverse_iterator(iterator it) : it_(it) {}
    SElement& operator*() const { iterator tmp = it_; --tmp; return *tmp; }
    SElement* operator->() const { iterator tmp = it_; --tmp; return &*tmp; }
    reverse_iterator& operator++() { --it_; return *this; }
    reverse_iterator operator++(int) { reverse_iterator tmp = *this; --it_; return tmp; }
    bool operator==(const reverse_iterator& o) const { return it_ == o.it_; }
    bool operator!=(const reverse_iterator& o) const { return it_ != o.it_; }
    // Access to the base iterator (one past current element).
    iterator base() const { return it_; }
    // Raw index of the element this reverse iterator points at.
    int index() const { iterator tmp = it_; --tmp; return tmp.index(); }
  };
  reverse_iterator rbegin() { return reverse_iterator(iterator(this, count.load(std::memory_order_relaxed))); }
  reverse_iterator rend()   { return reverse_iterator(iterator(this, 0)); }

  // Binary search for sorted insertion position (replaces free posInS).
  // Scans the whole of S. Returns an iterator to the position the
  // element would occupy if inserted (end() if larger than everything).
  iterator find_pos(const poly p, int ecart_p);

  // Bounded form: restricts the search to [0..end_bound). Only used
  // internally by find_divisor_search_bound.
  iterator find_pos(const poly p, int ecart_p, iterator end_bound);

  // Binary search for monfirst insertion position (replaces free posInSMonFirst).
  iterator find_pos_monfirst(const poly p);

  // Encapsulates the kFindDivisibleByInS search-range narrowing.
  // For monomial-ordered (non-Ring, non-component, non-lex) S, returns an
  // upper bound on the index of any element of S whose leading monomial
  // can divide p, restricted to indices [0..max_ind]. Returns max_ind for
  // rings or other cases where narrowing is not safe; the caller must
  // still check divisibility against each element up to the bound.
  iterator find_divisor_search_bound(poly p, iterator max, kStrategy strat);

private:
  SOrderMode order_;
  int live_count_;
  bool pairtest_any_;  // sentinel: true if any SElement.pairtest was set
  pthread_rwlock_t rwlock_;  // parallel drain path serialization (task 299; run 506)

  // Erase / compact instrumentation (cumulative across the lifetime of
  // the sBasisSet). deleted_count_ is the current tombstone population
  // (decremented by compact()); peak_deleted_count_ tracks the largest
  // value deleted_count_ reached before a compact call;
  // erase_call_count_ / compact_call_count_ are monotonically increasing.
  int deleted_count_;
  int peak_deleted_count_;
  int erase_call_count_;
  int compact_call_count_;

  // Internal: insert at a specific position (for reorder, etc.)
  iterator insert_at(int pos, const SElement& val) {
    BlockArray<SElement>::insert(pos, val);
    live_count_++;
    return iterator(this, pos);
  }

  // Internal raw element access (used by iterators and member functions)
  SElement& elem(int i) { return (*static_cast<BlockArray<SElement>*>(this))[i]; }
  const SElement& elem(int i) const { return (*static_cast<const BlockArray<SElement>*>(this))[i]; }

  friend class skStrategy;
};

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

  // Publication flag (task t-iterator-published, 510).
  //
  // In the parallel worker-side drain (future task), enterT will claim
  // a T-slot under an S-exclusive lock, fill all fields (p, t_p, sev,
  // ecart, length, pLength, ...), and then release-store `published`
  // to true.  Readers on other threads iterate T via a skipping
  // iterator that acquire-loads `published` and skips unpublished
  // slots.  This gives writer-filled-fields-visible-to-readers
  // semantics without a reader-side lock.
  //
  // Serial callers (THREADS=1, and the main-thread drain) observe
  // `published` synchronously: enterT publishes the slot before
  // returning, so every T reader sees it as published by the time it
  // matters.  Access is always via tobject_published_load /
  // tobject_publish helpers (matching SElement.deleted pattern; the
  // struct stays copyable for the T[i] = T[i-1] shift in enterT and
  // other copy-assign sites).
  bool published;


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

// Atomic helpers for sTObject.published (task 303; run 510).
//
// The published flag transitions at most once from false to true (in
// enterT, after all other fields of the T-slot are written).  Readers
// use acquire-load so that every field written before the release-store
// is visible once the flag is observed true.  Writers use release-store
// so prior field stores are flushed to other cores before the flag
// flips.  The field is stored as a plain `bool` (not std::atomic<bool>)
// so that sTObject remains copy-assignable — enterT performs
// T[i] = T[i-1] shifts, and BlockArray<TObject>::insert() does the same
// under the hood, and std::atomic<bool> is not copy-assignable.
static inline bool tobject_published_load(const sTObject &t) {
  // Serial-mode elision: outside bba_parallel_loop, kt_current_ctx is
  // NULL (set/cleared in kthread.cc:2070/:2375). enterT publishes
  // synchronously on serial paths, so every reader observes
  // published=true by the time it matters; the byte load + acquire
  // fence is dead weight and adds a 2nd-cache-line miss per probe in
  // hot loops like kFindDivisibleByInT. Short-circuit on the
  // kt_current_ctx test (L1-resident global, well-predicted branch).
  return kt_current_ctx == NULL
      || __atomic_load_n(&t.published, __ATOMIC_ACQUIRE);
}
static inline void tobject_publish(sTObject &t) {
  __atomic_store_n(&t.published, true, __ATOMIC_RELEASE);
}
static inline void tobject_unpublish(sTObject &t) {
  // Used when shifting slots inside enterT (the slot being overwritten
  // must revert to unpublished until the new writer publishes it).
  // Relaxed is sufficient: the caller holds S-exclusive and there is
  // no reader that could observe the stale published=true between the
  // shift and the new publish.
  __atomic_store_n(&t.published, false, __ATOMIC_RELAXED);
}

// Skipping iterator for BlockArray<TObject> (task 303; run 510).
//
// Mirrors the sBasisSet iterator's skip_deleted_forward shape, but T is
// append-only (no tombstones), so the iterator only needs to skip
// unpublished slots.  operator++ advances, then acquire-loads
// published on each slot until it finds a published one or reaches end.
//
// Usage (future worker-side drain):
//   for (auto it = begin_T(strat->T); it != end_T(strat->T); ++it) {
//       TObject& t = *it;
//       ...
//   }
//
// In this task only the infrastructure is added; no readers migrated.
// Serial code continues to use strat->T[i] / strat->T.size() unchanged.
//
// Iterators are stable across BlockArray growth (blocks never move), so
// a saved iterator remains valid across concurrent enterT appends.  An
// iterator captured at end_T() advances naturally once the appended
// slot publishes (matches the sBasisSet cbegin()/cend() contract).
class TObjectIterator {
  BlockArray<TObject>* set_;
  int pos_;

  // Acquire-load on published — pairs with tobject_publish's release
  // store in enterT.  Once a reader observes published=true, all other
  // field writes that preceded the publish are visible.
  //
  // The skip is gated at the call site (constructor / operator++) by
  // a kt_current_ctx == NULL check, so this function is only entered
  // when running under bba_parallel_loop.  Without that gating the
  // PLT-based dispatch to this function dominates serial-mode wall
  // (libSingular.so is built with default visibility, so even
  // self-calls to skip_unpublished_forward go through PLT).
  void skip_unpublished_forward() {
    while (pos_ < set_->size()
           && !tobject_published_load((*set_)[pos_])) pos_++;
  }

public:
  TObjectIterator() : set_(NULL), pos_(0) {}
  TObjectIterator(BlockArray<TObject>* s, int pos) : set_(s), pos_(pos) {
    if (kt_current_ctx != NULL) skip_unpublished_forward();
  }

  TObject& operator*()  const { return (*set_)[pos_]; }
  TObject* operator->() const { return &(*set_)[pos_]; }

  // Raw index — available for interoperability with integer-based APIs
  // (posInT family, R[i_r] lookups, etc.).  Prefer iterator-based
  // access in new code.
  int index() const { return pos_; }

  TObjectIterator& operator++() {
    pos_++;
    if (kt_current_ctx != NULL) skip_unpublished_forward();
    return *this;
  }
  TObjectIterator operator++(int) {
    TObjectIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const TObjectIterator& o) const { return pos_ == o.pos_; }
  bool operator!=(const TObjectIterator& o) const { return pos_ != o.pos_; }
};

class ConstTObjectIterator {
  const BlockArray<TObject>* set_;
  int pos_;

  // See TObjectIterator::skip_unpublished_forward for the call-site
  // gating rationale.
  void skip_unpublished_forward() {
    while (pos_ < set_->size()
           && !tobject_published_load((*set_)[pos_])) pos_++;
  }

public:
  ConstTObjectIterator() : set_(NULL), pos_(0) {}
  ConstTObjectIterator(const BlockArray<TObject>* s, int pos)
      : set_(s), pos_(pos) {
    if (kt_current_ctx != NULL) skip_unpublished_forward();
  }

  const TObject& operator*()  const { return (*set_)[pos_]; }
  const TObject* operator->() const { return &(*set_)[pos_]; }

  int index() const { return pos_; }

  ConstTObjectIterator& operator++() {
    pos_++;
    if (kt_current_ctx != NULL) skip_unpublished_forward();
    return *this;
  }
  ConstTObjectIterator operator++(int) {
    ConstTObjectIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const ConstTObjectIterator& o) const { return pos_ == o.pos_; }
  bool operator!=(const ConstTObjectIterator& o) const { return pos_ != o.pos_; }
};

// Free-function helpers: match begin()/end()/iterator_at(i) style of
// sBasisSet.  Taking BlockArray<TObject>& (rather than adding methods
// to the BlockArray template) keeps the TObject-specific iteration
// policy out of the generic BlockArray used for sevT
// (BlockArray<unsigned long>) and R (BlockArray<TObject*>), which
// don't have a published flag.
static inline TObjectIterator begin_T(BlockArray<TObject>& T) {
  return TObjectIterator(&T, 0);
}
static inline TObjectIterator end_T(BlockArray<TObject>& T) {
  return TObjectIterator(&T, T.size());
}
static inline TObjectIterator iterator_at_T(BlockArray<TObject>& T, int i) {
  return TObjectIterator(&T, i);
}
static inline ConstTObjectIterator begin_T(const BlockArray<TObject>& T) {
  return ConstTObjectIterator(&T, 0);
}
static inline ConstTObjectIterator end_T(const BlockArray<TObject>& T) {
  return ConstTObjectIterator(&T, T.size());
}
static inline ConstTObjectIterator iterator_at_T(const BlockArray<TObject>& T, int i) {
  return ConstTObjectIterator(&T, i);
}

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
  unsigned seq;       // the sequence number of the LSetChunk when this LObject was inserted
                      // used to determine LSetChunk ordering for equal LObjects
  size_t flat_index;      // index in writable_set's flat array for unordered iteration
  // Iterator pointing "one past the last S element against which this
  // pair was already checked" at creation time.  When the pair is later
  // pulled from L for reduction, rewCrit2 resumes scanning from this
  // iterator forward, skipping S elements already considered.
  //
  // Storing an iterator (rather than a raw unsigned S-size snapshot)
  // makes this value robust against S mutations that preserve prior
  // positions (appends at the end, lazy-mode tombstone erase).  In
  // particular, an iterator captured at S.end() remains valid across
  // later enter_bba() appends — it sticks to the captured position and
  // subsequent appends become the "new" range the resume-scan walks.
  // See sBasisSet's iterator invalidation contract above.
  //
  // Correctness depends on the lazy-mode append-only invariant: S
  // slots are never reused for new logical elements, only tombstoned.
  // sLObject is memset-zeroed in Init(); the resulting null iterator
  // (set_=NULL, pos_=0) is safe as long as it is never dereferenced
  // before being assigned a real value (non-SBA paths never read it).
  sBasisSet::const_iterator checked;
  BOOLEAN prod_crit;
                      // NOTE: If prod_crit = TRUE then the corresponding pair is
                      // detected by Buchberger's Product Criterion and can be
                      // deleted
  int fromQ = 0;      // from quotient ideal; copied into SElement.fromQ by
                      // sBasisSet::enter_bba / enter_sba.
  bool deleted = false; // Lazy-erase tombstone flag (LSetChunk tombstone-on-erase).
                      // Set by LSetChunk::erase(); LSetChunk iterators (ordered,
                      // unordered, filtered) skip entries with deleted==true.
                      // Cleanup of the LObject's polys (lcm/sig/p), pair_index
                      // and sev_flat_ sentinels all happen at erase time;
                      // compact() removes the tombstoned entry from the
                      // multiset tree and deallocates the LObject.
                      // Initialized to false by sLObject::Init() via memset.
                      //
                      // Parallel phase-1 (task 299; run 506): accessed via atomic ops.
                      // Writers under shared S-lock CAS deleted false->true
                      // to claim cleanup ownership; readers use a relaxed
                      // load.  Helper functions below (lobject_deleted_load,
                      // lobject_deleted_cas) perform these ops while keeping
                      // the struct copyable.

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

// Atomic helpers for LObject.deleted.  Serial code can read/write the
// field directly; parallel phase-1 drainers must use these.
static inline bool lobject_deleted_load(const LObject &o) {
  return __atomic_load_n(&o.deleted, __ATOMIC_RELAXED);
}
// CAS: attempt to flip deleted from false->true.  Returns true if the
// caller won the race (is "owner" of cleanup).  Returns false if another
// thread already tombstoned this element.
static inline bool lobject_deleted_cas(LObject &o) {
  bool expected = false;
  return __atomic_compare_exchange_n(&o.deleted, &expected, true,
                                     /*weak=*/false,
                                     __ATOMIC_ACQ_REL, __ATOMIC_RELAXED);
}
static inline void lobject_deleted_clear(LObject &o) {
  __atomic_store_n(&o.deleted, false, __ATOMIC_RELAXED);
}

// Specialize writable_set's tombstone predicate for LObject.  writable_set
// iterators (ordered + unordered) check this to skip tombstoned entries.
// Inlined so that the hot iteration loops remain cheap.  Uses atomic load
// so phase-1 readers under shared lock see a consistent value even while
// peer drainers CAS-tombstone entries.
template<>
inline bool writable_set_is_deleted<LObject>(const LObject& o) {
  return __atomic_load_n(&o.deleted, __ATOMIC_RELAXED);
}

EXTERN_VAR int HCord;

/** @class LSetChunk
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

// Hash function for (poly, poly) pairs used in LSetChunk pair index
struct PolyPairHash {
  std::size_t operator()(const std::pair<poly, poly>& p) const {
    // Combine the two pointer hashes
    std::size_t h1 = std::hash<poly>{}(p.first);
    std::size_t h2 = std::hash<poly>{}(p.second);
    // Simple combination: XOR and shift
    return h1 ^ (h2 << 1);
  }
};

class LSetChunk : public writable_set<LObject, CompareLObject> {
public:
  // Intrusive linked-list link, used by the LSet wrapper to chain
  // chunks together (task 360, step 6).  `next` is an atomic pointer
  // because workers under the wrapper's read-lock CAS-append new
  // chunks concurrently; readers walk the chain via
  // next.load(memory_order_acquire) and observe a published chunk's
  // contents via release/acquire pairing on next.
  //
  // Lifetime: chunks are heap-allocated by the wrapper and freed
  // either at compact time (under the wrapper's writer lock) or at
  // wrapper-destructor time.  The wrapper owns the head chunk by
  // value (chunk_); successor chunks live on the heap and are linked
  // through this pointer.
  std::atomic<LSetChunk*> next{nullptr};

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

  // --- sev_flat_ as the SOLE deletion record for the hot scan (task 379) ---
  //
  // Task 379 (A): sev_flat_ == 0 is the single deletion signal read by
  // the hot chainCritNormal chain-criterion scan (filtered_iterator::
  // advance).  Previously advance() additionally chased flat_ptr(pos_)
  // and did an atomic byte load of LObject.deleted per element; that
  // atomic load was the measured +1.0 s of the T=1 next-opt gap
  // (~/project/reports/parallel-bba-serial-lset-physical-erase-profile.md
  // F3/F6).  To trust sev_flat_ alone, its zeroing stores are promoted
  // to RELEASE and the scan load to ACQUIRE:
  //   * Deletion / pop / physical-erase / compact zero the slot with a
  //     release store (sev_flat_store_zero_rel_), publishing the prior
  //     pair_index / tree mutation to the acquiring scanner.
  //   * The scan reads each slot with an acquire load
  //     (sev_flat_is_deleted_acq_).  On x86 an aligned acquire load is
  //     the same instruction as a plain load — the hot loop pays nothing.
  //   * sevSig_flat_ is NOT a deleted-detector (it can legitimately be
  //     0); its store stays plain — it is only the filter array.  del is
  //     ALWAYS sev_flat_ even when the scan filters on sevSig_flat_.
  //
  // sev_flat_ == 0 covers tombstones (erased), initial generators
  // (sev_lcm == 0, processed by the drain but skipped by the scan),
  // popped-but-not-compacted slots, AND physical holes (flat_[i] ==
  // data_.end(), produced only by erase_at in physical_erase/compact,
  // both of which zero sev_flat_ under wrlock/serial — see physical-hole
  // analysis in the task-379 commit message).  The drain still
  // distinguishes generator (deleted==false → process) from tombstone
  // (deleted==true → skip) via LObject.deleted, which this task leaves
  // in place as the cold drain's record (deferred to (B)).
  static inline bool sev_flat_is_deleted_acq_(const unsigned long* p) {
    return __atomic_load_n(p, __ATOMIC_ACQUIRE) == 0;
  }
  // Erase arbiter: atomically swap the slot to 0 (acq_rel) and report
  // the prior value.  Exactly one thread sees a nonzero prior value
  // (the winner); concurrent peers attempting the same erase see 0.
  unsigned long sev_flat_xchg_zero_(size_t i) {
    return __atomic_exchange_n(&sev_flat_[i], 0UL, __ATOMIC_ACQ_REL);
  }
  // Non-arbiter zeroing (erase, pop, physical_erase, compact precursor):
  // release store so a concurrent acquire reader that observes the 0
  // also observes prior cleanup writes.
  void sev_flat_store_zero_rel_(size_t i) {
    __atomic_store_n(&sev_flat_[i], 0UL, __ATOMIC_RELEASE);
  }

  // --- Tombstone / erase instrumentation ---
  // live_count_: number of non-tombstoned entries (i.e. size()).
  // deleted_count_: number of currently-tombstoned entries (decremented on
  // compact()).  peak_deleted_count_: high-water mark of deleted_count_ seen
  // before a compact.  erase_/compact_/pop_skip_ counters are cumulative.
  int live_count_ = 0;
  int deleted_count_ = 0;
  int peak_deleted_count_ = 0;
  long erase_call_count_ = 0;
  long compact_call_count_ = 0;
  long pop_skip_count_ = 0;  // cumulative tombstones skipped by pop()

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
  using writable_set<LObject, CompareLObject>::size_type;
  using writable_set<LObject, CompareLObject>::unordered_iterator;
  using writable_set<LObject, CompareLObject>::ubegin;
  using writable_set<LObject, CompareLObject>::uend;

  // size()/empty() return LIVE-element counts (excluding tombstones), which
  // is what all downstream GB-engine callers expect.  Physical storage
  // (including tombstones) is exposed via physical_size() / sev_flat_size().
  size_type size() const { return static_cast<size_type>(live_count_); }
  bool empty() const { return live_count_ == 0; }

  // --- Tombstone statistics / diagnostics (not thread-safe) ---
  int deleted_count() const { return deleted_count_; }
  int peak_deleted_count() const { return peak_deleted_count_; }
  long erase_call_count() const { return erase_call_count_; }
  long compact_call_count() const { return compact_call_count_; }
  long pop_skip_count() const { return pop_skip_count_; }
  size_type physical_size() const {
    return writable_set<LObject, CompareLObject>::size();
  }
  void debug_print_stats(const char *tag = NULL) const {
    fprintf(stderr,
            "LSetChunk[%s] live=%d physical=%zu deleted=%d peak_deleted=%d "
            "erase_calls=%ld compact_calls=%ld pop_skips=%ld\n",
            tag ? tag : "", live_count_, physical_size(),
            deleted_count_, peak_deleted_count_,
            erase_call_count_, compact_call_count_, pop_skip_count_);
  }

  // compact(): remove tombstoned entries from the multiset tree, free
  // polys/pair_index entries that were left in-place by lazy erase, and
  // repack sev_flat_ / sevSig_flat_ / flat_index.  Called at cleanup
  // points (completeReduce, end of bba, cleanL).  Invalidates ALL
  // outstanding iterators.
  void compact();

  // Filtered unordered iterator: scans a contiguous sev array
  // for cache-friendly pre-filtering, only visiting elements whose sev
  // passes the filter.  Two filter modes:
  //   sev2_=0: divisibility — skip where (sev1_ & ~sev_array[i]) != 0
  //   sev2_!=0: incomparability — skip where both
  //     (sev_array[i] & ~sev1_) and (sev2_ & ~sev_array[i]) are nonzero
  // The same iterator type is used for both sev_flat_ (lcm) and
  // sevSig_flat_ (signature) scans — only the array pointer differs.
  class filtered_iterator {
    LSetChunk* owner_;
    size_t pos_;
    unsigned long sev1_;
    unsigned long sev2_;
    const std::vector<unsigned long>* sev_array_;
    friend class LSetChunk;

    void advance() {
      // Always use sev_flat_ to detect deleted entries (sev_lcm sentinel 0).
      // Use sev_array_ (which may be sev_flat_ or sevSig_flat_) for the
      // actual filter check. This avoids the problem that sevSig can
      // legitimately be 0, which would be confused with the deleted sentinel.
      //
      // Correctness invariant (task 379 (A)): the chainCritNormal scan
      // trusts sev_flat_ == 0 ALONE — no per-element flat_ptr() chase or
      // LObject.deleted atomic byte load.  Sound iff:
      //   (i) every deletion path zeroes sev_flat_ with RELEASE before
      //       the entry could be re-scanned (erase / pop / physical_erase
      //       / compact all use sev_flat_store_zero_rel_); and
      //   (ii) a live entry never transiently reads sev_flat_ == 0 (a new
      //        entry's live sev_lcm is published by insert() before the
      //        entry becomes scan-visible — the wrlock release that makes
      //        it visible release/acquire-syncs with the scanner's rdlock).
      // sev_flat_ == 0 also catches physical holes (flat_[i] ==
      // data_.end()): the only hole-producing paths are erase_at, used
      // exclusively by physical_erase/compact, both of which zero
      // sev_flat_ under wrlock/serial (exclusive of this rdlock scanner).
      // A stale-LIVE read (scanner sees old nonzero sev_flat_ for a
      // just-erased entry) is harmless: the entry's monomials are not
      // freed until compact under wrlock, which excludes the rdlock
      // scanner, so the scanner at worst considers a soon-dead pair,
      // never a freed one.
      const unsigned long* del = owner_->sev_flat_.data();
      const unsigned long* sev = sev_array_->data();
      const size_t sz = sev_array_->size();
      // Prototype knob: serial-mode plain-load scan.  The dispatch is
      // ONCE per advance() call (outside the inner per-element loop), so
      // the plain-load while-loop below has no atomic in its body and
      // the compiler is free to unroll it (matching next-opt's advance).
      if (g_bench_serial_plain_scan) {
        while (pos_ < sz) {
          if (del[pos_] == 0) { ++pos_; continue; }   // sole deletion signal (plain load)
          unsigned long s = sev[pos_];
          if (sev2_ == 0) {
            if (sev1_ & ~s) { ++pos_; continue; }       // divisibility: skip
          } else {
            if ((s & ~sev1_) && (sev2_ & ~s)) { ++pos_; continue; }  // incomp: skip
          }
          break;
        }
        return;
      }
      while (pos_ < sz) {
        if (sev_flat_is_deleted_acq_(&del[pos_])) { ++pos_; continue; } // sole deletion signal
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
    filtered_iterator(LSetChunk* owner, size_t pos, unsigned long sev1, unsigned long sev2,
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
  // physical_erase via filtered_iterator (task 370): the pre-task-360
  // erase contract reached through the cache-friendly sev scan.
  // Frees the entry's polys (kLSet_free_polys) AND physically removes
  // its multiset tree node + pair_index / sev_flat_ / sevSig_flat_
  // slot — no tombstone, no deferred compact.  Returns the next valid
  // filtered position (same iteration contract as erase()).  Same
  // safety contract as physical_erase(iterator): caller must hold
  // wrlock or be single-threaded (it asserts kt_current_ctx == NULL).
  filtered_iterator physical_erase(filtered_iterator it);

  // Direct access to sev_flat_ for index-based pair iteration (B dedup)
  const unsigned long* sev_flat_data() const { return sev_flat_.data(); }
  size_t sev_flat_size() const { return sev_flat_.size(); }

  // Default constructor
  LSetChunk() = default;

  // Copy constructor: base class copy rebuilds flat_ with new indices,
  // so we must rebuild sev_flat_, sevSig_flat_, and pair_index to match.
  // The copy ctor uses insert() for each element, so no tombstones are
  // carried over; live_count_ will match the iterated count.
  // `next` is NOT copied — copies are independent chunks; the wrapper
  // re-establishes the chain via append.
  LSetChunk(const LSetChunk& other)
    : writable_set<LObject, CompareLObject>(other), seq(other.seq),
      live_count_(0), deleted_count_(0), peak_deleted_count_(0),
      erase_call_count_(0), compact_call_count_(0), pop_skip_count_(0) {
    rebuild_sev_flat();
    rebuild_sevSig_flat();
    rebuild_pair_index();
    // Base class ::insert copies only live entries on construction; count.
    live_count_ = static_cast<int>(physical_size());
    next.store(nullptr, std::memory_order_relaxed);
  }

  // Move constructor: base class move rebuilds flat_ with new indices.
  // `next` is move-loaded from the source so that move-construction
  // preserves the chain link if the source had one.  The source's
  // `next` is reset to nullptr.
  LSetChunk(LSetChunk&& other) noexcept
    : writable_set<LObject, CompareLObject>(std::move(other)), seq(other.seq),
      live_count_(other.live_count_), deleted_count_(other.deleted_count_),
      peak_deleted_count_(other.peak_deleted_count_),
      erase_call_count_(other.erase_call_count_),
      compact_call_count_(other.compact_call_count_),
      pop_skip_count_(other.pop_skip_count_) {
    rebuild_sev_flat();
    rebuild_sevSig_flat();
    rebuild_pair_index();
    other.live_count_ = 0;
    other.deleted_count_ = 0;
    LSetChunk* n = other.next.load(std::memory_order_relaxed);
    other.next.store(nullptr, std::memory_order_relaxed);
    next.store(n, std::memory_order_relaxed);
  }

  // Copy assignment: same issue — base class rebuilds flat_ from scratch.
  // `next` is NOT touched: copy-assignment leaves *this's existing chain
  // intact (a copy creates an independent chunk).
  LSetChunk& operator=(const LSetChunk& other) {
    if (this != &other) {
      writable_set<LObject, CompareLObject>::operator=(other);
      seq = other.seq;
      live_count_ = static_cast<int>(physical_size());
      deleted_count_ = 0;
      peak_deleted_count_ = 0;
      erase_call_count_ = 0;
      compact_call_count_ = 0;
      pop_skip_count_ = 0;
      rebuild_sev_flat();
      rebuild_sevSig_flat();
      rebuild_pair_index();
    }
    return *this;
  }

  // Move assignment: data members are moved; `next` is taken over from
  // the source (and the source's reset to nullptr) so move-assigned
  // objects preserve the chain.
  LSetChunk& operator=(LSetChunk&& other) noexcept {
    if (this != &other) {
      writable_set<LObject, CompareLObject>::operator=(std::move(other));
      seq = other.seq;
      live_count_ = other.live_count_;
      deleted_count_ = other.deleted_count_;
      peak_deleted_count_ = other.peak_deleted_count_;
      erase_call_count_ = other.erase_call_count_;
      compact_call_count_ = other.compact_call_count_;
      pop_skip_count_ = other.pop_skip_count_;
      other.live_count_ = 0;
      other.deleted_count_ = 0;
      rebuild_sev_flat();
      rebuild_sevSig_flat();
      rebuild_pair_index();
      LSetChunk* n = other.next.load(std::memory_order_relaxed);
      other.next.store(nullptr, std::memory_order_relaxed);
      next.store(n, std::memory_order_relaxed);
    }
    return *this;
  }

  // Override insert to maintain pair_index, sev_flat_, and sevSig_flat_
  iterator insert(const LObject& lobject) {
    iterator it = writable_set<LObject, CompareLObject>::insert(lobject);
    live_count_++;
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
    live_count_ = 0;
    deleted_count_ = 0;
  }

  // Override reorder to rebuild pair_index, sev_flat_, and sevSig_flat_
  // reorder() compacts first — tombstones are pruned then rebuilt.
  void reorder() {
    // Compact tombstones first so reorder sees a clean array.
    if (deleted_count_ > 0) compact();
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
  // physical_erase: the pre-task-360 erase contract.  Removes the
  // tree node and frees the LObject's polys.  Returns the iterator
  // following the erased entry (matches std::multiset::erase /
  // LSetChunk::erase convention).  Out-of-line for access to
  // file-static kLSet_free_polys.  Caller must hold wrlock or be
  // single-threaded — physical removal is unsafe under rdlock-shared
  // (UAF in concurrent reader scans).  The post-task-360 erase()
  // overloads above are tombstone-only because they're called from
  // worker rdlock paths; physical_erase is the named primitive for
  // serial code paths and bulk-cleanup that want immediate poly-free.
  iterator physical_erase(iterator it);

  // Free polys for every entry (live and tombstoned), then drop the
  // multiset, flat_, sev_flat_, sevSig_flat_, and pair_index.  O(N) in
  // chunk size — replaces the O(N^2) loop of repeated tombstone-only
  // pop_and_erase()s at end-of-bba cleanup paths.  Restores the
  // poly-free contract that the original (pre-chunked) LSet::erase
  // provided.  Caller must hold wrlock or be single-threaded.
  void clear_and_erase();

  // Read the chunk's `seq` counter (for LSet::would_be_top to stamp
  // a probe LObject without a full push).  Task 360, step 6c.
  unsigned peek_seq() const { return seq; }

  // Absorb all live entries from `src` into `*this` (task 360, step 6e).
  // Used by LSet::compact to coalesce chunks.  Walks src's live entries
  // (skipping tombstones) and pushes each into *this via the regular
  // insert path, which maintains pair_index, sev_flat_, sevSig_flat_.
  // After the call, `src` is emptied via clear() (its polys have been
  // transferred by-value into *this through the LObject struct copy
  // during insert).
  //
  // Note: this is NOT std::multiset::merge — that would be a node
  // splice but data_ is private to writable_set and the splice would
  // not update flat_/sev_flat_/sevSig_flat_/pair_index either way.
  // Iterate-and-insert at compact time costs O(n log n) which is
  // acceptable inside the brief writer-lock window.  See report.
  void absorb_live(LSetChunk& src);
};

/** @class LSet
 *
 * Wrapper around a linked list of LSetChunk instances (task 360,
 * step 6).  The head chunk is owned by value (`chunk_`); successor
 * chunks live on the heap and are linked via the atomic
 * `LSetChunk::next` pointer.
 *
 * Concurrency model (steady state under parallel-bba):
 *
 *   - Read-lock (`rdlock`): held by chainCritNormal scans, B-merge
 *     CAS-appends, pop, and the wrapper's own queries (size, empty,
 *     etc.).  Multiple readers run concurrently.  CAS on the tail
 *     chunk's `next` pointer admits new chunks while readers walk
 *     the chain.
 *   - Write-lock (`wrlock`): held only by `compact()`, which
 *     coalesces all chunks back into the head and trims tombstones.
 *     Compact runs ~once or twice per second on staging-9454; the
 *     window is brief (~5 ms estimated).
 *
 * The wrapper's `iterator` and `filtered_iterator` types span chunks
 * — `++` advances within the current chunk and crosses to
 * `chunk->next` when at chunk end.  `end()` is the sentinel past
 * the last chunk.  Iterator validity is bounded by the read-lock:
 * compact (writer) frees successor chunks, so callers must not hold
 * a wrapper iterator past the read-lock-release.  The read-lock
 * window is held by the caller side (worker phase 1: chainCritNormal
 * scan + erases + kMergeBintoL all under one rdlock window in
 * kthread.cc); single-call ops (pop, empty, ...) acquire internally.
 *
 * Memory ownership: head chunk is `chunk_` (by value).  Successor
 * chunks are raw `LSetChunk*` allocated by `new` and freed either
 * during `compact()` (after their contents are absorbed into head)
 * or in `~LSet()`.
 *
 * Public interface mirrors what the single-chunk LSet exposed;
 * callers in kthread.cc / kstdfac.cc / kutil.cc need no change.
 * The wrapper-level `iterator` type is distinct from
 * `LSetChunk::iterator` (it carries a chunk pointer); call sites
 * that took `LSetChunk::iterator` typed bvecs / lt locals had their
 * types updated to `LSet::iterator` in step 6a.
 */
class LSet {
public:
  // Compact-policy thresholds (task 360, step 6f).
  // - tomb-fraction trigger preserves the pre-step-6 policy:
  //   tombstones > live AND tombstones > 1024.
  // - chunk-count trigger bounds the O(k) pop scan and per_chunk
  //   pair_index walk.  Picked from the prompt.
  static constexpr int COMPACT_TOMB_MIN_ABS = 1024;
  static constexpr size_t COMPACT_CHUNK_COUNT_THRESHOLD = 256;

  // ----- Wrapper iterator type (task 360, step 6a) ---------------
  //
  // Spans chunks: holds a chunk pointer + an in-chunk iterator.
  // operator++ advances within the chunk; on inner==chunk->end()
  // jumps to the next chunk's begin (which is tombstone-aware) and
  // skips empty chunks.  end() is the (nullptr, default) sentinel.
  //
  // operator+ and operator+= are forwarded to repeated ++ — same
  // semantics as the underlying writable_set iterator, but
  // chunk-aware.
  //
  // Validity: a wrapper iterator is invalidated when the chunk it
  // points into is freed by compact().  The read-lock excludes
  // compact, so iterators are stable for the duration of any
  // caller's read-lock window.
  class iterator {
    friend class LSet;
    LSetChunk*           chunk_ = nullptr;
    LSetChunk::iterator  inner_;
    // Skip past empty/end chunks: if inner_ is at chunk_->end(),
    // walk forward through the chain until either inner_ points at
    // a live entry or chunk_ becomes nullptr (end of wrapper).
    void cross_chunk_advance() {
      while (chunk_ != nullptr && inner_ == chunk_->end()) {
        chunk_ = chunk_->next.load(std::memory_order_acquire);
        if (chunk_ != nullptr) inner_ = chunk_->begin();
      }
    }
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type        = LObject;
    using difference_type   = std::ptrdiff_t;
    using pointer           = LObject*;
    using reference         = LObject&;
    iterator() = default;
    iterator(LSetChunk* chunk, LSetChunk::iterator inner)
      : chunk_(chunk), inner_(inner) {
      cross_chunk_advance();
    }
    LSetChunk* chunk() const { return chunk_; }
    LSetChunk::iterator inner() const { return inner_; }
    reference operator*()  const { return *inner_; }
    pointer   operator->() const { return inner_.operator->(); }
    iterator& operator++() {
      ++inner_;
      cross_chunk_advance();
      return *this;
    }
    iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
    iterator operator+(difference_type n) const {
      iterator r = *this;
      for (difference_type i = 0; i < n; ++i) ++r;
      return r;
    }
    iterator& operator+=(difference_type n) {
      for (difference_type i = 0; i < n; ++i) ++(*this);
      return *this;
    }
    bool operator==(const iterator& o) const {
      // end() vs end() compare equal regardless of inner.  Two
      // non-end iterators must be in the same chunk and the same
      // inner position.
      if (chunk_ == nullptr && o.chunk_ == nullptr) return true;
      if (chunk_ != o.chunk_) return false;
      return inner_ == o.inner_;
    }
    bool operator!=(const iterator& o) const { return !(*this == o); }
  };

  // ----- Wrapper filtered_iterator (task 360, step 6a) -----------
  //
  // Crosses chunks like `iterator`.  Each chunk has its own
  // sev_flat_/sevSig_flat_ array; the wrapper iterator advances
  // within a chunk via the per-chunk filtered_iterator and crosses
  // at chunk boundary.  Erase dispatches into the underlying chunk.
  class filtered_iterator {
    friend class LSet;
    LSetChunk*                   chunk_     = nullptr;
    LSetChunk::filtered_iterator inner_;
    // Filter parameters (so we can build a fresh inner when
    // crossing into a new chunk): mode is encoded by sev2_==0 vs !=0.
    unsigned long                sev1_      = 0;
    unsigned long                sev2_      = 0;
    bool                         use_sig_   = false;  // false: sev_lcm; true: sevSig
    // Cross-chunk advance: inner_ has reached its chunk's end (per
    // pos_ == sev_array_->size()), advance to next chunk's begin
    // with the same filter, skip empty chunks.
    void cross_chunk_advance() {
      while (chunk_ != nullptr && inner_ == chunk_end_iter()) {
        chunk_ = chunk_->next.load(std::memory_order_acquire);
        if (chunk_ != nullptr) {
          inner_ = chunk_begin_iter();
        }
      }
    }
    LSetChunk::filtered_iterator chunk_begin_iter() const {
      if (use_sig_)         return chunk_->ufbegin_sig(sev1_);
      if (sev2_ == 0)       return chunk_->ufbegin_lcm(sev1_);
      return                       chunk_->ufbegin_lcm(sev1_, sev2_);
    }
    LSetChunk::filtered_iterator chunk_end_iter() const {
      return use_sig_ ? chunk_->ufend_sig() : chunk_->ufend_lcm();
    }
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type        = LObject;
    using difference_type   = std::ptrdiff_t;
    using pointer           = LObject*;
    using reference         = LObject&;
    filtered_iterator() = default;
    filtered_iterator(LSetChunk* chunk, LSetChunk::filtered_iterator inner,
                      unsigned long sev1, unsigned long sev2, bool use_sig)
      : chunk_(chunk), inner_(inner), sev1_(sev1), sev2_(sev2),
        use_sig_(use_sig) {
      cross_chunk_advance();
    }
    LSetChunk* chunk() const { return chunk_; }
    LSetChunk::filtered_iterator inner() const { return inner_; }
    reference operator*()  const { return *inner_; }
    pointer   operator->() const { return inner_.operator->(); }
    filtered_iterator& operator++() {
      ++inner_;
      cross_chunk_advance();
      return *this;
    }
    filtered_iterator operator++(int) {
      filtered_iterator tmp = *this; ++(*this); return tmp;
    }
    bool operator==(const filtered_iterator& o) const {
      if (chunk_ == nullptr && o.chunk_ == nullptr) return true;
      if (chunk_ != o.chunk_) return false;
      return inner_ == o.inner_;
    }
    bool operator!=(const filtered_iterator& o) const { return !(*this == o); }
  };

private:
  // Head chunk: owned by value.  Successor chunks (if any) live on
  // the heap and are reached via head_chunk_.next, ...->next, etc.
  // After every compact(), only chunk_ (head) holds entries; all
  // successor chunks have been freed.
  LSetChunk chunk_;

  // Tail-pointer hint for CAS-append.  Updated best-effort by
  // appenders; readers that need the precise tail walk the chain.
  // After compact(): tail_ == &chunk_.
  std::atomic<LSetChunk*> tail_{&chunk_};

  // Reader-writer lock guarding the chain.  Initialized
  // writer-preferring at construction so that compact (a writer)
  // doesn't get starved by a steady stream of readers.
  pthread_rwlock_t rwlock_;

  // Compact statistics (wrapper-level, task 360, step 6e).
  long compact_total_ns_ = 0;
  long compact_max_ns_   = 0;
  long compact_call_count_wrapper_ = 0;
  size_t last_compact_k_ = 1;

  // Helpers -----------------------------------------------------------

  // Walk the chunk chain starting at &chunk_, count the chunks.
  size_t chunk_count_unlocked() const {
    size_t k = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      ++k;
    }
    return k;
  }

  // Free all successor chunks (everything after head).  After this
  // call, chunk_.next == nullptr and tail_ == &chunk_.  Caller must
  // ensure chunks have been emptied (their polys absorbed into
  // head) before calling — otherwise polys leak.  Used only by
  // compact() under wrlock.
  void free_successors_unlocked() {
    LSetChunk* c = chunk_.next.load(std::memory_order_acquire);
    chunk_.next.store(nullptr, std::memory_order_release);
    while (c != nullptr) {
      LSetChunk* nxt = c->next.load(std::memory_order_acquire);
      delete c;
      c = nxt;
    }
    tail_.store(&chunk_, std::memory_order_release);
  }

  // O(k) global-min scan over chunk tops (task 360, step 6c).
  // Returns wrapper iterator to the global minimum, or end() if
  // empty.  Caller must hold rdlock or wrlock.
  //
  // Bench toggle (task 363; run 587): when g_bench_elide_lset_wrapper
  // is set AND we're in serial mode (kt_current_ctx == NULL), skip
  // the chunk-walk and go straight to the head chunk.  At T=1 there
  // is only ever one chunk (compact runs in refill_and_publish which
  // is parallel-only; serial bba's compact at exitBuchMora has
  // nothing to do because we already physical-erase), so this is
  // semantically identical and bypasses the iterator + acquire-load
  // indirection.
  iterator find_global_min_unlocked() {
    if (g_bench_elide_lset_wrapper && kt_current_ctx == NULL) {
      LSetChunk::iterator it = chunk_.begin();
      if (it == chunk_.end()) return iterator();
      return iterator(&chunk_, it);
    }
    LSetChunk*           winner = nullptr;
    LSetChunk::iterator  winner_it;
    for (LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      auto it = c->begin();           // tombstone-skipping
      if (it == c->end()) continue;   // empty chunk
      if (winner == nullptr ||
          chunk_.key_comp()(*it, *winner_it)) {
        winner = c; winner_it = it;
      }
    }
    return iterator(winner, winner ? winner_it : LSetChunk::iterator());
  }

public:
  using iterator_t         = iterator;       // disambiguate inner alias
  using const_iterator     = iterator;       // const reads use the same wrapper
  using reverse_iterator   = LSetChunk::reverse_iterator;
  using unordered_iterator = LSetChunk::unordered_iterator;
  using size_type          = LSetChunk::size_type;
  using value_type         = LSetChunk::value_type;
  using difference_type    = LSetChunk::difference_type;

  // Pair-index map type forwarded for cleanTSbaRing's direct
  // wrapper-level access.  After step 6d cleanTSbaRing uses
  // pair_index_find() / pair_index_walk_size() instead of indexing
  // a single chunk.
  using pair_index_t = decltype(LSetChunk::pair_index);

  // ----- Construction / destruction -----
  LSet() {
    pthread_rwlockattr_t attr;
    pthread_rwlockattr_init(&attr);
#if defined(__GLIBC__) && defined(PTHREAD_RWLOCK_PREFER_WRITER_NONRECURSIVE_NP)
    pthread_rwlockattr_setkind_np(&attr,
        PTHREAD_RWLOCK_PREFER_WRITER_NONRECURSIVE_NP);
#endif
    pthread_rwlock_init(&rwlock_, &attr);
    pthread_rwlockattr_destroy(&attr);
  }
  ~LSet() {
    // Free any successor chunks that survived past last compact (if
    // we shut down before compact ran).  Each LSetChunk destructor
    // frees its multiset's polys via the writable_set destructor's
    // clear() call; we don't re-call free here.
    LSetChunk* c = chunk_.next.load(std::memory_order_relaxed);
    while (c != nullptr) {
      LSetChunk* nxt = c->next.load(std::memory_order_relaxed);
      delete c;
      c = nxt;
    }
    pthread_rwlock_destroy(&rwlock_);
  }
  // Non-copyable / non-movable: pthread_rwlock_t is not movable, and
  // strat->L is logically a fixed-position container anyway.  Use
  // copy_contents_from for the (rare) case that callers want to
  // duplicate L's contents.
  LSet(const LSet&)            = delete;
  LSet& operator=(const LSet&) = delete;
  LSet(LSet&&)                 = delete;
  LSet& operator=(LSet&&)      = delete;

  // Copy chunk contents from another LSet, preserving rwlock_
  // identity.  Used by kstdfac.cc:copyL to clone a strategy's L.
  // Only the head chunk is copied — the source must have been
  // compacted (single chunk) for this to be a complete copy.
  // Serial code paths only.
  void copy_contents_from(const LSet& other) {
    chunk_ = other.chunk_;
    chunk_.next.store(nullptr, std::memory_order_relaxed);
    tail_.store(&chunk_, std::memory_order_relaxed);
  }

  // ----- Pair-index lookup (task 360, step 6d) ---------------------
  //
  // Walk the chunk chain, returning a per-chunk iterator + chunk
  // pointer for the first hit.  Returns (nullptr, default) if no
  // chunk has the key.  Caller must hold rdlock or wrlock.
  // Wrapped to bump pair_index_lookup_ns under KTHREAD_INSTRUMENT
  // (see kutil.cc:isInPairsetL / cleanTSbaRing).
  iterator pair_index_find(const std::pair<poly, poly>& key) {
    for (LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      auto found = c->pair_index.find(key);
      if (found != c->pair_index.end())
        return iterator(c, found->second);
    }
    return iterator();  // end()
  }

  // Erase a pair_index entry under wrapper, dispatching to the
  // owning chunk.  Used by cleanTSbaRing after it has identified
  // a key to forget.
  void pair_index_erase(const std::pair<poly, poly>& key) {
    for (LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      auto found = c->pair_index.find(key);
      if (found != c->pair_index.end()) {
        c->pair_index.erase(found);
        return;
      }
    }
  }

  // Aggregate pair_index size (sum across chunks).  Diagnostic.
  size_t pair_index_walk_size() const {
    size_t n = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      n += c->pair_index.size();
    }
    return n;
  }

  // ----- Iterator factories ----------------------------------------
  iterator begin() { return iterator(&chunk_, chunk_.begin()); }
  iterator end()   { return iterator(); }   // (nullptr, default)
  // Const begin/end use the same wrapper (read-only deref).
  const_iterator cbegin() const {
    // const-correctness: chunk_'s next.load is allowed; we cast away
    // const for the wrapper iterator since LObject reads still go
    // through a non-const path internally (legacy behaviour).
    return const_cast<LSet*>(this)->begin();
  }
  const_iterator cend() const { return const_iterator(); }

  // Filtered iteration (lcm divisibility / incomparability and
  // signature divisibility) — wrapper-level, crosses chunks.
  filtered_iterator ufbegin_lcm(unsigned long sev) {
    return filtered_iterator(&chunk_, chunk_.ufbegin_lcm(sev),
                             sev, 0, /*use_sig=*/false);
  }
  filtered_iterator ufbegin_lcm(unsigned long sev1, unsigned long sev2) {
    return filtered_iterator(&chunk_, chunk_.ufbegin_lcm(sev1, sev2),
                             sev1, sev2, /*use_sig=*/false);
  }
  filtered_iterator ufend_lcm() { return filtered_iterator(); }
  filtered_iterator ufbegin_sig(unsigned long sev) {
    return filtered_iterator(&chunk_, chunk_.ufbegin_sig(sev),
                             sev, 0, /*use_sig=*/true);
  }
  filtered_iterator ufend_sig() { return filtered_iterator(); }

  // ----- Push (single-element insert into head) -------------------
  //
  // Used by serial code paths and by ring strategies (enterOneStrongPoly
  // etc.) that push one entry at a time into L.  Always goes into
  // the head chunk to keep the wrapper interface intact.  In
  // parallel mode the worker phase-1 path uses the move-append
  // (`append_chunk`) instead.
  iterator push(LObject& lobject) {
    auto inner = chunk_.push(lobject);
    return iterator(&chunk_, inner);
  }

  // would_be_top: would `lobject` be at the global minimum if
  // pushed?  Compute by: if any existing entry is smaller than
  // `lobject` under the comparator, the answer is no.  The check
  // mutates lobject.seq (matches LSetChunk::would_be_top semantics).
  bool would_be_top(LObject& lobject) {
    if (empty()) {
      // Defer to chunk's would_be_top so seq is bumped consistently.
      return chunk_.would_be_top(lobject);
    }
    // Find the global min top, compare against lobject.
    iterator min_it = find_global_min_unlocked();
    // mimic chunk's seq stamping (LSetChunk::would_be_top sets seq).
    lobject.seq = chunk_.peek_seq();
    return chunk_.key_comp()(lobject, *min_it);
  }

  // Serial compact-on-pop hook (task 368).  Called at the END of a
  // serial-mode pop()/pop_and_erase(), AFTER the physical erase has
  // completed and every local iterator (winner_it, etc.) has gone out
  // of scope, before returning.  compact() invalidates ALL outstanding
  // strat->L iterators — every serial pop()/pop_and_erase() caller
  // (bba/sba serial driver in kstd2.cc, pop_and_prepare in kthread.cc)
  // does `P = top(); pop();` with NO strat->L iterator held across the
  // pop, so this is safe.  Serial mode is single-threaded (no
  // concurrent rdlock/wrlock holders), so calling LSet::compact()
  // directly without the wrlock is safe (the wrlock contract exists
  // only to exclude concurrent readers, of which there are none at
  // T=1; find_global_min_unlocked's comment already treats serial
  // compact as safe, just never-triggered).  No-op unless the knob is
  // enabled AND we are in serial mode.
  inline void maybe_serial_compact_on_pop(void) {
    if (g_bench_serial_compact == 0) return;          // off: true no-op
    if (kt_current_ctx != NULL) return;               // serial only
    if (g_bench_serial_compact == 1) {                // threshold
      if (!needs_compact()) return;
    }
    // mode 2 (everypop): always — compact() early-returns when
    // deleted_count_ == 0, so this is "compact if anything to compact".
    compact();
  }

  // ----- Top / Pop -------------------------------------------------
  const LObject& top(void) {
    iterator it = find_global_min_unlocked();
    return *it;
  }

  void pop(void) {
    // Pop the global minimum.  Walk to find which chunk hosts it,
    // then call that chunk's pop() to physically remove.
    //
    // Bench toggle (task 363; run 587): serial-mode short-circuit
    // to head chunk — see find_global_min_unlocked for rationale.
    if (g_bench_elide_lset_wrapper && kt_current_ctx == NULL) {
      if (!chunk_.empty()) chunk_.pop();
      maybe_serial_compact_on_pop();   // task 368 (serial only)
      return;
    }
    {
      LSetChunk*           winner = nullptr;
      LSetChunk::iterator  winner_it;
      for (LSetChunk* c = &chunk_; c != nullptr;
           c = c->next.load(std::memory_order_acquire)) {
        auto it = c->begin();
        if (it == c->end()) continue;
        if (winner == nullptr || chunk_.key_comp()(*it, *winner_it)) {
          winner = c; winner_it = it;
        }
      }
      if (winner != nullptr) winner->pop();
    }   // winner_it out of scope before compact()
    maybe_serial_compact_on_pop();     // task 368 (serial only)
  }

  void pop_and_erase(void) {
    // Tombstone-erase the global minimum: same as pop but uses
    // erase() so polys are freed at compact time (vs immediately
    // by pop, which transfers them to the caller).  Used at error-
    // path L-clear (kthread.cc parallel_shutdown).
    //
    // Bench toggle (task 363; run 587): serial-mode short-circuit.
    if (g_bench_elide_lset_wrapper && kt_current_ctx == NULL) {
      if (!chunk_.empty()) chunk_.pop_and_erase();
      maybe_serial_compact_on_pop();   // task 368 (serial only)
      return;
    }
    {
      LSetChunk*           winner = nullptr;
      LSetChunk::iterator  winner_it;
      for (LSetChunk* c = &chunk_; c != nullptr;
           c = c->next.load(std::memory_order_acquire)) {
        auto it = c->begin();
        if (it == c->end()) continue;
        if (winner == nullptr || chunk_.key_comp()(*it, *winner_it)) {
          winner = c; winner_it = it;
        }
      }
      if (winner != nullptr) winner->pop_and_erase();
    }   // winner_it out of scope before compact()
    maybe_serial_compact_on_pop();     // task 368 (serial only)
  }

  // ----- Erase via wrapper iterator -------------------------------
  //
  // Dispatches between physical erase (serial mode) and tombstone-only
  // erase (parallel mode, where worker rdlock holders would UAF on
  // physical removal).  kt_current_ctx is non-NULL only inside
  // bba_parallel_loop (kthread.cc:2070 / :2373).  Single load + branch
  // per call.
  iterator erase(iterator it) {
    if (kt_current_ctx == NULL) return physical_erase(it);
    if (it.chunk_ == nullptr) return it;
    LSetChunk* c = it.chunk_;
    auto next_inner = c->erase(it.inner_);
    // The new wrapper iterator points at the next entry in the
    // SAME chunk; if that's chunk-end, cross_chunk_advance moves
    // forward through the chain.
    return iterator(c, next_inner);
  }

  // physical_erase: pre-task-360 erase semantics (physical tree
  // removal + immediate poly-free).  Used directly by serial-only
  // callers that want to bypass the dispatch in erase().  The
  // dispatch in erase() routes here automatically when
  // kt_current_ctx == NULL.
  iterator physical_erase(iterator it) {
    if (it.chunk_ == nullptr) return it;
    LSetChunk* c = it.chunk_;
    auto next_inner = c->physical_erase(it.inner_);
    return iterator(c, next_inner);
  }

  // erase via wrapper filtered_iterator.  Default (tombstone) path
  // stays tombstone-only — filtered_iterator is the worker-rdlock
  // chunked-scan pattern, so by convention even when kt_current_ctx
  // == NULL a filtered_iterator caller is inside a chunked-LSet scan
  // that doesn't want physical removal.
  //
  // Task 370 knob: when SINGULAR_BENCH_SERIAL_LSET_ERASE=physical
  // (g_bench_serial_lset_erase == PHYSICAL) AND we are in serial mode
  // (kt_current_ctx == NULL), route through physical_erase instead —
  // the A/B comparator for per-delete physical erase (next-opt
  // contract) vs deferred tombstone + everypop compact (task 368).
  // The gate is one global load + one branch on the hot path; the
  // default (g_bench_serial_lset_erase == 0) compiles to exactly the
  // pre-task-370 tombstone code (true no-op for non-bench / default
  // builds).  Parallel mode (kt_current_ctx != NULL) is never routed
  // to physical_erase, so T>1 is byte- and TSan-identical.
  filtered_iterator erase(filtered_iterator fit) {
    if (fit.chunk_ == nullptr) return fit;
    if (g_bench_serial_lset_erase == KT_SERIAL_LSET_ERASE_PHYSICAL
        && kt_current_ctx == NULL)
      return physical_erase(fit);
    LSetChunk* c = fit.chunk_;
    auto next_inner = c->erase(fit.inner_);
    fit.inner_ = next_inner;
    fit.cross_chunk_advance();
    return fit;
  }

  // physical_erase via wrapper filtered_iterator (task 370): the
  // pre-task-360 physical-erase contract reached through the
  // chunk-aware sev scan.  Forwards to the chunk's
  // physical_erase(filtered_iterator) (immediate tree-node removal +
  // poly-free), then re-establishes the wrapper filtered_iterator
  // (cross-chunk-aware) at the next position — mirroring exactly how
  // the tombstone erase(filtered_iterator) forwards to the chunk.
  // Serial-only / wrlock-only (the chunk method asserts
  // kt_current_ctx == NULL); reached only via the gated erase() above
  // or by an explicit serial caller.
  filtered_iterator physical_erase(filtered_iterator fit) {
    if (fit.chunk_ == nullptr) return fit;
    LSetChunk* c = fit.chunk_;
    auto next_inner = c->physical_erase(fit.inner_);
    fit.inner_ = next_inner;
    fit.cross_chunk_advance();
    return fit;
  }

  // ----- Clear / compact / reorder --------------------------------
  void clear() {
    // Free successor chunks then clear the head.
    free_successors_unlocked();
    chunk_.clear();
  }

  // O(N) bulk cleanup: walk every chunk, free each entry's polys via
  // kLSet_free_polys, drop all multiset / flat_ / sev_flat_ / pair_index
  // state, then free heap-allocated successor chunks.  Replaces the
  // O(N^2) `while (!L.empty()) L.pop_and_erase()` end-of-bba loops.
  // Caller must hold wrlock or be single-threaded (the parallel-bba
  // post-join cleanup path is single-threaded by construction).
  void clear_and_erase();

  // compact(): coalesce all chunks into the head chunk and drop
  // tombstones.  Caller must hold the writer lock.
  // Implementation in kutil.cc.
  void compact();

  void reorder() {
    // reorder() can only sensibly run when there's a single chunk
    // (its semantics rebuild flat_index for one set).  If we have
    // multiple chunks, compact first.
    if (chunk_.next.load(std::memory_order_acquire) != nullptr) compact();
    chunk_.reorder();
  }

  // ----- Append a finalized chunk (task 360, step 6b) -------------
  //
  // Take ownership of `src` (move-from) and link it at the tail
  // via CAS.  Returns the heap pointer of the newly-published
  // chunk (or nullptr if src was empty so nothing was published).
  // Caller must hold rdlock; the rdlock excludes only compact,
  // not other appenders, so the CAS handles appender-vs-appender
  // races.
  LSetChunk* append_chunk(LSetChunk&& src) {
    if (src.size() == 0) {
      // Even if empty, src may carry tombstones; reset for caller's
      // convenience.
      src.clear();
      return nullptr;
    }
    // Heap-allocate via move.  `src` is left empty after the move.
    LSetChunk* fresh = new LSetChunk(std::move(src));
    fresh->next.store(nullptr, std::memory_order_relaxed);
    // CAS the tail-chunk's next.  If a peer appended ahead of us,
    // walk forward to a chunk whose `next` is null and CAS there.
    LSetChunk* cur = tail_.load(std::memory_order_acquire);
    while (true) {
      // Walk to a chunk with next==nullptr (the actual tail; the
      // tail_ hint may be stale).
      while (true) {
        LSetChunk* nx = cur->next.load(std::memory_order_acquire);
        if (nx == nullptr) break;
        cur = nx;
      }
      LSetChunk* expected = nullptr;
      if (cur->next.compare_exchange_weak(
            expected, fresh,
            std::memory_order_release,
            std::memory_order_acquire)) {
        // Successfully linked.  Update the hint best-effort —
        // readers tolerate a stale hint (they walk the chain).
        tail_.store(fresh, std::memory_order_release);
        return fresh;
      }
      // CAS failed: another appender beat us.  cur stays the same
      // — the next iteration's inner loop walks past the new
      // successor before retrying.
    }
  }

  // ----- Aggregate counters (walk-on-demand) ----------------------
  //
  // For all the "size-like" queries we walk the chunk chain and sum
  // per-chunk values.  k bounded to 256 by compact policy, so each
  // query is O(k) atomic loads — under 1 µs on modern CPU.
  size_type size() const {
    size_type s = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      s += c->size();
    }
    return s;
  }
  bool empty() const {
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      if (!c->empty()) return false;
    }
    return true;
  }
  size_type physical_size() const {
    size_type s = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      s += c->physical_size();
    }
    return s;
  }
  size_t sev_flat_size() const {
    size_t s = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      s += c->sev_flat_size();
    }
    return s;
  }

  // Comparator access — used by callers in chainCritNormal /
  // kMergeBintoL_and_return_iterators.  All chunks share the same
  // comparator (CompareLObject::strat == strat), so head's
  // key_comp is representative.  Note: const overload of
  // writable_set::key_comp returns by value (matches std::set), so
  // we mirror that here.
  CompareLObject&  key_comp()        { return chunk_.key_comp(); }
  CompareLObject   key_comp() const  { return chunk_.key_comp(); }

  // Tombstone diagnostics (aggregates across chunks).
  int  deleted_count() const {
    int d = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      d += c->deleted_count();
    }
    return d;
  }
  int  peak_deleted_count() const  { return chunk_.peak_deleted_count(); }
  long erase_call_count() const {
    long e = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      e += c->erase_call_count();
    }
    return e;
  }
  long compact_call_count() const  { return compact_call_count_wrapper_; }
  long compact_total_ns() const    { return compact_total_ns_; }
  long compact_max_ns() const      { return compact_max_ns_; }
  size_t last_compact_k() const    { return last_compact_k_; }
  long pop_skip_count() const {
    long p = 0;
    for (const LSetChunk* c = &chunk_; c != nullptr;
         c = c->next.load(std::memory_order_acquire)) {
      p += c->pop_skip_count();
    }
    return p;
  }
  void debug_print_stats(const char *tag = NULL) const { chunk_.debug_print_stats(tag); }

  // Chunk count (for the diagnostic dump).  O(k) walk.
  size_t chunk_count() const { return chunk_count_unlocked(); }

  // ----- Compact threshold ----------------------------------------
  // (task 360, step 6f).  Two trigger conditions, OR'd:
  //   - tomb-fraction:  tombstones > live  AND tombstones > 1024
  //   - chunk-count:    k > COMPACT_CHUNK_COUNT_THRESHOLD
  // Caller checks this AFTER releasing the rdlock; on true, takes
  // wrlock and calls compact() (which re-checks under wrlock).
  bool needs_compact() const {
    int d = deleted_count();
    int l = (int)size();
    if (d > l && d > COMPACT_TOMB_MIN_ABS) return true;
    if (chunk_count_unlocked() > COMPACT_CHUNK_COUNT_THRESHOLD) return true;
    return false;
  }

  // ----- Reader-writer lock interface -----
  //
  // Concurrency discipline (steady state under parallel-bba):
  //   - Multi-call L-touching blocks (worker phase 1: chainCritNormal
  //     scan + erases + kMergeBintoL): the **caller** wraps the block
  //     in rdlock() / unlock().  Wrapper methods invoked inside the
  //     block do NOT acquire — iterator validity is preserved across
  //     the block.
  //   - Single-call ops (empty, size, top, pop, ...): caller may invoke
  //     directly, but for safety against a concurrent compact under
  //     write-lock, the caller-side wrapper IS expected to hold the
  //     rdlock.  The current call-site discipline (in kthread.cc and
  //     kutil.cc parallel paths) is "every site that touched L_lock
  //     pre-step-5 now takes rdlock"; non-parallel paths (exitBuchMora,
  //     bba serial driver) need no locking because no concurrency.
  //   - compact(): caller must hold the wrlock.  The natural pattern
  //     is rdlock-pop-unlock, then check needs_compact(), then
  //     wrlock-compact-unlock (drop-and-reacquire).  See
  //     fill_active_slots in kthread.cc.
  void rdlock()   { pthread_rwlock_rdlock(&rwlock_); }
  void wrlock()   { pthread_rwlock_wrlock(&rwlock_); }
  void unlock()   { pthread_rwlock_unlock(&rwlock_); }
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
  // atS: iterator to the position in S where h should be inserted.
  //   Pass strat->S.end() to mean "compute position via find_pos internally"
  //   (previously encoded as int atS = -1).
  //   Other iterator values pin the insertion position explicitly
  //   (previously encoded as int atS >= 0).
  sBasisSet::iterator (*enterS)(LObject &h, kStrategy strat, int atR /*= -1*/, sBasisSet::iterator atS) = NULL;
  void (*initEcartPair)(LObject * h, poly f, poly g, int ecartF, int ecartG) = NULL;
  void (*enterOnePair) (const SElement &si,poly p,int ecart, int isFromQ,kStrategy strat, int atR /*= -1*/) = NULL;
  void (*chainCrit) (poly p,int ecart,kStrategy strat) = NULL;
  BOOLEAN (*syzCrit) (poly sig, unsigned long not_sevSig, kStrategy strat) = NULL;
  // rewCrit{1,2,3}: rewritten-criterion scans of strat->S for SBA/F5.
  // start is an iterator pointing at the first element to check
  // (inclusive).  Callers pass strat->S.cbegin() for a full scan, or a
  // captured iterator (e.g. sLObject::checked) for a resume-from-snapshot
  // scan.  Converted from `int start` in Stage B of the sBasisSet
  // iterator migration so that stored snapshots (Lp.checked) survive
  // S appends without stale-index bugs under lazy-mode deletion.
  BOOLEAN (*rewCrit1) (poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start) = NULL;
  BOOLEAN (*rewCrit2) (poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start) = NULL;
  BOOLEAN (*rewCrit3) (poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start) = NULL;
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
  LSet      L;     // chunked-LSet wrapper (step 4: single-chunk no-op)
  LSetChunk B;     // serial-only B (kept as bare chunk; phase-1 drainers
                   // use thread-local local_B via t_local_B_override).
  // arrival_counter: monotonic counter incremented on every successful
  // enterS.  Used by the parallel phase-1 drainer to filter S iteration
  // to entries that arrived before the current survivor h — because
  // enterS places h at a sorted position, the insertion index does not
  // reflect arrival order.  Serial code also bumps the counter (keeps
  // the invariant that every S element has a unique arrival_id) but
  // never reads it.  Declared std::atomic<uint64_t> for parallel
  // fetch_add semantics; zero-cost in serial mode.
  std::atomic<uint64_t> arrival_counter{0};
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

  // S-to-T lookup is now sBasisSet::S_2_T / s_2_t — call via strat->S.
};

// Inline definition of strat_B.  Returns the thread-local override B if
// set (parallel phase-1 drain), else the shared strat->B (serial path).
static inline LSetChunk& strat_B(kStrategy strat) {
  return t_local_B_override ? *t_local_B_override : strat->B;
}

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
void cleanT (kStrategy strat);
sBasisSet::iterator enterSBba (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS);
sBasisSet::iterator enterSBbaShift (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS);
sBasisSet::iterator enterSSba (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS);
void initEcartPairBba (LObject* Lp,poly f,poly g,int ecartF,int ecartG);
void initEcartPairMora (LObject* Lp,poly f,poly g,int ecartF,int ecartG);
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


int posInSyz (const kStrategy strat, const poly sig);
// redtailBba overloads: the `end` iterator is an *exclusive* upper bound
// on the range of strat->S considered for reduction (STL-style:
// reduce against [strat->S.begin(), end)). Previous int end_pos was
// *inclusive* — translation for callers: old `pos-1` (inclusive) becomes
// new `pos` (exclusive iterator); old `strat->S.size()-1` becomes
// `strat->S.end()`.
KINLINE poly redtailBba (poly p,sBasisSet::const_iterator end,kStrategy strat,BOOLEAN normalize=FALSE);
KINLINE poly redtailBbaBound (poly p,sBasisSet::const_iterator end,kStrategy strat,int bound,BOOLEAN normalize=FALSE);
KINLINE poly redtailBba_Ring (poly p,sBasisSet::const_iterator end,kStrategy strat);
KINLINE poly redtailBba_Z (poly p,sBasisSet::const_iterator end,kStrategy strat);
poly redtailBba_NF (poly p, kStrategy strat );
poly redtailBba_Ring (LObject* L, sBasisSet::const_iterator end, kStrategy strat );
poly redtailBba_Z (LObject* L, sBasisSet::const_iterator end, kStrategy strat );
// Task 361 (run 581): the parallel-bba phase-0(b) drain calls these
// concurrently across worker threads.  The two strat fields below are
// pure outputs of these calls (not inputs), so we expose them as
// optional per-call out parameters; when non-null the caller gets a
// thread-local copy and can avoid the redtail_lock.  When null, the
// existing serial-caller behaviour is preserved (writes through to
// strat).  See ~/project/docs/parallel-bba-redtail-decouple.md.
void redtailBbaAlsoLC_Z (LObject* L, kStrategy strat,
                         bool *out_redTailChange = nullptr,
                         bool *out_completeReduce_retry = nullptr);
poly redtailBba (LObject *L, sBasisSet::const_iterator end,kStrategy strat,
                 BOOLEAN withT = FALSE,BOOLEAN normalize=FALSE,
                 bool *out_redTailChange = nullptr,
                 bool *out_completeReduce_retry = nullptr);
poly redtailBbaBound (LObject *L, sBasisSet::const_iterator end,kStrategy strat,int bound,
                 BOOLEAN withT = FALSE,BOOLEAN normalize=FALSE);
poly redtailSba (LObject *L, sBasisSet::const_iterator end,kStrategy strat,
                 BOOLEAN withT = FALSE,BOOLEAN normalize=FALSE);
poly redtailBba (TObject *T, sBasisSet::const_iterator end,kStrategy strat);
poly redtail (poly p,sBasisSet::const_iterator end,kStrategy strat);
poly redtail (LObject *L,sBasisSet::const_iterator end,kStrategy strat);
poly redNF (poly h,int nonorm,kStrategy strat);
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
void superenterpairs (poly h,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR = -1);
void superenterpairsSig (poly h,poly hSig,int hFrom,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR = -1);
int redLazy (LObject* h,kStrategy strat);
int redHomog (LObject* h,kStrategy strat);
int redSig (LObject* h,kStrategy strat);
int redSigRing (LObject* h,kStrategy strat);
//adds hSig to be able to check with F5's criteria when entering pairs!
void enterpairsSig (poly h, poly hSig, int from, int k, int ec, sBasisSet::iterator pos,kStrategy strat, int atR = -1);
void enterpairs (poly h, int k, int ec, sBasisSet::iterator pos,kStrategy strat, int atR = -1);
void entersets (LObject h);
void pairs ();
BOOLEAN sbaCheckGcdPair (LObject* h,kStrategy strat);
void message (int i,int* olddeg,LSetChunk::size_type* reduc,kStrategy strat,int red_result);
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
KINLINE BOOLEAN arriRewDummy(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start);
BOOLEAN arriRewCriterion(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start);
BOOLEAN arriRewCriterionPre(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start);
BOOLEAN faugereRewCriterion(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator start);
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
/// return strat->S.end() if no divisor is found;
/// iterator to first divisor in S, otherwise
sBasisSet::iterator kFindDivisibleByInS(const kStrategy strat, LObject* L);

sBasisSet::iterator kFindNextDivisibleByInS(const kStrategy strat, sBasisSet::iterator start, LObject* L);
TObject* kFindDivisibleByInS_T(kStrategy strat, sBasisSet::const_iterator end, LObject* L, TObject *T, long ecart = LONG_MAX);

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
sBasisSet::iterator enterSMora (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS);
sBasisSet::iterator enterSMoraNF (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS);


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
          int& hilbcount, int& srmax, LSetChunk::size_type& reduc, ideal Q,
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
BOOLEAN kCheckStrongCreation(int atR, poly m1, sBasisSet::const_iterator atS, poly m2, kStrategy strat);
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

BOOLEAN enterOnePairShift (poly q, poly p, int ecart, int isFromQ, kStrategy strat, int atR, int ecartq, int qisFromQ, int shiftcount, sBasisSet::const_iterator ifromS);

void enterpairsShift (poly h,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR);

void superenterpairsShift (poly h,int k,int ecart,int pos,kStrategy strat, int atR);

poly redtailBbaShift (LObject* L, sBasisSet::const_iterator end, kStrategy strat, BOOLEAN withT, BOOLEAN normalize);

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
