/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/*
* ABSTRACT: kernel: utils for kStd
*/

// #define PDEBUG 2
// #define PDIV_DEBUG
#define KUTIL_CC

#define MYTEST 0

//All vs Just strategy over rings:
// 1 - Just
// 0 - All
#define ALL_VS_JUST 0
//Extended Spoly Strategy:
// 0 - new gen sig
// 1 - ann*old sig
#define EXT_POLY_NEW 0

#include "kernel/mod2.h"

#include "misc/mylimits.h"
#include "misc/options.h"
#include "polys/nc/nc.h"
#include "polys/nc/sca.h"
#include "polys/weight.h" /* for kDebugPrint: maxdegreeWecart*/

#include <stdlib.h>
#include <string.h>
#include <algorithm>

#ifdef KDEBUG
#undef KDEBUG
#define KDEBUG 2
#endif

#ifdef DEBUGF5
#undef DEBUGF5
//#define DEBUGF5 1
#endif

// define if enterT should use memmove instead of doing it manually
// on topgun, this is slightly faster (see monodromy_l.tst, homog_gonnet.sing)
#ifndef SunOS_4
#define ENTER_USE_MEMMOVE
#endif

// define, if the my_memmove inlines should be used instead of
// system memmove -- it does not seem to pay off, though
// #define ENTER_USE_MYMEMMOVE

#include "kernel/GBEngine/kutil.h"
#include "kernel/GBEngine/kthread.h"
#include "polys/kbuckets.h"
#include "coeffs/numbers.h"
#include "kernel/polys.h"
#include "polys/monomials/ring.h"
#include "kernel/ideals.h"
#include "kernel/combinatorics/stairc.h"
#include "kernel/GBEngine/kstd1.h"
#include "polys/operations/pShallowCopyDelete.h"

#ifdef HAVE_SHIFTBBA
#include "polys/shiftop.h"
#endif

#include "polys/prCopy.h"

#ifdef HAVE_RATGRING
#include "kernel/GBEngine/ratgring.h"
#endif

/* Task 325 parallel-bba-event-log: defer-frees wrapper for pLmFree /
 * pDelete / p_LmFree / p_Delete.  MUST be included AFTER all poly
 * headers so the real inline functions / macros are declared first;
 * this header then #define's them to the kt_* wrappers.  Untraced
 * runs: wrappers are zero-cost (a single load of g_defer_frees). */
#include "kernel/GBEngine/kevlog.h"
#include "kernel/GBEngine/kevlog_wrap.h"

#ifdef DEBUGF5
#undef DEBUGF5
#define DEBUGF5 2
#endif

VAR denominator_list DENOMINATOR_LIST=NULL;

// Thread-local override for strat->B.  NULL in serial mode; set by the
// parallel phase-1 drain to a stack-allocated LSet so concurrent drainers
// each build their own B.  See kutil.h for full rationale.
__thread LSet* t_local_B_override = NULL;

// Thread-local "my arrival" filter for the parallel phase-1 drain.
// UINT64_MAX in serial mode (no filter); set to the drainer's
// my_arrival before phase 1 so S-iteration skips entries with
// arrival_id >= my_arrival.  See kutil.h for rationale.
__thread uint64_t t_local_my_arrival = UINT64_MAX;

// Thread-local pairtest-hit list.  NULL in serial mode (falls back to
// SElement.pairtest); set by the parallel drain to a stack-allocated
// vector before phase 1.  See kutil.h for rationale.
__thread std::vector<SElement*>* t_local_pairtest_hits = NULL;

// Arrival-id filter for S iteration during the parallel phase-1 drain.
// Returns true if the entry should be processed (its arrival_id is
// strictly less than the current drainer's my_arrival).  In serial mode
// t_local_my_arrival == UINT64_MAX so this is always true.
static inline bool arrival_id_ok(const SElement &s) {
  return s.arrival_id < t_local_my_arrival;
}

// Record a pairtest hit.  In serial mode, sets SElement.pairtest
// (existing behaviour).  In parallel phase-1 mode, pushes the
// SElement* onto the thread-local hit list so that chainCritNormal
// can iterate the current drainer's hits only, without racing on
// SElement.pairtest across drainers.
static inline void record_pairtest_hit(const SElement &si, kStrategy strat) {
  if (t_local_pairtest_hits != NULL) {
    t_local_pairtest_hits->push_back(const_cast<SElement*>(&si));
  } else {
    selement_pairtest_set(si);
    strat->S.set_pairtest_any();
  }
}


#ifdef ENTER_USE_MYMEMMOVE
inline void _my_memmove_d_gt_s(unsigned long* d, unsigned long* s, long l)
{
  REGISTER unsigned long* _dl = (unsigned long*) d;
  REGISTER unsigned long* _sl = (unsigned long*) s;
  REGISTER long _i = l - 1;

  do
  {
    _dl[_i] = _sl[_i];
    _i--;
  }
  while (_i >= 0);
}

inline void _my_memmove_d_lt_s(unsigned long* d, unsigned long* s, long l)
{
  REGISTER long _ll = l;
  REGISTER unsigned long* _dl = (unsigned long*) d;
  REGISTER unsigned long* _sl = (unsigned long*) s;
  REGISTER long _i = 0;

  do
  {
    _dl[_i] = _sl[_i];
    _i++;
  }
  while (_i < _ll);
}

inline void _my_memmove(void* d, void* s, long l)
{
  unsigned long _d = (unsigned long) d;
  unsigned long _s = (unsigned long) s;
  unsigned long _l = ((l) + SIZEOF_LONG - 1) >> LOG_SIZEOF_LONG;

  if (_d > _s) _my_memmove_d_gt_s(_d, _s, _l);
  else _my_memmove_d_lt_s(_d, _s, _l);
}

#undef memmove
#define memmove(d,s,l) _my_memmove(d, s, l)
#endif

static poly redMora (poly h,sBasisSet::const_iterator end,kStrategy strat);
static poly redBba (poly h,sBasisSet::const_iterator end,kStrategy strat);

#define pDivComp_EQUAL 2
#define pDivComp_LESS 1
#define pDivComp_GREATER -1
#define pDivComp_INCOMP 0
/* Checks the relation of LM(p) and LM(q)
     LM(p) = LM(q) => return pDivComp_EQUAL
     LM(p) | LM(q) => return pDivComp_LESS
     LM(q) | LM(p) => return pDivComp_GREATER
     else return pDivComp_INCOMP */
static inline int pDivCompRing(poly p, poly q)
{
  if ((currRing->pCompIndex < 0)
  || (__p_GetComp(p,currRing) == __p_GetComp(q,currRing)))
  {
    BOOLEAN a=FALSE, b=FALSE;
    int i;
    unsigned long la, lb;
    unsigned long divmask = currRing->divmask;
    for (i=0; i<currRing->VarL_Size; i++)
    {
      la = p->exp[currRing->VarL_Offset[i]];
      lb = q->exp[currRing->VarL_Offset[i]];
      if (la != lb)
      {
        if (la < lb)
        {
          if (b) return pDivComp_INCOMP;
          if (((la & divmask) ^ (lb & divmask)) != ((lb - la) & divmask))
            return pDivComp_INCOMP;
          a = TRUE;
        }
        else
        {
          if (a) return pDivComp_INCOMP;
          if (((la & divmask) ^ (lb & divmask)) != ((la - lb) & divmask))
            return pDivComp_INCOMP;
          b = TRUE;
        }
      }
    }
    if (a) return pDivComp_LESS;
    if (b) return pDivComp_GREATER;
    if (!a & !b) return pDivComp_EQUAL;
  }
  return pDivComp_INCOMP;
}

static inline int pDivComp(poly p, poly q)
{
  if ((currRing->pCompIndex < 0)
  || (__p_GetComp(p,currRing) == __p_GetComp(q,currRing)))
  {
#ifdef HAVE_RATGRING
    if (rIsRatGRing(currRing))
    {
      if (_p_LmDivisibleByPart(p,currRing,
                           q,currRing,
                           currRing->real_var_start, currRing->real_var_end))
        return 0;
      return pLmCmp(q,p); // ONLY FOR GLOBAL ORDER!
    }
#endif
    BOOLEAN a=FALSE, b=FALSE;
    int i;
    unsigned long la, lb;
    unsigned long divmask = currRing->divmask;
    for (i=0; i<currRing->VarL_Size; i++)
    {
      la = p->exp[currRing->VarL_Offset[i]];
      lb = q->exp[currRing->VarL_Offset[i]];
      if (la != lb)
      {
        if (la < lb)
        {
          if (b) return 0;
          if (((la & divmask) ^ (lb & divmask)) != ((lb - la) & divmask))
            return 0;
          a = TRUE;
        }
        else
        {
          if (a) return 0;
          if (((la & divmask) ^ (lb & divmask)) != ((la - lb) & divmask))
            return 0;
          b = TRUE;
        }
      }
    }
    if (a) { /*assume(pLmCmp(q,p)==1);*/ return 1; }
    if (b) { /*assume(pLmCmp(q,p)==-1);*/return -1; }
    /*assume(pLmCmp(q,p)==0);*/
  }
  return 0;
}

#ifdef HAVE_SHIFTBBA
static inline int pLPDivComp(poly p, poly q)
{
  if ((currRing->pCompIndex < 0) || (__p_GetComp(p,currRing) == __p_GetComp(q,currRing)))
  {
    // maybe there is a more performant way to do this? This will get called quite often in bba.
    if (_p_LPLmDivisibleByNoComp(p, q, currRing)) return 1;
    if (_p_LPLmDivisibleByNoComp(q, p, currRing)) return -1;
  }

  return 0;
}
#endif


VAR int     HCord;
VAR int     Kstd1_deg;
VAR int     Kstd1_mu=INT_MAX;

static void deleteHCBucket(LObject *L, kStrategy strat)
{
  if ((strat->kNoether!=NULL)
  && (L->bucket != NULL))
  {
    for (int i=1; i<= (int) L->bucket->buckets_used; i++)
    {
      poly p=L->bucket->buckets[i];
      if(p!=NULL)
      {
        if (p_Cmp(p,strat->kNoetherTail(), L->tailRing) == -1)
        {
          L->bucket->buckets[i]=NULL;
          L->bucket->buckets_length[i]=0;
        }
        else
        {
          do
          {
            if (p_Cmp(pNext(p),strat->kNoetherTail(), L->tailRing) == -1)
            {
              p_Delete(&pNext(p), L->tailRing);
              L->bucket->buckets_length[i]=pLength(L->bucket->buckets[i]);
              break;
            }
            pIter(p);
          } while(p!=NULL);
        }
      }
    }
    int i=L->bucket->buckets_used;
    while ((i>0)&&(L->bucket->buckets[i]==NULL))
    {
      i--;
      L->bucket->buckets_used=i;
    }
  }
}

/*2
*deletes higher monomial of p, re-compute ecart and length
*works only for orderings with ecart =pFDeg(end)-pFDeg(start)
*/
void deleteHC(LObject *L, kStrategy strat, BOOLEAN fromNext)
{
  if (strat->kNoether!=NULL)
  {
    kTest_L(L,strat);
    poly p1;
    poly p = L->GetLmTailRing();
    int l = 1;

    if (!fromNext && p_Cmp(p,strat->kNoetherTail(), L->tailRing) == -1)
    {
      if (L->bucket != NULL) kBucketDestroy(&L->bucket);
      L->Delete();
      L->Clear();
      L->ecart = -1;
      return;
    }
    if (L->bucket != NULL)
    {
      deleteHCBucket(L,strat);
      return;
    }
    BOOLEAN cut=FALSE;
    p1 = p;
    while (pNext(p1)!=NULL)
    {
      if (p_LmCmp(pNext(p1), strat->kNoetherTail(), L->tailRing) == -1)
      {
        cut=(pNext(p1)!=NULL);
        if (cut)
        {
          p_Delete(&pNext(p1), L->tailRing);

          if (p1 == p)
          {
            if (L->t_p != NULL)
            {
              assume(L->p != NULL && p == L->t_p);
              pNext(L->p) = NULL;
            }
            L->max_exp  = NULL;
          }
          else if (fromNext)
            L->max_exp  = p_GetMaxExpP(pNext(L->p), L->tailRing ); // p1;
          //if (L->pLength != 0)
          L->pLength = l;
          // Hmmm when called from updateT, then only
          // reset ecart when cut
          if (fromNext)
            L->ecart = L->pLDeg() - L->GetpFDeg();
        }
        break;
      }
      l++;
      pIter(p1);
    }
    if ((!fromNext) && cut)
    {
      L->SetpFDeg();
      L->ecart = L->pLDeg(strat->LDegLast) - L->GetpFDeg();
    }
    kTest_L(L,strat);
  }
}

void deleteHC(poly* p, int* e, int* l,kStrategy strat)
{
  LObject L(*p, currRing, strat->tailRing);

  deleteHC(&L, strat);
  *p = L.p;
  *e = L.ecart;
  *l = L.length;
  if (L.t_p != NULL) p_LmFree(L.t_p, strat->tailRing);
}

/*2
*tests if p.p=monomial*unit and cancels the unit
*/
void cancelunit (LObject* L,BOOLEAN inNF)
{
  if(rHasGlobalOrdering (currRing)) return;
  if(TEST_OPT_CANCELUNIT) return;

  ring r = L->tailRing;
  poly p = L->GetLmTailRing();
  if(p_GetComp(p, r) != 0 && !p_OneComp(p, r)) return;

  number lc=NULL; /*dummy, is always set if rField_is_Ring(r) */
  if (rField_is_Ring(r) /*&& (rHasLocalOrMixedOrdering(r))*/)
    lc = pGetCoeff(p);

  // Leading coef have to be a unit
  // example 2x+4x2 should be simplified to 2x*(1+2x)
  // and 2 is not a unit in Z
  //if ( !(n_IsUnit(pGetCoeff(p), r->cf)) ) return;

  poly h = pNext(p);
  int  i;

  if(rField_is_Ring(currRing))
  {
    loop
    {
      if (h==NULL)
      {
        p_Delete(&pNext(p), r);
        if (!inNF)
        {
          number eins= nCopy(lc);
          if (L->p != NULL)
          {
            pSetCoeff(L->p,eins);
            if (L->t_p != NULL)
              pSetCoeff0(L->t_p,eins);
          }
          else
            pSetCoeff(L->t_p,eins);
          /* p and t_p share the same coeff, if both are !=NULL */
          /* p==NULL==t_p cannot happen here */
        }
        L->ecart = 0;
        L->length = 1;
        //if (L->pLength > 0)
        L->pLength = 1;
        L->max_exp = NULL;

        if (L->t_p != NULL && pNext(L->t_p) != NULL)
          p_Delete(&pNext(L->t_p),r);
        if (L->p != NULL && pNext(L->p) != NULL)
          pNext(L->p) = NULL;
        return;
      }
      i = rVar(r);
      loop
      {
        if (p_GetExp(p,i,r) > p_GetExp(h,i,r)) return; // does not divide
        i--;
        if (i == 0) break; // does divide, try next monom
      }
      //wrp(p); PrintS(" divide ");wrp(h); PrintLn();
      // Note: As long as qring j forbidden if j contains integer (i.e. ground rings are
      //       domains), no zerodivisor test needed  CAUTION
      if (!n_DivBy(pGetCoeff(h),lc,r->cf))
      {
        return;
      }
      pIter(h);
    }
  }
  else
  {
    loop
    {
      if (h==NULL)
      {
        p_Delete(&pNext(p), r);
        if (!inNF)
        {
          number eins=nInit(1);
          if (L->p != NULL)
          {
            pSetCoeff(L->p,eins);
            if (L->t_p != NULL)
              pSetCoeff0(L->t_p,eins);
          }
          else
            pSetCoeff(L->t_p,eins);
          /* p and t_p share the same coeff, if both are !=NULL */
          /* p==NULL==t_p cannot happen here */
        }
        L->ecart = 0;
        L->length = 1;
        //if (L->pLength > 0)
        L->pLength = 1;
        L->max_exp = NULL;

        if (L->t_p != NULL && pNext(L->t_p) != NULL)
          p_Delete(&pNext(L->t_p),r);
        if (L->p != NULL && pNext(L->p) != NULL)
          pNext(L->p) = NULL;

        return;
      }
      i = rVar(r);
      loop
      {
        if (p_GetExp(p,i,r) > p_GetExp(h,i,r)) return; // does not divide
        i--;
        if (i == 0) break; // does divide, try next monom
      }
      //wrp(p); PrintS(" divide ");wrp(h); PrintLn();
      pIter(h);
    }
  }
}

/*2
*pp is the new element in s
*returns TRUE (in strat->kAllAxis) if
*-HEcke is allowed
*-we are in the last componente of the vector
*-on all axis are monomials (all elements in NotUsedAxis are FALSE)
*returns FALSE for pLexOrderings,
*assumes in module case an ordering of type c* !!
* HEckeTest is only called with strat->kAllAxis==FALSE !
*/
void HEckeTest (poly pp,kStrategy strat)
{
  int   j,p;

  if (currRing->pLexOrder
  || rHasMixedOrdering(currRing)
  || (strat->ak >1)
  || (rField_is_Ring(currRing) && (!n_IsUnit(pGetCoeff(pp),currRing->cf))))
  {
    return;
  }
  p=pIsPurePower(pp);
  if (p!=0)
    strat->NotUsedAxis[p] = FALSE;
  /*- the leading term of pp is a power of the p-th variable -*/
  for (j=(currRing->N);j>0; j--)
  {
    if (strat->NotUsedAxis[j])
    {
      strat->kAllAxis=FALSE;
      return;
    }
  }
  strat->kAllAxis=TRUE;
}

/*2
*utilities for TSet, LSet
*/
inline static intset initec (const int maxnr)
{
  return (intset)omAlloc(maxnr*sizeof(int));
}

inline static unsigned long* initsevS (const int maxnr)
{
  return (unsigned long*)omAlloc0(maxnr*sizeof(unsigned long));
}
inline static int* initS_2_R (const int maxnr)
{
  return (int*)omAlloc0(maxnr*sizeof(int));
}

// enlargeT removed — BlockArray grows automatically via ensure_capacity

void cleanT (kStrategy strat)
{
  int i,j;
  poly  p;
  assume(currRing == strat->tailRing || strat->tailRing != NULL);

  // Optional tombstone-pile-up instrumentation for SORDER_STANDARD /
  // SORDER_APPEND tombstone-on-erase (task 503). Set env var
  // SINGULAR_SBASIS_STATS to see the stats at cleanT entry.
  if (getenv("SINGULAR_SBASIS_STATS") != NULL) {
    strat->S.debug_print_stats("cleanT");
  }

  pShallowCopyDeleteProc p_shallow_copy_delete =
    (strat->tailRing != currRing ?
     pGetShallowCopyDeleteProc(strat->tailRing, currRing) :
     NULL);

#ifdef KDEBUG
  // Validate all S entries before cleanT processing
  for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
  {
    if (sit->p != NULL && sit->p->coef == NULL)
    {
      fprintf(stderr, "cleanT ENTRY: S[%d].p=%p already has NULL coef!\n", sit.index(), (void*)sit->p);
    }
  }
#endif

  for (j=0; j < strat->T.size(); j++)
  {
    p = strat->T[j].p;
    strat->T[j].p=NULL;
    if (strat->T[j].max_exp != NULL)
    {
      p_LmFree(strat->T[j].max_exp, strat->tailRing);
    }
    bool found_in_S = false;
    for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
    {
      if (p == sit->p)
      {
        if (strat->T[j].t_p != NULL)
        {
          if (p_shallow_copy_delete!=NULL)
          {
            pNext(p) = p_shallow_copy_delete(pNext(p),strat->tailRing,currRing,
                                           currRing->PolyBin);
          }
          p_LmFree(strat->T[j].t_p, strat->tailRing);
        }
        found_in_S = true;
        break;
      }
    }
    if (!found_in_S)
    {
      if (strat->T[j].t_p != NULL)
      {
        p_Delete(&(strat->T[j].t_p), strat->tailRing);
        p_LmFree(p, currRing);
      }
      else
      {
#ifdef HAVE_SHIFTBBA
        if (currRing->isLPring && strat->T[j].shift > 0)
        {
          pNext(p) = NULL; // pNext(p) points to the unshifted tail, don't try to delete it here
        }
#endif
        pDelete(&p);
      }
    }
#ifdef KDEBUG
    // Check all S entries after processing T[j]
    for (auto skit = strat->S.begin(); skit != strat->S.end(); ++skit)
    {
      if (skit->p != NULL && skit->p->coef == NULL)
      {
        fprintf(stderr, "cleanT: after T[%d] (p=%p, t_p=%p, matched S=%s), S[%d].p=%p got NULL coef!\n",
                j, (void*)p, (void*)strat->T[j].t_p,
                found_in_S ? "found" : "not_found",
                skit.index(), (void*)skit->p);
        abort();
      }
    }
#endif
  }
  strat->T.setsize(0);
}

void cleanTSbaRing (kStrategy strat)
{
  int i,j;
  poly  p;
  assume(currRing == strat->tailRing || strat->tailRing != NULL);

  pShallowCopyDeleteProc p_shallow_copy_delete =
    (strat->tailRing != currRing ?
     pGetShallowCopyDeleteProc(strat->tailRing, currRing) :
     NULL);
  for (j=0; j < strat->T.size(); j++)
  {
    p = strat->T[j].p;
    strat->T[j].p=NULL;
    if (strat->T[j].max_exp != NULL)
    {
      p_LmFree(strat->T[j].max_exp, strat->tailRing);
    }
    bool found_in_S = false;
    for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
    {
      if (p == sit->p)
      {
        if (strat->T[j].t_p != NULL)
        {
          assume(p_shallow_copy_delete != NULL);
          pNext(p) = p_shallow_copy_delete(pNext(p),strat->tailRing,currRing,
                                           currRing->PolyBin);
          p_LmFree(strat->T[j].t_p, strat->tailRing);
        }
        found_in_S = true;
        break;
      }
    }
    if (!found_in_S)
    {
      if (strat->T[j].t_p != NULL)
      {
        p_Delete(&(strat->T[j].t_p), strat->tailRing);
        p_LmFree(p, currRing);
      }
      else
      {
        //pDelete(&p);
        p = NULL;
      }
    }
  }
  strat->T.setsize(0);
}

// initPairtest is no longer needed — pairtest is per-SElement.
// clear_pairtest() on sBasisSet resets all flags.


/*2
*test whether (p1,p2) or (p2,p1) is in L at or after iterator it
*it returns TRUE if yes and modifies the iterator to point to the match
*Uses pair_index for O(1) lookup, then checks position constraint
*/
BOOLEAN isInPairsetL(LSet::iterator &it,poly p1,poly p2,kStrategy strat)
{
  if (it == strat->L.end()) return FALSE;
  if (p1 == NULL || p2 == NULL) return FALSE;
  auto key = LSet::canonicalize_pair(p1, p2);
  auto found = strat->L.pair_index.find(key);
  if (found != strat->L.pair_index.end()) {
    LSet::iterator candidate = found->second;
    // Check position constraint: candidate must be at or after 'it'
    if (candidate == it || !strat->L.key_comp()(*candidate, *it)) {
      it = candidate;
      return TRUE;
    }
  }
  return FALSE;
}

int kFindInT(poly p, const BlockArray<TObject> &T, int tlength)
{
  int i;

  for (i=0; i<=tlength; i++)
  {
    if (T[i].p == p) return i;
  }
  return -1;
}

int kFindInT(poly p, kStrategy strat)
{
  int i;
  do
  {
    i = kFindInT(p, strat->T, strat->T.size()-1);
    if (i >= 0) return i;
    strat = strat->next;
  }
  while (strat != NULL);
  return -1;
}

#ifdef HAVE_SHIFTBBA
int kFindInTShift(poly p, const BlockArray<TObject> &T, int tlength)
{
  int i;

  for (i=0; i<=tlength; i++)
  {
    // in the Letterplace ring the LMs in T and L are copies thus we have to use pEqualPolys() instead of ==
    if (pEqualPolys(T[i].p, p)) return i;
  }
  return -1;
}
#endif

#ifdef HAVE_SHIFTBBA
int kFindInTShift(poly p, kStrategy strat)
{
  int i;
  do
  {
    i = kFindInTShift(p, strat->T, strat->T.size()-1);
    if (i >= 0) return i;
    strat = strat->next;
  }
  while (strat != NULL);
  return -1;
}
#endif

#ifdef KDEBUG
void sTObject::wrp()
{
  if (t_p != NULL) p_wrp(t_p, tailRing);
  else if (p != NULL) p_wrp(p, currRing, tailRing);
  else ::wrp(NULL);
}
#endif

#define kFalseReturn(x) do { if (!x) return FALSE;} while (0)

#ifdef KDEBUG
// check that Lm's of a poly from T are "equal"
static const char* kTest_LmEqual(poly p, poly t_p, ring tailRing)
{
  int i;
  for (i=1; i<=tailRing->N; i++)
  {
    if (p_GetExp(p, i, currRing) != p_GetExp(t_p, i, tailRing))
      return "Lm[i] different";
  }
  if (p_GetComp(p, currRing) != p_GetComp(t_p, tailRing))
    return "Lm[0] different";
  if (pNext(p) != pNext(t_p))
    return "Lm.next different";
  if (pGetCoeff(p) != pGetCoeff(t_p))
    return "Lm.coeff different";
  return NULL;
}
#endif

#ifdef KDEBUG
STATIC_VAR BOOLEAN sloppy_max = FALSE;
BOOLEAN kTest_T(TObject * T, kStrategy strat, int i, char TN)
{
  ring tailRing = T->tailRing;
  ring strat_tailRing = strat->tailRing;
  if (strat_tailRing == NULL) strat_tailRing = tailRing;
  r_assume(strat_tailRing == tailRing);

  poly p = T->p;
  // ring r = currRing;

  if (T->p == NULL && T->t_p == NULL && i >= 0)
    return dReportError("%c[%d].poly is NULL", TN, i);

  if (T->p!=NULL)
  {
    nTest(pGetCoeff(T->p));
    if ((T->t_p==NULL)&&(pNext(T->p)!=NULL)) p_Test(pNext(T->p),currRing);
  }
  if (T->t_p!=NULL)
  {
    nTest(pGetCoeff(T->t_p));
    if (pNext(T->t_p)!=NULL) p_Test(pNext(T->t_p),strat_tailRing);
  }
  if ((T->p!=NULL)&&(T->t_p!=NULL)) assume(pGetCoeff(T->p)==pGetCoeff(T->t_p));

  if (T->tailRing != currRing)
  {
    if (T->t_p == NULL && i > 0)
      return dReportError("%c[%d].t_p is NULL", TN, i);
    pFalseReturn(p_Test(T->t_p, T->tailRing));
    if (T->p != NULL) pFalseReturn(p_LmTest(T->p, currRing));
    if ((T->p != NULL) && (T->t_p != NULL))
    {
      const char* msg = kTest_LmEqual(T->p, T->t_p, T->tailRing);
      if (msg != NULL)
        return dReportError("%c[%d] %s", TN, i, msg);
      // r = T->tailRing;
      p = T->t_p;
    }
    if (T->p == NULL)
    {
      p = T->t_p;
      // r = T->tailRing;
    }
    if (T->t_p != NULL && i >= 0 && TN == 'T')
    {
      if (pNext(T->t_p) == NULL)
      {
        if (T->max_exp != NULL)
          return dReportError("%c[%d].max_exp is not NULL as it should be", TN, i);
      }
      else
      {
        if (T->max_exp == NULL)
          return dReportError("%c[%d].max_exp is NULL", TN, i);
        if (pNext(T->max_exp) != NULL)
          return dReportError("pNext(%c[%d].max_exp) != NULL", TN, i);

        pFalseReturn(p_CheckPolyRing(T->max_exp, tailRing));
        omCheckBinAddrSize(T->max_exp, (omSizeWOfBin(tailRing->PolyBin))*SIZEOF_LONG);
#if KDEBUG > 0
        if (! sloppy_max)
        {
          poly test_max = p_GetMaxExpP(pNext(T->t_p), tailRing);
          p_Setm(T->max_exp, tailRing);
          p_Setm(test_max, tailRing);
          BOOLEAN equal = p_ExpVectorEqual(T->max_exp, test_max, tailRing);
          if (! equal)
            return dReportError("%c[%d].max out of sync", TN, i);
          p_LmFree(test_max, tailRing);
        }
#endif
      }
    }
  }
  else
  {
    if (T->p == NULL && i > 0)
      return dReportError("%c[%d].p is NULL", TN, i);
#ifdef HAVE_SHIFTBBA
    if (currRing->isLPring && T->shift > 0)
    {
      // in this case, the order is not correct. test LM and tail separately
      pFalseReturn(p_LmTest(T->p, currRing));
      pFalseReturn(p_Test(pNext(T->p), currRing));
    }
    else
#endif
    {
      pFalseReturn(p_Test(T->p, currRing));
    }
  }

  if ((i >= 0) && (T->pLength != 0)
  && (! rIsSyzIndexRing(currRing)) && (T->pLength != pLength(p)))
  {
    int l=T->pLength;
    T->pLength=pLength(p);
    return dReportError("%c[%d] pLength error: has %d, specified to have %d",
                        TN, i , pLength(p), l);
  }

  // check FDeg,  for elements in L and T
  if (i >= 0 && (TN == 'T' || TN == 'L'))
  {
    // FDeg has ir element from T of L set
    if (strat->homog && (T->FDeg  != T->pFDeg()))
    {
      int d=T->FDeg;
      T->FDeg=T->pFDeg();
      return dReportError("%c[%d] FDeg error: has %d, specified to have %d",
                          TN, i , T->pFDeg(), d);
    }
  }

  // check is_normalized for elements in T
  if (i >= 0 && TN == 'T')
  {
    if (T->is_normalized && ! nIsOne(pGetCoeff(p)))
      return dReportError("T[%d] is_normalized error", i);

  }
  return TRUE;
}
#endif

#ifdef KDEBUG
BOOLEAN kTest_L(LObject *L, kStrategy strat,
                BOOLEAN testp, int lpos, BlockArray<TObject> *T, int tlength)
{
  ring strat_tailRing=strat->tailRing;
  if (L->p!=NULL)
  {
    if ((L->t_p==NULL)
    &&(pNext(L->p)!=NULL)
    &&(pGetCoeff(pNext(L->p))!=NULL)) /* !=strat->tail*/
    {
      p_Test(pNext(L->p),currRing);
      nTest(pGetCoeff(L->p));
    }
  }
  if (L->t_p!=NULL)
  {
    if ((pNext(L->t_p)!=NULL)
    &&(pGetCoeff(pNext(L->t_p))!=NULL)) /* !=strat->tail*/
    {
      p_Test(pNext(L->t_p),strat_tailRing);
      nTest(pGetCoeff(L->t_p));
    }
  }
  if ((L->p!=NULL)&&(L->t_p!=NULL)) assume(pGetCoeff(L->p)==pGetCoeff(L->t_p));

  if (testp)
  {
    if (L->bucket != NULL)
    {
      kFalseReturn(kbTest(L->bucket));
      r_assume(L->bucket->bucket_ring == L->tailRing);
      // Previously: save pNext(L->p), set to NULL, call kTest_T,
      // restore.  That transient mutation was not thread-safe — it
      // raced with concurrent readers of L->p in parallel Groebner
      // basis code (ksCreateSpoly captures a2 = pNext(p2) and l2 =
      // T[i_r2].pLength in separate reads; a transient NULL makes
      // them inconsistent, firing kBucketInit's length assertion).
      // After PrepareRed puts the tail in the bucket, pNext(L->p) is
      // already NULL, so this dance was redundant for parallel-bba
      // callers.  If some caller reaches here with pNext(L->p) !=
      // NULL we'd rather catch that than hide it with a mutation.
    }
    kFalseReturn(kTest_T(L, strat, lpos, 'L'));

    ring r;
    poly p;
    L->GetLm(p, r);
    if (L->sev != 0L)
    {
      if (p_GetShortExpVector(p, r) != L->sev)
      {
        return dReportError("L[%d] wrong sev: has %lo, specified to have %lo",
                          lpos, p_GetShortExpVector(p, r), L->sev);
      }
    }
  }
  if (L->p1 == NULL)
  {
    // L->p2 either NULL or "normal" poly
    pFalseReturn(pp_Test(L->p2, currRing, L->tailRing));
  }
  else if (tlength > 0 && T != NULL && (lpos >=0))
  {
    // now p1 and p2 must be != NULL and must be contained in T
    int i;
#ifdef HAVE_SHIFTBBA
    if (rIsLPRing(currRing))
      i = kFindInTShift(L->p1, *T, tlength);
    else
#endif
      i = kFindInT(L->p1, *T, tlength);
    if (i < 0)
      return dReportError("L[%d].p1 not in T",lpos);
#ifdef HAVE_SHIFTBBA
    if (rIsLPRing(currRing))
    {
      if (rField_is_Ring(currRing)) return TRUE; // m*shift(q) is not in T
      i = kFindInTShift(L->p2, *T, tlength);
    }
    else
#endif
      i = kFindInT(L->p2, *T, tlength);
    if (i < 0)
      return dReportError("L[%d].p2 not in T",lpos);
  }
  return TRUE;
}
#endif

#ifdef KDEBUG
BOOLEAN kTest (kStrategy strat)
{
  int i;
  // test P
  kFalseReturn(kTest_L(&(strat->P), strat,
                       (strat->P.p != NULL && pNext(strat->P.p)!=strat->tail),
                       -1, &strat->T, strat->T.size()-1));

  // test T
  if ((!strat->T.empty()))
  {
    for (i=0; i < strat->T.size(); i++)
    {
      kFalseReturn(kTest_T(&(strat->T[i]), strat, i, 'T'));
      if (strat->sevT[i] != pGetShortExpVector(strat->T[i].p))
        return dReportError("strat->sevT[%d] out of sync", i);
    }
  }

  // test L
  i=0;
  for (auto& Lp: strat->L) {
    kFalseReturn(kTest_L(&Lp, strat,
                         Lp.Next() != strat->tail, i,
                         &strat->T, strat->T.size()-1));
    i++;
      // may be unused
      //if (strat->use_buckets && Lp.Next() != strat->tail &&
      //    Lp.Next() != NULL && Lp.p1 != NULL)
      //{
      //  assume(strat->L[i].bucket != NULL);
      //}
  }

  // test S
  if (!strat->S.empty())
    kFalseReturn(kTest_S(strat));

  return TRUE;
}
#endif

#ifdef KDEBUG
BOOLEAN kTest_S(kStrategy strat)
{
  BOOLEAN ret = TRUE;
  for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
  {
    if (sit->p != NULL &&
        sit->sev != pGetShortExpVector(sit->p))
    {
      return dReportError("S[%d] wrong sev: has %o, specified to have %o",
                          sit.index() , pGetShortExpVector(sit->p), sit->sev);
    }
  }
  return ret;
}
#endif

#ifdef KDEBUG
BOOLEAN kTest_TS(kStrategy strat)
{
  int i, j;
  // BOOLEAN ret = TRUE;
  kFalseReturn(kTest(strat));

  // test strat->R, strat->T[i].i_r
  for (i=0; i < strat->T.size(); i++)
  {
    if (strat->T[i].i_r < 0 || strat->T[i].i_r > strat->T.size()-1)
      return dReportError("strat->T[%d].i_r == %d out of bounds", i,
                          strat->T[i].i_r);
    if (strat->R[strat->T[i].i_r] != strat->T.addr(i))
      return dReportError("T[%d].i_r with R out of sync", i);
  }
  // test containment of S inT
  if ((!strat->S.empty())&&((!strat->T.empty())))
  {
    for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
    {
      j = kFindInT(sit->p, strat->T, strat->T.size()-1);
      if (j < 0)
        return dReportError("S[%d] not in T", sit.index());
      if (sit->s_2_r != strat->T[j].i_r)
        return dReportError("S_2_R[%d]=%d != T[%d].i_r=%d\n",
                            sit.index(), sit->s_2_r, j, strat->T[j].i_r);
    }
  }
  // test strat->L[i].i_r1
  i=0;
  #ifdef HAVE_SHIFTBBA
  if (!rIsLPRing(currRing)) // in the Letterplace ring we currently don't set/use i_r1 and i_r2
  #endif
  for (auto& Lp: strat->L) {
    if (Lp.p1 != NULL && Lp.p2)
    {
      if (Lp.i_r1 < 0 ||
          Lp.i_r1 > strat->T.size()-1 ||
          Lp.T_1(strat)->p != Lp.p1)
        return dReportError("L[%d].i_r1 out of sync", i);
      if (Lp.i_r2 < 0 ||
          Lp.i_r2 > strat->T.size()-1 ||
          Lp.T_2(strat)->p != Lp.p2)
        return dReportError("L[%d].i_r2 out of sync", i);
    }
    else
    {
      if (Lp.i_r1 != -1)
        return dReportError("L[%d].i_r1 out of sync", i);
      if (Lp.i_r2 != -1)
        return dReportError("L[%d].i_r2 out of sync", i);
    }
    if (Lp.i_r != -1)
      return dReportError("L[%d].i_r out of sync", i);
    i++;
  }
  return TRUE;
}
#endif // KDEBUG

/*2
* sBasisSet::erase_and_next — delete element at iterator, return next valid iterator.
*/
sBasisSet::iterator sBasisSet::erase_and_next(iterator it)
{
  int pos = it.pos_;
  erase(it);
  // In non-lazy mode, erase() shifts elements down, so pos now points
  // to the next element (which slid into position pos). In lazy mode,
  // the element is flagged deleted and we need to advance past it.
  iterator next(this, pos);
  next.skip_deleted_forward();
  return next;
}

/*
* sBasisSet::S_2_T — return TObject corresponding to the S element at the iterator.
* Assumes it exists (no error checking).
*/
TObject* sBasisSet::S_2_T(const_iterator it, kStrategy strat)
{
  int i = it.index();
  assume(i >= 0 && i < size());
  assume(elem(i).s_2_r >= 0 && elem(i).s_2_r < strat->T.size());
  TObject* TT = strat->R[elem(i).s_2_r];
  assume(TT != NULL && TT->p == elem(i).p);
  return TT;
}

/*
* sBasisSet::s_2_t — return TObject corresponding to the S element at the iterator.
* Returns NULL if it cannot be found.
*/
TObject* sBasisSet::s_2_t(const_iterator it, kStrategy strat)
{
  int i = it.index();
  if (i >= 0 && i < size())
  {
    int sri = elem(i).s_2_r;
    if ((sri >= 0) && (sri < strat->T.size()))
    {
      TObject* t = strat->R[sri];
      if ((t != NULL) && (t->p == elem(i).p))
        return t;
    }
    // last but not least, try kFindInT
    sri = kFindInT(elem(i).p, strat->T, strat->T.size() - 1);
    if (sri >= 0)
      return &(strat->T[sri]);
  }
  return NULL;
}

/*
* sBasisSet::clear_if_divisible — if p divides *at, delete *at and adjust the iterator.
* Replaces the clearS inline in kInline.h.
*/
void sBasisSet::clear_if_divisible(poly p, unsigned long p_sev,
                                   iterator &at, kStrategy strat)
{
  assume(p_sev == pGetShortExpVector(p));
  if (strat->noClearS) return;
  int i = at.index();
  if (rField_is_Ring(currRing))
  {
    if (!pLmShortDivisibleBy(p, p_sev, elem(i).p, ~ elem(i).sev))
      return;
    if (!n_DivBy(pGetCoeff(elem(i).p), pGetCoeff(p), currRing->cf))
      return;
  }
  else
  {
    if (!pLmShortDivisibleBy(p, p_sev, elem(i).p, ~ elem(i).sev))
      return;
  }
  // Erase and back up the iterator (like the old clearS did *at-- and *k--)
  at = erase_and_next(at);
  // erase_and_next returns the next valid iterator, but the old clearS
  // decremented *at so the caller's loop would re-increment it.
  // We need to back up one position so the caller's ++at lands on the
  // element that erase_and_next would return.
  --at;
}

sBasisSet::iterator sBasisSet::simple_find_pos(poly p, int len, wlen_type wlen, kStrategy strat)
{
  if (deleted_count_ > 0) compact();
  if (empty()) return begin();
  int length = size() - 1;
  int i, an = 0, en = length;
  if (strat->use_lenSw) {
    wlen_type wl = wlen;
    if ((wl > elem(length).wlength)
        || ((wl == elem(length).wlength) && (pLmCmp(elem(length).p, p) == -1)))
      return end();
    loop {
      if (an >= en - 1) {
        if ((wl < elem(an).wlength)
            || ((wl == elem(an).wlength) && (pLmCmp(elem(an).p, p) == 1)))
          return iterator(this, an);
        return iterator(this, en);
      }
      i = (an + en) / 2;
      if ((wl < elem(i).wlength)
          || ((wl == elem(i).wlength) && (pLmCmp(elem(i).p, p) == 1))) en = i;
      else an = i;
    }
  } else {
    if ((len > elem(length).length)
        || ((len == elem(length).length) && (pLmCmp(elem(length).p, p) == -1)))
      return end();
    loop {
      if (an >= en - 1) {
        if ((len < elem(an).length)
            || ((len == elem(an).length) && (pLmCmp(elem(an).p, p) == 1)))
          return iterator(this, an);
        return iterator(this, en);
      }
      i = (an + en) / 2;
      if ((len < elem(i).length)
          || ((len == elem(i).length) && (pLmCmp(elem(i).p, p) == 1))) en = i;
      else an = i;
    }
  }
}

#ifdef HAVE_SHIFTBBA
static BOOLEAN is_shifted_p1(const kStrategy strat)
{
  if (rIsLPRing(currRing)
  && (strat->P.p1!=NULL))
  {
    // clean up strat->P.p1: may be shifted
    poly p=strat->P.p1;
    int lv=currRing->isLPring;
    BOOLEAN is_shifted=TRUE;
    for (int i=lv;i>0;i--)
    {
      if (pGetExp(p,i)!=0) { is_shifted=FALSE; break;}
    }
    if (is_shifted
    && (! kExistsInL1(p, strat))
    && (kFindInT(p, strat->T, strat->T.size()-1) < 0)
    )
    {
      return TRUE;
    }
  }
  return FALSE;
}
#endif

// Free the polynomials referenced by an LObject (lcm, sig, p) using the
// same cleanup rules as the legacy LSet::erase: lcm via kDeleteLcm, sig
// via pLmDelete/pLmFree depending on coeff, p via pLmDelete/pLmFree when
// attached to strat->tail or via Lp.Delete() when not found in T.
// Shared between erase (called at chainCritNormal / tombstone time) and
// compact (called to flush accumulated tombstones).
// NB: deliberately does not null the poly fields; for tombstone-on-erase
// we defer physical removal from the multiset tree until compact(), and
// the tree comparator must continue to see consistent Lp.p / Lp.lcm /
// Lp.p1 / Lp.p2 pointers for any live entry that happens to be compared
// against this tombstone during subsequent inserts.  compact() is the
// one caller that actually removes the entry from the tree, at which
// point poly pointers become moot.
static void kLSet_free_polys(LObject& Lp, kStrategy strat) {
  if (Lp.lcm != NULL) {
    kDeleteLcm(&Lp);
  }
  if (Lp.sig != NULL) {
    if (pGetCoeff(Lp.sig) != NULL)
      pLmDelete(Lp.sig);
    else
      pLmFree(Lp.sig);
    Lp.sig = NULL;
  }
  if (Lp.p != NULL) {
    if (pNext(Lp.p) == strat->tail) {
      if (pGetCoeff(Lp.p) != NULL)
        pLmDelete(Lp.p);
      else
        pLmFree(Lp.p);
      Lp.p = NULL;
      /*- tail belongs to several int spolys -*/
    } else {
      // search p in T; if it is there, do not delete it
      if (rHasGlobalOrdering(currRing) || (kFindInT(Lp.p, strat) < 0)) {
        Lp.Delete();
      }
    }
  }
  #ifdef HAVE_SHIFTBBA
  /* this logic was added in commit 95b7138 to fix a memory leak */
  if (is_shifted_p1(/*strat->P.p1,*/strat)) {
    // clean up strat->P.p1: may be shifted
    pLmDelete(strat->P.p1);
    strat->P.p1 = NULL;
  }
  #endif
}

// Tombstone an L entry: record the erase, mark the LObject as deleted,
// remove from the dedup pair_index, and zero sev_flat_/sevSig_flat_ so
// the cache-friendly unordered and filtered scans skip the slot.  Does
// NOT touch the multiset tree and does NOT free polys — compact() will
// do both in one pass when the pile-up gets inconvenient.
LSet::iterator LSet::erase(LSet::iterator it) {
  LObject& Lp = *it;
  if (Lp.deleted) {
    // Idempotent: a second erase of the same tombstone is a no-op.
    ++it;
    return it;
  }
  // Remove from pair_index
  if (Lp.p1 != NULL && Lp.p2 != NULL) {
    auto key = canonicalize_pair(Lp.p1, Lp.p2);
    pair_index.erase(key);
  }
  // Mark the sev_flat_ and sevSig_flat_ entries as sentinel (0)
  if (Lp.flat_index < sev_flat_.size()) sev_flat_[Lp.flat_index] = 0;
  if (Lp.flat_index < sevSig_flat_.size()) sevSig_flat_[Lp.flat_index] = 0;
  Lp.deleted = true;
  ++erase_call_count_;
  ++deleted_count_;
  if (deleted_count_ > peak_deleted_count_)
    peak_deleted_count_ = deleted_count_;
  --live_count_;
  // Advance to the next live entry (skip_deleted_forward is triggered
  // by operator++ on the ordered iterator).
  ++it;
  return it;
}

LSet::unordered_iterator LSet::erase(LSet::unordered_iterator it) {
  LObject& Lp = *it;
  if (Lp.deleted) {
    ++it;  // skip_deleted() inside will advance past it
    return it;
  }
  if (Lp.p1 != NULL && Lp.p2 != NULL) {
    auto key = canonicalize_pair(Lp.p1, Lp.p2);
    pair_index.erase(key);
  }
  if (Lp.flat_index < sev_flat_.size()) sev_flat_[Lp.flat_index] = 0;
  if (Lp.flat_index < sevSig_flat_.size()) sevSig_flat_[Lp.flat_index] = 0;
  Lp.deleted = true;
  ++erase_call_count_;
  ++deleted_count_;
  if (deleted_count_ > peak_deleted_count_)
    peak_deleted_count_ = deleted_count_;
  --live_count_;
  ++it;
  return it;
}

LSet::filtered_iterator LSet::erase(LSet::filtered_iterator fit) {
  LObject* lp = flat_ptr(fit.pos_);
  if (lp != nullptr && !lp->deleted) {
    if (lp->p1 != NULL && lp->p2 != NULL) {
      auto key = canonicalize_pair(lp->p1, lp->p2);
      pair_index.erase(key);
    }
    if (fit.pos_ < sev_flat_.size()) sev_flat_[fit.pos_] = 0;
    if (fit.pos_ < sevSig_flat_.size()) sevSig_flat_[fit.pos_] = 0;
    lp->deleted = true;
    ++erase_call_count_;
    ++deleted_count_;
    --live_count_;
  }
  // Build a new filtered_iterator from the advanced position — the
  // filtered_iterator::advance inside the ctor handles tombstone/sev
  // skip-forward.
  return filtered_iterator(this, fit.pos_ + 1, fit.sev1_, fit.sev2_, fit.sev_array_);
}

// compact(): physically remove all tombstoned entries from the multiset
// tree, free their polys, and rebuild sev_flat_/sevSig_flat_/pair_index
// / flat_index against the compacted tree.  Invalidates ALL outstanding
// iterators.  Fast path: if deleted_count_ == 0, update peak and return.
void LSet::compact() {
  ++compact_call_count_;
  if (deleted_count_ > peak_deleted_count_)
    peak_deleted_count_ = deleted_count_;
  if (deleted_count_ == 0) return;

  const kStrategy strat = key_comp().strat;

  // Walk the flat_ array (physical order).  For each tombstoned entry,
  // free polys and remove its underlying multiset iterator (dealloc
  // the LObject).  Use writable_set<>::erase on each to go through the
  // base-class cleanup (flat_[fi] = data_.end(); delete *it; data_.erase).
  //
  // We collect the unordered_iterator positions of tombstoned entries
  // first, then erase them.  This avoids iterator-invalidation issues
  // while walking flat_ (erase() does flat_[fi] = end() at the iteration
  // index — safe — but for clarity we batch).
  size_t n = flat_size();
  std::vector<size_t> tomb_positions;
  tomb_positions.reserve(deleted_count_);
  for (size_t i = 0; i < n; i++) {
    LObject* lp = flat_ptr(i);
    if (lp != nullptr && lp->deleted) {
      tomb_positions.push_back(i);
    }
  }
  for (size_t pos : tomb_positions) {
    // flat_ptr may have been cleared by a previous erase (if two
    // tombstones happened to share a flat index, which shouldn't
    // happen, but be defensive).
    LObject* lp = flat_ptr(pos);
    if (lp == nullptr) continue;
    kLSet_free_polys(*lp, strat);
    // Use erase_at (raw, no skip) — uiter_at + erase(unordered) would
    // first call skip_deleted which advances past our tombstone to a
    // live entry, and we'd erase the wrong slot.
    writable_set<LObject, CompareLObject>::erase_at(pos);
  }

  // Rebuild flat_index for remaining entries.  writable_set stores
  // multiset iterators in flat_; we need to rebuild sev_flat_ /
  // sevSig_flat_ / pair_index.  After the erase above, flat_[pos] is
  // data_.end() (i.e. flat_ptr(pos) returns nullptr) for tombstones.
  // We repack flat_ so only live entries remain, updating each live
  // LObject's flat_index to its new position.
  //
  // We do this by calling reorder-ish code directly.  writable_set
  // exposes a reorder() that rebuilds flat_ from data_; reuse that.
  // But reorder() also invalidates comparator-ordered positions — fine
  // here because no tree keys changed, we just need the flat_index
  // repacked.  writable_set::reorder does: data_.swap(old) then
  // re-insert pointers; flat_ cleared and rebuilt.  That's overkill
  // since we don't need to re-sort the tree — but it's correct and
  // cheap (tree has log-n inserts).
  //
  // Alternative: implement a tighter repack.  For clarity/correctness
  // we just delegate to the base reorder.  Note: that is safe because
  // CompareLObject is deterministic given the same (p, lcm, p1, p2,
  // ecart, length, seq) — pointers in tree still compare the same
  // way after reorder.
  writable_set<LObject, CompareLObject>::reorder();
  deleted_count_ = 0;

  // Rebuild sev_flat_, sevSig_flat_, pair_index against the compacted
  // flat_ layout.  Live entries now sit at contiguous positions
  // 0..physical_size()-1 with flat_index set correctly.
  rebuild_sev_flat();
  rebuild_sevSig_flat();
  rebuild_pair_index();
}

/*2
* computes the normal ecart;
* used in mora case and if pLexOrder & sugar in bba case
*/
void initEcartNormal (TObject* h)
{
  h->FDeg = h->pFDeg();
  h->ecart = h->pLDeg() - h->FDeg;
  // h->length is set by h->pLDeg
  h->length=h->pLength=pLength(h->p);
}

void initEcartBBA (TObject* h)
{
  h->FDeg = h->pFDeg();
  (*h).ecart = 0;
  h->length=h->pLength=pLength(h->p);
}

void initEcartPairBba (LObject* Lp,poly /*f*/,poly /*g*/,int /*ecartF*/,int /*ecartG*/)
{
  Lp->FDeg = Lp->pFDeg();
  (*Lp).ecart = 0;
  (*Lp).length = 0;
}

void initEcartPairMora (LObject* Lp,poly /*f*/,poly /*g*/,int ecartF,int ecartG)
{
  Lp->FDeg = Lp->pFDeg();
  (*Lp).ecart = si_max(ecartF,ecartG);
  (*Lp).ecart = (*Lp).ecart- (Lp->FDeg -p_FDeg((*Lp).lcm,currRing));
  (*Lp).length = 0;
}

/*2
*if ecart1<=ecart2 it returns TRUE
*/
static inline BOOLEAN sugarDivisibleBy(int ecart1, int ecart2)
{
  return (ecart1 <= ecart2);
}

/*2
* put the pair (s[i],p)  into the set B, ecart=ecart(p) (ring case)
*/
static void enterOnePairRing (const SElement &si,poly p,int /*ecart*/, int isFromQ,kStrategy strat, int atR)
{
  assume(atR >= 0);
  assume(p!=NULL);
  assume(rField_is_Ring(currRing));
  #if ALL_VS_JUST
  //Over rings, if we construct the strong pair, do not add the spair
  if(rField_is_Ring(currRing))
  {
    number s,t,d;
    d = n_ExtGcd(pGetCoeff(p), pGetCoeff(si.p), &s, &t, currRing->cf);

    if (!nIsZero(s) && !nIsZero(t))  // evtl. durch divBy tests ersetzen
    {
      nDelete(&d);
      nDelete(&s);
      nDelete(&t);
      return;
    }
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
  }
  #endif
  int      compare,compareCoeff;
  LObject  h;

#ifdef KDEBUG
  h.ecart=0; h.length=0;
#endif
  /*- computes the lcm(s[i],p) -*/
  if(pHasNotCFRing(p,si.p))
  {
      strat->cp++;
      return;
  }
  h.lcm = p_Lcm(p,si.p,currRing);
  h.sev_lcm = p_GetShortExpVector(h.lcm, currRing);
  pSetCoeff0(h.lcm, n_Lcm(pGetCoeff(p), pGetCoeff(si.p), currRing->cf));
  if (nIsZero(pGetCoeff(h.lcm)))
  {
      strat->cp++;
      pLmDelete(h.lcm);
      return;
  }
  // basic chain criterion
  /*
  *the set B collects the pairs of type (S[j],p)
  *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p) != lcm(r,p)
  *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
  *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
  */

  for(auto jt = strat_B(strat).ubegin(); jt != strat_B(strat).uend(); )
  {
    bool j_deleted = false;
    compare=pDivCompRing(jt->lcm,h.lcm);
    compareCoeff = n_DivComp(pGetCoeff(jt->lcm), pGetCoeff(h.lcm), currRing->cf);
    if(compare == pDivComp_EQUAL)
    {
      //They have the same LM
      if(compareCoeff == pDivComp_LESS)
      {
        if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
        {
          strat->c3++;
          pLmDelete(h.lcm);
          return;
        }
        break;
      }
      if(compareCoeff == pDivComp_GREATER)
      {
        jt = strat_B(strat).erase(jt);
        j_deleted = true;
        strat->c3++;
      }
      if(compareCoeff == pDivComp_EQUAL)
      {
        if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
        {
          strat->c3++;
          pLmDelete(h.lcm);
          return;
        }
        break;
      }
    }
    if(compareCoeff == compare || compareCoeff == pDivComp_EQUAL)
    {
      if(compare == pDivComp_LESS)
      {
        if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
        {
          strat->c3++;
          pLmDelete(h.lcm);
          return;
        }
        break;
      }
      if(compare == pDivComp_GREATER)
      {
        jt = strat_B(strat).erase(jt);
        j_deleted = true;
        strat->c3++;
      }
    }
    if (! j_deleted) ++jt;
  }
  number s, t;
  poly m1, m2, gcd = NULL;
  s = pGetCoeff(si.p);
  t = pGetCoeff(p);
  k_GetLeadTerms(p,si.p,currRing,m1,m2,currRing);
  ksCheckCoeff(&s, &t, currRing->cf);
  pSetCoeff0(m1, s);
  pSetCoeff0(m2, t);
  m2 = pNeg(m2);
  p_Test(m1,strat->tailRing);
  p_Test(m2,strat->tailRing);
  poly si_copy = pCopy(si.p);
  poly pm1 = pp_Mult_mm(pNext(p), m1, strat->tailRing);
  poly sim2 = pp_Mult_mm(pNext(si_copy), m2, strat->tailRing);
  pDelete(&si_copy);
  p_LmDelete(m1, currRing);
  p_LmDelete(m2, currRing);
  if(sim2 == NULL)
  {
    if(pm1 == NULL)
    {
      if(h.lcm != NULL)
      {
        pLmDelete(h.lcm);
        h.lcm=NULL;
      }
      h.Clear();
      record_pairtest_hit(si, strat);
      return;
    }
    else
    {
      gcd = pm1;
      pm1 = NULL;
    }
  }
  else
  {
    if((pGetComp(si.p) == 0) && (0 != pGetComp(p)))
    {
      p_SetCompP(sim2, pGetComp(p), strat->tailRing);
      pSetmComp(sim2);
    }
    //p_Write(pm1,strat->tailRing);p_Write(sim2,strat->tailRing);
    gcd = p_Add_q(pm1, sim2, strat->tailRing);
  }
  p_Test(gcd, strat->tailRing);
#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    wrp(gcd);
    PrintLn();
  }
#endif
  h.p = gcd;
  h.i_r = -1;
  if(h.p == NULL)
  {
    record_pairtest_hit(si, strat);
    return;
  }
  h.tailRing = strat->tailRing;
  //h.pCleardenom();
  //pSetm(h.p);
  h.i_r1 = -1;h.i_r2 = -1;
  strat->initEcart(&h);
  #if 1
  h.p2 = si.p;
  h.p1 = p;
  #endif
  #if 1
  if (atR >= 0)
  {
    h.i_r1 = atR;
    h.i_r2 = si.s_2_r;
  }
  #endif
  h.sev = pGetShortExpVector(h.p);
  if (currRing!=strat->tailRing)
    h.t_p = k_LmInit_currRing_2_tailRing(h.p, strat->tailRing);
  if (strat->P.p!=NULL) strat->P.sev = pGetShortExpVector(strat->P.p);
  else strat->P.sev=0L;
  strat_B(strat).push(h);
  kTest_TS(strat);
}

/*2
* put the  lcm(s[i],p)  into the set B
*/

static BOOLEAN enterOneStrongPoly (const SElement &si,poly p,int /*ecart*/, int /*isFromQ*/,kStrategy strat, int atR, bool enterTstrong)
{
  number d, s, t;
  assume(atR >= 0);
  assume(rField_is_Ring(currRing));
  poly m1, m2, gcd,si_p;
  si_p = si.p;
  //printf("\n--------------------------------\n");
  //pWrite(p);pWrite(si_p);
  d = n_ExtGcd(pGetCoeff(p), pGetCoeff(si_p), &s, &t, currRing->cf);

  if (nIsZero(s) || nIsZero(t))  // evtl. durch divBy tests ersetzen
  {
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
    return FALSE;
  }

  k_GetStrongLeadTerms(p, si_p, currRing, m1, m2, gcd, strat->tailRing);

  if (!rHasLocalOrMixedOrdering(currRing))
  {
    unsigned long sev = pGetShortExpVector(gcd);

    for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
    {
      if (&*sjt == &si)
        continue;

      if (n_DivBy(d, pGetCoeff(sjt->p), currRing->cf)
      && !(sjt->sev & ~sev)
      && p_LmDivisibleBy(sjt->p, gcd, currRing))
      {
        nDelete(&d);
        nDelete(&s);
        nDelete(&t);
        return FALSE;
      }
    }
  }

  //p_Test(m1,strat->tailRing);
  //p_Test(m2,strat->tailRing);
  /*if(!enterTstrong)
  {
    while (! kCheckStrongCreation(atR, m1, i, m2, strat) )
    {
      memset(&(strat->P), 0, sizeof(strat->P));
      kStratChangeTailRing(strat);
      strat->P = *(strat->R[atR]);
      p_LmFree(m1, strat->tailRing);
      p_LmFree(m2, strat->tailRing);
      p_LmFree(gcd, currRing);
      k_GetStrongLeadTerms(p, si, currRing, m1, m2, gcd, strat->tailRing);
    }
  }*/
  pSetCoeff0(m1, s);
  pSetCoeff0(m2, t);
  pSetCoeff0(gcd, d);
  p_Test(m1,strat->tailRing);
  p_Test(m2,strat->tailRing);
  //printf("\n===================================\n");
  //pWrite(m1);pWrite(m2);pWrite(gcd);
#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    // Print("t = %d; s = %d; d = %d\n", nInt(t), nInt(s), nInt(d));
    PrintS("m1 = ");
    p_wrp(m1, strat->tailRing);
    PrintS(" ; m2 = ");
    p_wrp(m2, strat->tailRing);
    PrintS(" ; gcd = ");
    wrp(gcd);
    PrintS("\n--- create strong gcd poly: ");
    PrintS("\n p: ");
    wrp(p);
    PrintS("\n strat->S[].p: ");
    wrp(si_p);
    PrintS(" ---> ");
  }
#endif

  pNext(gcd) = p_Add_q(pp_Mult_mm(pNext(p), m1, strat->tailRing), pp_Mult_mm(pNext(si_p), m2, strat->tailRing), strat->tailRing);
  p_LmDelete(m1, strat->tailRing);
  p_LmDelete(m2, strat->tailRing);
#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    wrp(gcd);
    PrintLn();
  }
#endif

  LObject h;
  h.p = gcd;
  h.tailRing = strat->tailRing;
  strat->initEcart(&h);
  h.sev = pGetShortExpVector(h.p);
  h.i_r1 = -1;h.i_r2 = -1;
  if (currRing!=strat->tailRing)
    h.t_p = k_LmInit_currRing_2_tailRing(h.p, strat->tailRing);
  if(!enterTstrong)
  {
    #if 1
    h.p1 = p;h.p2 = si.p;
    #endif
    if (atR >= 0)
    {
      h.i_r2 = si.s_2_r;
      h.i_r1 = atR;
    }
    else
    {
      h.i_r1 = -1;
      h.i_r2 = -1;
    }
    strat->L.push(h);
  }
  else
  {
    if(h.IsNull()) return FALSE;
    //int red_result;
    //reduzieren ist teur!!!
    //if(strat->L != NULL)
      //red_result = strat->red(&h,strat);
    if(!h.IsNull())
    {
      enterT(h, strat,-1);
      //auto pos = strat->S.find_pos(h.p,h.ecart);
      //strat->enterS(h, strat, strat->T.size()-1, strat->S.end());
    }
  }
  return TRUE;
}

BOOLEAN sbaCheckGcdPair (LObject* h,kStrategy strat)
{
  if(strat->S.empty()) return FALSE;
  // iterate up to S.size()-1 (exclude last element)
  auto send = strat->S.end();
  if (send != strat->S.begin()) --send;
  for(auto sit = strat->S.begin(); sit != send; ++sit)
  {
    //Construct the gcd pair between h and S[i]
    number d, s, t;
    poly m1, m2, gcd;
    d = n_ExtGcd(pGetCoeff(h->p), pGetCoeff(sit->p), &s, &t, currRing->cf);
    if (nIsZero(s) || nIsZero(t))  // evtl. durch divBy tests ersetzen
    {
      nDelete(&d);
      nDelete(&s);
      nDelete(&t);
    }
    else
    {
      k_GetStrongLeadTerms(h->p, sit->p, currRing, m1, m2, gcd, strat->tailRing);
      pSetCoeff0(m1, s);
      pSetCoeff0(m2, t);
      pSetCoeff0(gcd, d);
      pNext(gcd) = p_Add_q(pp_Mult_mm(pNext(h->p), m1, strat->tailRing), pp_Mult_mm(pNext(sit->p), m2, strat->tailRing), strat->tailRing);
      poly pSigMult = p_Copy(h->sig,currRing);
      poly sSigMult = p_Copy(sit->sig,currRing);
      pSigMult = p_Mult_mm(pSigMult,m1,currRing);
      sSigMult = p_Mult_mm(sSigMult,m2,currRing);
      p_LmDelete(m1, strat->tailRing);
      p_LmDelete(m2, strat->tailRing);
      poly pairsig = p_Add_q(pSigMult,sSigMult,currRing);
      if(pairsig!= NULL && pLtCmp(pairsig,h->sig) == 0)
      {
        pDelete(&h->p);
        h->p = gcd;
        pDelete(&h->sig);
        h->sig = pairsig;
        pNext(h->sig) = NULL;
        strat->initEcart(h);
        h->sev = pGetShortExpVector(h->p);
        h->sevSig = pGetShortExpVector(h->sig);
        h->i_r1 = -1;h->i_r2 = -1;
        if(h->lcm != NULL)
        {
          pLmDelete(h->lcm);
          h->lcm = NULL;
        }
        if (currRing!=strat->tailRing)
          h->t_p = k_LmInit_currRing_2_tailRing(h->p, strat->tailRing);
        return TRUE;
      }
      //Delete what you didn't use
      pDelete(&gcd);
      pDelete(&pairsig);
    }
  }
  return FALSE;
}

static BOOLEAN enterOneStrongPolySig (const SElement &si_elem,poly p,poly sig,int /*ecart*/, int /*isFromQ*/,kStrategy strat, int atR)
{
  number d, s, t;
  assume(atR >= 0);
  poly m1, m2, gcd,si;
  si = si_elem.p;
  //printf("\n--------------------------------\n");
  //pWrite(p);pWrite(si);
  d = n_ExtGcd(pGetCoeff(p), pGetCoeff(si), &s, &t, currRing->cf);

  if (nIsZero(s) || nIsZero(t))  // evtl. durch divBy tests ersetzen
  {
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
    return FALSE;
  }

  k_GetStrongLeadTerms(p, si, currRing, m1, m2, gcd, strat->tailRing);
  //p_Test(m1,strat->tailRing);
  //p_Test(m2,strat->tailRing);
  /*if(!enterTstrong)
  {
    while (! kCheckStrongCreation(atR, m1, i, m2, strat) )
    {
      memset(&(strat->P), 0, sizeof(strat->P));
      kStratChangeTailRing(strat);
      strat->P = *(strat->R[atR]);
      p_LmFree(m1, strat->tailRing);
      p_LmFree(m2, strat->tailRing);
      p_LmFree(gcd, currRing);
      k_GetStrongLeadTerms(p, si, currRing, m1, m2, gcd, strat->tailRing);
    }
  }*/
  pSetCoeff0(m1, s);
  pSetCoeff0(m2, t);
  pSetCoeff0(gcd, d);
  p_Test(m1,strat->tailRing);
  p_Test(m2,strat->tailRing);
  //printf("\n===================================\n");
  //pWrite(m1);pWrite(m2);pWrite(gcd);
#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    // Print("t = %d; s = %d; d = %d\n", nInt(t), nInt(s), nInt(d));
    PrintS("m1 = ");
    p_wrp(m1, strat->tailRing);
    PrintS(" ; m2 = ");
    p_wrp(m2, strat->tailRing);
    PrintS(" ; gcd = ");
    wrp(gcd);
    PrintS("\n--- create strong gcd poly: ");
    PrintS("\n p: ");
    wrp(p);
    PrintS("\n strat->S[].p: ");
    wrp(si);
    PrintS(" ---> ");
  }
#endif

  pNext(gcd) = p_Add_q(pp_Mult_mm(pNext(p), m1, strat->tailRing), pp_Mult_mm(pNext(si), m2, strat->tailRing), strat->tailRing);

#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    wrp(gcd);
    PrintLn();
  }
#endif

  //Check and set the signatures
  poly pSigMult = p_Copy(sig,currRing);
  poly sSigMult = p_Copy(si_elem.sig,currRing);
  pSigMult = p_Mult_mm(pSigMult,m1,currRing);
  sSigMult = p_Mult_mm(sSigMult,m2,currRing);
  p_LmDelete(m1, strat->tailRing);
  p_LmDelete(m2, strat->tailRing);
  poly pairsig;
  if(pLmCmp(pSigMult,sSigMult) == 0)
  {
    //Same lm, have to add them
    pairsig = p_Add_q(pSigMult,sSigMult,currRing);
    //This might be zero
  }
  else
  {
    //Set the sig to either pSigMult or sSigMult
    if(pLtCmp(pSigMult,sSigMult)==1)
    {
      pairsig = pSigMult;
      pDelete(&sSigMult);
    }
    else
    {
      pairsig = sSigMult;
      pDelete(&pSigMult);
    }
  }

  LObject h;
  h.p = gcd;
  h.tailRing = strat->tailRing;
  h.sig = pairsig;
  strat->initEcart(&h);
  h.sev = pGetShortExpVector(h.p);
  h.i_r1 = -1;h.i_r2 = -1;
  if (currRing!=strat->tailRing)
    h.t_p = k_LmInit_currRing_2_tailRing(h.p, strat->tailRing);
  if(h.sig == NULL)
    {
      //sigdrop since we loose the signature
      strat->sigdrop = TRUE;
      //Try to reduce it as far as we can via redRing
      int red_result = redRing(&h,strat);
      if(red_result == 0)
      {
        // Cancel the sigdrop
        p_Delete(&h.sig,currRing);h.sig = NULL;
        strat->sigdrop = FALSE;
        return FALSE;
      }
      else
      {
        strat->enterS(strat->P, strat, strat->T.size(), strat->S.end());
        #if 1
        strat->enterS(h, strat, strat->T.size()-1, strat->S.begin());
        #endif
        return FALSE;
      }
    }
  if(!nGreaterZero(pGetCoeff(h.sig)))
  {
    h.sig = pNeg(h.sig);
    h.p = pNeg(h.p);
  }

    if(rField_is_Ring(currRing) &&  pLtCmp(h.sig,sig) == -1)
    {
      strat->sigdrop = TRUE;
      // Completely reduce it
      int red_result = redRing(&h,strat);
      if(red_result == 0)
      {
        // Reduced to 0
        strat->sigdrop = FALSE;
        p_Delete(&h.sig,currRing);h.sig = NULL;
        return FALSE;
      }
      else
      {
        strat->enterS(strat->P, strat, strat->T.size(), strat->S.end());
        // 0 - add just the original poly causing the sigdrop, 1 - add also this
        #if 1
        strat->enterS(h, strat, strat->T.size(), strat->S.begin());
        #endif
        return FALSE;
      }
    }
  //Check for sigdrop
  if(gcd != NULL && pLtCmp(sig,pairsig) > 0 && pLtCmp(si_elem.sig,pairsig) > 0)
  {
    strat->sigdrop = TRUE;
    //Enter this element to S
    strat->enterS(strat->P, strat, strat->T.size(), strat->S.end());
    strat->enterS(h, strat, strat->T.size(), strat->S.end());
  }
  #if 1
  h.p1 = p;h.p2 = si_elem.p;
  #endif
  if (atR >= 0)
  {
    h.i_r2 = si_elem.s_2_r;
    h.i_r1 = atR;
  }
  else
  {
    h.i_r1 = -1;
    h.i_r2 = -1;
  }
  strat->L.push(h);
  return TRUE;
}

/*2
* put the pair (s[i],p)  into the set B, ecart=ecart(p)
*/

// =====================================================================
// SINGULAR_TRACE_PAIRCRIT instrumentation.
//
// enterOnePairNormal is instrumented at every decision point:
//   ENTER           - function entry (after si is bound to p)
//   KILL by=fromT   - fromT ecart kill (sugarCrit branch and else)
//   KILL by=prodCrit- product criterion kill
//   KILL by=domByB  - dominated by existing B entry, pair rejected
//   BKILL           - existing B entry erased by new pair (dominates)
//   ZERO            - spoly(si, p) reduced to 0 at construction
//   KEEP            - pair accepted, pushed into B
// "atT" field = strat->S.size() at the call (monotone per thread).
// LMs are formatted with p_String(...,currRing) and freed.
//
// Callers MUST guard with g_paircrit_this_dispatch to skip formatting
// cost when the trace is disabled.
// =====================================================================

// Helper: format LM (leading monomial only, NOT the full polynomial)
// as a C string.  Caller must omFree() the returned buffer.  Safe
// with NULL.
//
// Writes into a caller-allocated char buffer.  We avoid p_Head /
// p_LmFree because ommalloc/omFree against per-ring bins is not
// guaranteed thread-safe in all parallel contexts (observed SEGV in
// smoke test when called from enterOnePairNormal across worker
// threads).  Instead we manually iterate exponents and format the
// coefficient, producing a stable LM-only string like "3xy2z" or
// "-5" for a constant.  No ring-bin allocation; just omStrDup at the
// end to match the omFree contract the callers already use.
//
// Exported as kt_lm_str for use by the reduce-tracer in kthread.cc
// (and any future tracer that needs a monomial-only LM format
// without touching ring bins).  paircrit_lm_str stays as a local
// alias so the chainCritNormal / enterOnePairNormal call sites below
// don't need to be renamed.
char *kt_lm_str(poly p) {
  if (p == NULL) return NULL;
  const ring r = currRing;
  // 64 bytes per variable is enough for any sane exponent + name.
  // We also reserve room for the coefficient (n_Write produces
  // unbounded strings for multiprecision rationals, so we fall back
  // to a placeholder when the coef is not 1/-1/small-int).  For this
  // diagnostic, a short "cN" tag is sufficient — we're primarily
  // comparing monomials, not coefficients.
  int nvars = rVar(r);
  size_t cap = 64 + (size_t)nvars * 32;
  char *out = (char *)omAlloc(cap);
  size_t off = 0;
  out[0] = 0;
  // Coefficient: for Q we want readable leading coefficient.  Try a
  // minimal formatter: long integer / rational check via n_IsOne /
  // n_IsMOne; otherwise print "c" + pointer-ish tag.
  number coef = p_GetCoeff(p, r);
  if (coef != NULL && r->cf != NULL) {
    if (n_IsOne(coef, r->cf)) {
      /* no sign */
    } else if (n_IsMOne(coef, r->cf)) {
      if (off + 1 < cap) out[off++] = '-';
    } else {
      // Generic: mark with "c" placeholder; not core to the diff.
      int n = snprintf(out + off, cap - off, "c?");
      if (n > 0 && (size_t)n < cap - off) off += (size_t)n;
    }
  }
  bool wrote_var = false;
  for (int v = 1; v <= nvars; v++) {
    long e = p_GetExp(p, v, r);
    if (e == 0) continue;
    const char *name = r->names[v - 1];
    if (!name) name = "?";
    int n;
    if (e == 1)
      n = snprintf(out + off, cap - off, "%s", name);
    else
      n = snprintf(out + off, cap - off, "%s%ld", name, e);
    if (n > 0 && (size_t)n < cap - off) off += (size_t)n;
    wrote_var = true;
  }
  if (!wrote_var) {
    // Constant LM — emit "1" if we haven't printed sign/coef.
    if (off == 0 || (off == 1 && out[0] == '-')) {
      int n = snprintf(out + off, cap - off, "1");
      if (n > 0 && (size_t)n < cap - off) off += (size_t)n;
    }
  }
  // Ensure NUL termination.
  if (off >= cap) off = cap - 1;
  out[off] = 0;
  return out;
}

// Local alias so the paircrit call sites (all in this file) don't need
// to be renamed.  static inline lets the compiler fold the call away.
static inline char *paircrit_lm_str(poly p) { return kt_lm_str(p); }

void enterOnePairNormal (const SElement &si,poly p,int ecart, int isFromQ,kStrategy strat, int atR = -1)
{

  int      compare;

  if (g_paircrit_this_dispatch)
  {
    char *h_lm = paircrit_lm_str(p);
    char *s_lm = paircrit_lm_str(si.p);
    kt_paircrit_logf(
      "ENTER tid=%d atT=%d h_lm=%s si_s2r=%d si_arrival=%lu si_lm=%s ecart=%d si_ecart=%d\n",
      kt_debug_tid, (int)strat->S.size(),
      h_lm ? h_lm : "NULL", (int)si.s_2_r, (unsigned long)si.arrival_id,
      s_lm ? s_lm : "NULL", ecart, si.ecart);
    if (h_lm) omFree(h_lm);
    if (s_lm) omFree(s_lm);
  }
  // Task 325 event log: ENTERPAIR.  arg_a = si.arrival_id, arg_c =
  // ecart, pointers = (h_poly, si.p).
  if (g_event_log_enabled) {
    kevlog_emit(EVT_ENTERPAIR, (uint16_t)kt_debug_tid,
                (uint16_t)strat->S.size(), 0,
                (uint32_t)si.arrival_id, 0, (uint32_t)ecart, 0,
                (const void *)p, (const void *)si.p);
  }

  /*- check product criterion and ecart BEFORE computing the lcm -*/
  if (strat->sugarCrit && ALLOW_PROD_CRIT(strat))
  {
    if (strat->fromT && (si.ecart>ecart))
    {
      if (g_paircrit_this_dispatch) {
        char *h_lm = paircrit_lm_str(p);
        char *s_lm = paircrit_lm_str(si.p);
        kt_paircrit_logf(
          "KILL tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s by=fromT_ecart\n",
          kt_debug_tid, (int)strat->S.size(),
          h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL");
        if (h_lm) omFree(h_lm);
        if (s_lm) omFree(s_lm);
      }
      return;
      /*the pair is (s[i],t[.]), discard it if the ecart is too big*/
    }
    if((!((si.ecart>0)&&(ecart>0)))
    && pHasNotCF(p,si.p))
    {
    /*
    *the product criterion has applied for (s,p),
    *i.e. lcm(s,p)=product of the leading terms of s and p.
    *Suppose (s,r) is in L and the leading term
    *of p divides lcm(s,r)
    *(==> the leading term of p divides the leading term of r)
    *but the leading term of s does not divide the leading term of r
    *(notice that tis condition is automatically satisfied if r is still
    *in S), then (s,r) can be cancelled.
    *This should be done here because the
    *case lcm(s,r)=lcm(s,p) is not covered by chainCrit.
    *
    *Moreover, skipping (s,r) holds also for the noncommutative case.
    */
      strat->cp++;
      if (g_paircrit_this_dispatch) {
        char *h_lm = paircrit_lm_str(p);
        char *s_lm = paircrit_lm_str(si.p);
        kt_paircrit_logf(
          "KILL tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s by=prodCrit (sugar)\n",
          kt_debug_tid, (int)strat->S.size(),
          h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL");
        if (h_lm) omFree(h_lm);
        if (s_lm) omFree(s_lm);
      }
      if (g_event_log_enabled) {
        kevlog_emit(EVT_KILL, (uint16_t)kt_debug_tid,
                    (uint16_t)strat->S.size(), 0,
                    (uint32_t)si.arrival_id, (uint32_t)KR_PROD_CRIT, 0, 0,
                    (const void *)p, (const void *)si.p);
      }
      return;
    }
  }
  else /*sugarcrit*/
  {
    if (ALLOW_PROD_CRIT(strat))
    {
      if (strat->fromT && (si.ecart>ecart))
      {
        if (g_paircrit_this_dispatch) {
          char *h_lm = paircrit_lm_str(p);
          char *s_lm = paircrit_lm_str(si.p);
          kt_paircrit_logf(
            "KILL tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s by=fromT_ecart\n",
            kt_debug_tid, (int)strat->S.size(),
            h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL");
          if (h_lm) omFree(h_lm);
          if (s_lm) omFree(s_lm);
        }
        return;
        /*the pair is (s[i],t[.]), discard it if the ecart is too big*/
      }
      // if currRing->nc_type!=quasi (or skew)
      // TODO: enable productCrit for super commutative algebras...
      if(/*(strat->ak==0) && productCrit(p,si.p)*/
      pHasNotCF(p,si.p))
      {
      /*
      *the product criterion has applied for (s,p),
      *i.e. lcm(s,p)=product of the leading terms of s and p.
      *Suppose (s,r) is in L and the leading term
      *of p divides lcm(s,r)
      *(==> the leading term of p divides the leading term of r)
      *but the leading term of s does not divide the leading term of r
      *(notice that tis condition is automatically satisfied if r is still
      *in S), then (s,r) can be canceled.
      *This should be done here because the
      *case lcm(s,r)=lcm(s,p) is not covered by chainCrit.
      */
          strat->cp++;
          if (g_paircrit_this_dispatch) {
            char *h_lm = paircrit_lm_str(p);
            char *s_lm = paircrit_lm_str(si.p);
            kt_paircrit_logf(
              "KILL tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s by=prodCrit\n",
              kt_debug_tid, (int)strat->S.size(),
              h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL");
            if (h_lm) omFree(h_lm);
            if (s_lm) omFree(s_lm);
          }
          if (g_event_log_enabled) {
            kevlog_emit(EVT_KILL, (uint16_t)kt_debug_tid,
                        (uint16_t)strat->S.size(), 0,
                        (uint32_t)si.arrival_id, (uint32_t)KR_PROD_CRIT, 0, 0,
                        (const void *)p, (const void *)si.p);
          }
          return;
      }
    }
  }

  /*- Only initialize Lp for pairs that survive the product criterion -*/
  LObject  Lp;
  Lp.i_r = -1;

#ifdef KDEBUG
  Lp.ecart=0; Lp.length=0;
#endif

  /*- computes the lcm(s[i],p) -*/
  Lp.lcm = pInit();

#ifndef HAVE_RATGRING
  pLcm(p,si.p,Lp.lcm);
#elif defined(HAVE_RATGRING)
  if (rIsRatGRing(currRing))
    pLcmRat(p,si.p,Lp.lcm, currRing->real_var_start); // int rat_shift
  else
    pLcm(p,si.p,Lp.lcm);
#endif
  pSetm(Lp.lcm);
  // sev(lcm(a,b)) ⊇ sev(a) | sev(b): the OR is a safe overapproximation
  // (more bits set = fewer pre-filter rejections in chainCritNormal, but
  // no false rejections since sev is only used for pre-filtering).
  // This avoids walking the LCM exponent vector via p_GetShortExpVector.
  Lp.sev_lcm = p_GetShortExpVector(p, currRing) | si.sev;

  if (strat->sugarCrit && ALLOW_PROD_CRIT(strat))
  {
    Lp.ecart = si_max(ecart,si.ecart);
    /*
    *the set B collects the pairs of type (S[j],p)
    *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p)#lcm(r,p)
    *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
    *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
    */
    {
      const unsigned long sev_lp = Lp.sev_lcm;
      for (auto it = strat_B(strat).ufbegin_lcm(sev_lp, sev_lp); it != strat_B(strat).ufend_lcm(); )
      {
        compare=pDivComp(it->lcm,Lp.lcm);
        if ((compare==1)
        &&(sugarDivisibleBy(it->ecart,Lp.ecart)))
        {
          strat->c3++;
          if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
          {
            if (g_paircrit_this_dispatch) {
              char *h_lm = paircrit_lm_str(p);
              char *s_lm = paircrit_lm_str(si.p);
              char *lp_lm = paircrit_lm_str(Lp.lcm);
              char *it_lm = paircrit_lm_str(it->lcm);
              kt_paircrit_logf(
                "KILL tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s lp_lcm=%s "
                "by=domByB victim_lcm=%s (sugar)\n",
                kt_debug_tid, (int)strat->S.size(),
                h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL",
                lp_lm ? lp_lm : "NULL", it_lm ? it_lm : "NULL");
              if (h_lm) omFree(h_lm);
              if (s_lm) omFree(s_lm);
              if (lp_lm) omFree(lp_lm);
              if (it_lm) omFree(it_lm);
            }
            pLmFree(Lp.lcm);
            return;
          }
          break;
        }
        else
        if ((compare ==-1)
        && sugarDivisibleBy(Lp.ecart,it->ecart))
        {
          if (g_paircrit_this_dispatch) {
            char *lp_lm = paircrit_lm_str(Lp.lcm);
            char *it_lm = paircrit_lm_str(it->lcm);
            kt_paircrit_logf(
              "BKILL tid=%d atT=%d victim_lcm=%s killer_lcm=%s (sugar)\n",
              kt_debug_tid, (int)strat->S.size(),
              it_lm ? it_lm : "NULL", lp_lm ? lp_lm : "NULL");
            if (lp_lm) omFree(lp_lm);
            if (it_lm) omFree(it_lm);
          }
          it = strat_B(strat).erase(it);
          strat->c3++;
        }
        else
          ++it;
      }
    }
  }
  else /*sugarcrit*/
  {
    if (ALLOW_PROD_CRIT(strat))
    {
      /*
      *the set B collects the pairs of type (S[j],p)
      *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p)#lcm(r,p)
      *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
      *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
      */
      const unsigned long sev_lp = Lp.sev_lcm;
      for (auto it = strat_B(strat).ufbegin_lcm(sev_lp, sev_lp); it != strat_B(strat).ufend_lcm(); )
      {
        compare=pDivComp(it->lcm,Lp.lcm);
        if (compare==1)
        {
          strat->c3++;
          if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
          {
            if (g_paircrit_this_dispatch) {
              char *h_lm = paircrit_lm_str(p);
              char *s_lm = paircrit_lm_str(si.p);
              char *lp_lm = paircrit_lm_str(Lp.lcm);
              char *it_lm = paircrit_lm_str(it->lcm);
              kt_paircrit_logf(
                "KILL tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s lp_lcm=%s "
                "by=domByB victim_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL",
                lp_lm ? lp_lm : "NULL", it_lm ? it_lm : "NULL");
              if (h_lm) omFree(h_lm);
              if (s_lm) omFree(s_lm);
              if (lp_lm) omFree(lp_lm);
              if (it_lm) omFree(it_lm);
            }
            pLmFree(Lp.lcm);
            return;
          }
          break;
        }
        else
        if (compare ==-1)
        {
          if (g_paircrit_this_dispatch) {
            char *lp_lm = paircrit_lm_str(Lp.lcm);
            char *it_lm = paircrit_lm_str(it->lcm);
            kt_paircrit_logf(
              "BKILL tid=%d atT=%d victim_lcm=%s killer_lcm=%s\n",
              kt_debug_tid, (int)strat->S.size(),
              it_lm ? it_lm : "NULL", lp_lm ? lp_lm : "NULL");
            if (lp_lm) omFree(lp_lm);
            if (it_lm) omFree(it_lm);
          }
          it = strat_B(strat).erase(it);
          strat->c3++;
        }
        else
          ++it;
      }
    }
  }
  /*
  *the pair (S[i],p) enters B if the spoly != 0
  */
  /*-  compute the short s-polynomial -*/
  if (strat->fromT && !TEST_OPT_INTSTRATEGY)
    pNorm(p);

  if ((si.p==NULL) || (p==NULL))
    return;

  if ((strat->hasFromQ) && (isFromQ!=0) && (si.fromQ!=0))
    Lp.p=NULL;
  else
  {
    #ifdef HAVE_PLURAL
    if ( rIsPluralRing(currRing) )
    {
      if(pHasNotCF(p, si.p))
      {
        if(ncRingType(currRing) == nc_lie)
        {
          // generalized prod-crit for lie-type
          strat->cp++;
          Lp.p = nc_p_Bracket_qq(pCopy(p),si.p, currRing);
        }
        else
        if( ALLOW_PROD_CRIT(strat) )
        {
          // product criterion for homogeneous case in SCA
          strat->cp++;
          Lp.p = NULL;
        }
        else
        {
          Lp.p = // nc_CreateSpoly(strat->S[i].p,p,currRing);
               nc_CreateShortSpoly(si.p, p, currRing);
          assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
          pNext(Lp.p) = strat->tail; // !!!
        }
      }
      else
      {
        Lp.p = // nc_CreateSpoly(strat->S[i].p,p,currRing);
              nc_CreateShortSpoly(si.p, p, currRing);

        assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
        pNext(Lp.p) = strat->tail; // !!!
      }
    }
    else
    #endif
    {
      assume(!rIsPluralRing(currRing));
      Lp.p = ksCreateShortSpoly(si.p, p, strat->tailRing);
    }
  }
  if (Lp.p == NULL)
  {
    if (g_paircrit_this_dispatch) {
      char *h_lm = paircrit_lm_str(p);
      char *s_lm = paircrit_lm_str(si.p);
      char *lp_lm = paircrit_lm_str(Lp.lcm);
      kt_paircrit_logf(
        "ZERO tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s lcm=%s\n",
        kt_debug_tid, (int)strat->S.size(),
        h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL",
        lp_lm ? lp_lm : "NULL");
      if (h_lm) omFree(h_lm);
      if (s_lm) omFree(s_lm);
      if (lp_lm) omFree(lp_lm);
    }
    /*- the case that the s-poly is 0 -*/
    record_pairtest_hit(si, strat);
    /*hint for spoly(S[i],p) == 0 for some i,0 <= i <= sl*/
    /*
    *suppose we have (s,r),(r,p),(s,p) and spoly(s,p) == 0 and (r,p) is
    *still in B (i.e. lcm(r,p) == lcm(s,p) or the leading term of s does not
    *divide lcm(r,p)). In the last case (s,r) can be canceled if the leading
    *term of p divides the lcm(s,r)
    *(this canceling should be done here because
    *the case lcm(s,p) == lcm(s,r) is not covered in chainCrit)
    *the first case is handled in chainCrit
    */
    if (Lp.lcm!=NULL) pLmFree(Lp.lcm);
  }
  else
  {
    /*- the pair (S[i],p) enters B -*/
    Lp.p1 = si.p;
    Lp.p2 = p;

    if (
        (!rIsPluralRing(currRing))
//      ||  (rIsPluralRing(currRing) && (ncRingType(currRing) != nc_lie))
       )
    {
      assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
      pNext(Lp.p) = strat->tail; // !!!
    }

    if (atR >= 0)
    {
      Lp.i_r1 = si.s_2_r;
      Lp.i_r2 = atR;
      // Creation-time CONSISTENCY check.
      // atR is an R-SLOT INDEX; the authoritative check is R[atR]->p.
      // (Historically this compared against T[atR].p, but T-array
      // positions shift under posInT-driven sorted insertion while
      // R slots are stable.  See progress 20Apr2026-1030.md.)
      if (atR < strat->T.size())
      {
        TObject *r_entry = strat->R[atR];
        poly R_atR_p = (r_entry != NULL) ? r_entry->p : NULL;
        poly T_atR_p = strat->T[atR].p;
        if (R_atR_p != NULL && R_atR_p != p)
          kt_debug_tag(
            "enterOnePair:R[atR]->p!=p (BORN_R_INCONSISTENT, REAL)",
            (void*)p, atR, si.s_2_r);
        if (T_atR_p != p && R_atR_p == p)
          kt_debug_tag(
            "enterOnePair:T[atR]!=p (cosmetic; R consistent)",
            (void*)p, atR, si.s_2_r);
      }
      else
      {
        // atR out of bounds at creation time — pair is created with
        // i_r2 pointing beyond current T array.  This is a race with
        // T not yet grown to atR+1.
        kt_debug_tag("enterOnePair:atR>=T.size (SKIPPED_CHECK)",
                     (void*)p, atR, (int)strat->T.size());
      }
    }
    else
    {
      Lp.i_r1 = -1;
      Lp.i_r2 = -1;
    }
    strat->initEcartPair(&Lp,si.p,p,si.ecart,ecart);

    if (TEST_OPT_INTSTRATEGY)
    {
      if (!rIsPluralRing(currRing)
      && !rField_is_Ring(currRing)
      && (Lp.p->coef!=NULL))
        nDelete(&(Lp.p->coef));
    }

    // Fingerprint experiment (off-by-one pair-mutation probe).
    // Stamp Lp.dbg_fp from (p1, p2, i_r1, i_r2) just before the push
    // into strat->L.  If at pop / L-scan time the recomputed fp still
    // matches, the four fields were preserved as a set (bug must be
    // elsewhere — e.g., T-side shift).  If the fp differs, at least
    // one of the four fields was mutated in flight inside strat->L.
    Lp.dbg_fp = kt_debug_pair_fp((void*)Lp.p1, (void*)Lp.p2,
                                 Lp.i_r1, Lp.i_r2);

    if (g_paircrit_this_dispatch) {
      char *h_lm = paircrit_lm_str(p);
      char *s_lm = paircrit_lm_str(si.p);
      char *lp_lm = paircrit_lm_str(Lp.lcm);
      kt_paircrit_logf(
        "KEEP tid=%d atT=%d h_lm=%s si_s2r=%d si_lm=%s lcm=%s ecart=%d\n",
        kt_debug_tid, (int)strat->S.size(),
        h_lm ? h_lm : "NULL", (int)si.s_2_r, s_lm ? s_lm : "NULL",
        lp_lm ? lp_lm : "NULL", (int)Lp.ecart);
      if (h_lm) omFree(h_lm);
      if (s_lm) omFree(s_lm);
      if (lp_lm) omFree(lp_lm);
    }
    // Task 325 event log: KEEP.  Pair survived all pre-insert kills
    // and is about to enter B.  arg_a = si.arrival_id, arg_c = ecart.
    // poly_ptr_2 = si.p (second base), poly_ptr_1 = h.
    if (g_event_log_enabled) {
      kevlog_emit(EVT_KEEP, (uint16_t)kt_debug_tid,
                  (uint16_t)strat->S.size(), 0,
                  (uint32_t)si.arrival_id, 0, (uint32_t)Lp.ecart, 0,
                  (const void *)p, (const void *)si.p);
      kevlog_register_poly(Lp.lcm);
    }
    strat_B(strat).push(Lp);
  }
}

/// p_HasNotCF for the IDLIFT case and syzComp==1: ignore component
static inline BOOLEAN p_HasNotCF_Lift(poly p1, poly p2, const ring r)
{
  int i = rVar(r);
  loop
  {
    if ((p_GetExp(p1, i, r) > 0) && (p_GetExp(p2, i, r) > 0))
      return FALSE;
    i--;
    if (i == 0)
      return TRUE;
  }
}

/*2
* put the pair (s[i],p)  into the set B, ecart=ecart(p) for idLift(I,T)
*  (in the special case: idLift for ideals, i.e. strat->syzComp==1)
*  (prod.crit applies)
*/

static void enterOnePairLift (const SElement &si,poly p,int ecart, int isFromQ,kStrategy strat, int atR = -1)
{
  assume(ALLOW_PROD_CRIT(strat));
  assume(!rIsPluralRing(currRing));

  assume(strat->syzComp==1);

  if ((si.p==NULL) || (p==NULL))
    return;

  int      compare;
  LObject  Lp;
  Lp.i_r = -1;

#ifdef KDEBUG
  Lp.ecart=0; Lp.length=0;
#endif
  /*- computes the lcm(s[i],p) -*/
  Lp.lcm = p_Lcm(p,si.p,currRing);
  Lp.sev_lcm = p_GetShortExpVector(Lp.lcm, currRing);

  if (strat->sugarCrit)
  {
    if((!((si.ecart>0)&&(ecart>0)))
    && p_HasNotCF_Lift(p,si.p,currRing))
    {
    /*
    *the product criterion has applied for (s,p),
    *i.e. lcm(s,p)=product of the leading terms of s and p.
    *Suppose (s,r) is in L and the leading term
    *of p divides lcm(s,r)
    *(==> the leading term of p divides the leading term of r)
    *but the leading term of s does not divide the leading term of r
    *(notice that tis condition is automatically satisfied if r is still
    *in S), then (s,r) can be cancelled.
    *This should be done here because the
    *case lcm(s,r)=lcm(s,p) is not covered by chainCrit.
    *
    *Moreover, skipping (s,r) holds also for the noncommutative case.
    */
      strat->cp++;
      pLmFree(Lp.lcm);
      return;
    }
    else
      Lp.ecart = si_max(ecart,si.ecart);
    if (strat->fromT && (si.ecart>ecart))
    {
      pLmFree(Lp.lcm);
      return;
      /*the pair is (s[i],t[.]), discard it if the ecart is too big*/
    }
    /*
    *the set B collects the pairs of type (S[j],p)
    *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p)#lcm(r,p)
    *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
    *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
    */
    {
      for (auto jt = strat_B(strat).ubegin(); jt != strat_B(strat).uend(); )
      {
        compare=pDivComp(jt->lcm,Lp.lcm);
        if ((compare==1)
        &&(sugarDivisibleBy(jt->ecart,Lp.ecart)))
        {
          strat->c3++;
          if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
          {
            pLmFree(Lp.lcm);
            return;
          }
          break;
        }
        else
        if ((compare ==-1)
        && sugarDivisibleBy(Lp.ecart,jt->ecart))
        {
          jt = strat_B(strat).erase(jt);
          strat->c3++;
        }
        else
          ++jt;
      }
    }
  }
  else /*sugarcrit*/
  {
    if(/*(strat->ak==0) && productCrit(p,si.p)*/
    p_HasNotCF_Lift(p,si.p,currRing))
    {
    /*
    *the product criterion has applied for (s,p),
    *i.e. lcm(s,p)=product of the leading terms of s and p.
    *Suppose (s,r) is in L and the leading term
    *of p divides lcm(s,r)
    *(==> the leading term of p divides the leading term of r)
    *but the leading term of s does not divide the leading term of r
    *(notice that tis condition is automatically satisfied if r is still
    *in S), then (s,r) can be canceled.
    *This should be done here because the
    *case lcm(s,r)=lcm(s,p) is not covered by chainCrit.
    */
      strat->cp++;
      pLmFree(Lp.lcm);
      return;
    }
    if (strat->fromT && (si.ecart>ecart))
    {
      pLmFree(Lp.lcm);
      return;
      /*the pair is (s[i],t[.]), discard it if the ecart is too big*/
    }
    /*
    *the set B collects the pairs of type (S[j],p)
    *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p)#lcm(r,p)
    *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
    *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
    */
    for(auto jt = strat_B(strat).ubegin(); jt != strat_B(strat).uend(); )
    {
      compare=pDivComp(jt->lcm,Lp.lcm);
      if (compare==1)
      {
        strat->c3++;
        if ((!strat->hasFromQ) || (isFromQ==0) || (si.fromQ==0))
        {
          pLmFree(Lp.lcm);
          return;
        }
        break;
      }
      else
      if (compare ==-1)
      {
        jt = strat_B(strat).erase(jt);
        strat->c3++;
      }
      else
        ++jt;
    }
  }
  /*
  *the pair (S[i],p) enters B if the spoly != 0
  */
  /*-  compute the short s-polynomial -*/
  if (strat->fromT && !TEST_OPT_INTSTRATEGY)
    pNorm(p);

  if ((strat->hasFromQ) && (isFromQ!=0) && (si.fromQ!=0))
    Lp.p=NULL;
  else
  {
    assume(!rIsPluralRing(currRing));
    Lp.p = ksCreateShortSpoly(si.p, p, strat->tailRing);
  }
  if (Lp.p == NULL)
  {
    /*- the case that the s-poly is 0 -*/
    record_pairtest_hit(si, strat);
    /*hint for spoly(S[i],p) == 0 for some i,0 <= i <= sl*/
    /*
    *suppose we have (s,r),(r,p),(s,p) and spoly(s,p) == 0 and (r,p) is
    *still in B (i.e. lcm(r,p) == lcm(s,p) or the leading term of s does not
    *divide lcm(r,p)). In the last case (s,r) can be canceled if the leading
    *term of p divides the lcm(s,r)
    *(this canceling should be done here because
    *the case lcm(s,p) == lcm(s,r) is not covered in chainCrit)
    *the first case is handled in chainCrit
    */
    if (Lp.lcm!=NULL) pLmFree(Lp.lcm);
  }
  else
  {
    /*- the pair (S[i],p) enters B -*/
    Lp.p1 = si.p;
    Lp.p2 = p;

    pNext(Lp.p) = strat->tail; // !!!

    if (atR >= 0)
    {
      Lp.i_r1 = si.s_2_r;
      Lp.i_r2 = atR;
      // Creation-time CONSISTENCY check.
      // atR is an R-SLOT INDEX; the authoritative check is R[atR]->p.
      // (Historically this compared against T[atR].p, but T-array
      // positions shift under posInT-driven sorted insertion while
      // R slots are stable.  See progress 20Apr2026-1030.md.)
      if (atR < strat->T.size())
      {
        TObject *r_entry = strat->R[atR];
        poly R_atR_p = (r_entry != NULL) ? r_entry->p : NULL;
        poly T_atR_p = strat->T[atR].p;
        if (R_atR_p != NULL && R_atR_p != p)
          kt_debug_tag(
            "enterOnePair:R[atR]->p!=p (BORN_R_INCONSISTENT, REAL)",
            (void*)p, atR, si.s_2_r);
        if (T_atR_p != p && R_atR_p == p)
          kt_debug_tag(
            "enterOnePair:T[atR]!=p (cosmetic; R consistent)",
            (void*)p, atR, si.s_2_r);
      }
      else
      {
        // atR out of bounds at creation time — pair is created with
        // i_r2 pointing beyond current T array.  This is a race with
        // T not yet grown to atR+1.
        kt_debug_tag("enterOnePair:atR>=T.size (SKIPPED_CHECK)",
                     (void*)p, atR, (int)strat->T.size());
      }
    }
    else
    {
      Lp.i_r1 = -1;
      Lp.i_r2 = -1;
    }
    strat->initEcartPair(&Lp,si.p,p,si.ecart,ecart);

    if (TEST_OPT_INTSTRATEGY)
    {
      if (!rIsPluralRing(currRing)
      && !rField_is_Ring(currRing)
      && (Lp.p->coef!=NULL))
        nDelete(&(Lp.p->coef));
    }

    strat_B(strat).push(Lp);
  }
}

/*2
* put the pair (s[i],p)  into the set B, ecart=ecart(p)
* NOTE: here we need to add the signature-based criteria
*/

#ifdef DEBUGF5
static void enterOnePairSig (const SElement &si, sBasisSet::const_iterator si_it, poly p, poly pSig, int from, int ecart, int isFromQ, kStrategy strat, int atR = -1)
#else
static void enterOnePairSig (const SElement &si, sBasisSet::const_iterator si_it, poly p, poly pSig, int, int ecart, int isFromQ, kStrategy strat, int atR = -1)
#endif
{


  poly m1 = NULL,m2 = NULL; // we need the multipliers for the s-polynomial to compute
              // the corresponding signatures for criteria checks
  LObject  Lp;
  poly pSigMult = p_Copy(pSig,currRing);
  poly sSigMult = p_Copy(si.sig,currRing);
  unsigned long pSigMultNegSev,sSigMultNegSev;
  Lp.i_r = -1;

#ifdef KDEBUG
  Lp.ecart=0; Lp.length=0;
#endif
  /*- computes the lcm(s[i],p) -*/
  Lp.lcm = pInit();
  k_GetLeadTerms(p,si.p,currRing,m1,m2,currRing);
#ifndef HAVE_RATGRING
  pLcm(p,si.p,Lp.lcm);
#elif defined(HAVE_RATGRING)
  if (rIsRatGRing(currRing))
    pLcmRat(p,si.p,Lp.lcm, currRing->real_var_start); // int rat_shift
  else
    pLcm(p,si.p,Lp.lcm);
#endif
  pSetm(Lp.lcm);
  Lp.sev_lcm = p_GetShortExpVector(Lp.lcm, currRing);

  // set coeffs of multipliers m1 and m2
  pSetCoeff0(m1, nInit(1));
  pSetCoeff0(m2, nInit(1));
//#if 1
#ifdef DEBUGF5
  PrintS("P1  ");
  pWrite(pHead(p));
  PrintS("P2  ");
  pWrite(pHead(si.p));
  PrintS("M1  ");
  pWrite(m1);
  PrintS("M2  ");
  pWrite(m2);
#endif
  // get multiplied signatures for testing
  pSigMult = currRing->p_Procs->pp_Mult_mm(pSigMult,m1,currRing);
  pSigMultNegSev = ~p_GetShortExpVector(pSigMult,currRing);
  sSigMult = currRing->p_Procs->pp_Mult_mm(sSigMult,m2,currRing);
  sSigMultNegSev = ~p_GetShortExpVector(sSigMult,currRing);

//#if 1
#ifdef DEBUGF5
  PrintS("----------------\n");
  pWrite(pSigMult);
  pWrite(sSigMult);
  PrintS("----------------\n");
  Lp.checked  = strat->S.cbegin();
#endif
  int sigCmp = p_LmCmp(pSigMult,sSigMult,currRing);
//#if 1
#if DEBUGF5
  Print("IN PAIR GENERATION - COMPARING SIGS: %d\n",sigCmp);
  pWrite(pSigMult);
  pWrite(sSigMult);
#endif
  if(sigCmp==0)
  {
    // printf("!!!!   EQUAL SIGS   !!!!\n");
    // pSig = sSig, delete element due to Rewritten Criterion
    pDelete(&pSigMult);
    pDelete(&sSigMult);
    if (rField_is_Ring(currRing))
      pLmDelete(Lp.lcm);
    else
      pLmFree(Lp.lcm);
    pDelete (&m1);
    pDelete (&m2);
    return;
  }
  // testing by syzCrit = F5 Criterion
  // testing by rewCrit1 = Rewritten Criterion
  // NOTE: Arri's Rewritten Criterion is tested below, we need Lp.p for it!
  auto si_next = si_it; ++si_next;
  if  ( strat->syzCrit(pSigMult,pSigMultNegSev,strat) ||
        strat->syzCrit(sSigMult,sSigMultNegSev,strat)
        || strat->rewCrit1(sSigMult,sSigMultNegSev,Lp.lcm,strat,si_next)
      )
  {
    pDelete(&pSigMult);
    pDelete(&sSigMult);
    if (rField_is_Ring(currRing))
      pLmDelete(Lp.lcm);
    else
      pLmFree(Lp.lcm);
    pDelete (&m1);
    pDelete (&m2);
    return;
  }
  /*
  *the pair (S[i],p) enters B if the spoly != 0
  */
  /*-  compute the short s-polynomial -*/
  if (strat->fromT && !TEST_OPT_INTSTRATEGY)
    pNorm(p);

  if ((si.p==NULL) || (p==NULL))
    return;

  if ((strat->hasFromQ) && (isFromQ!=0) && (si.fromQ!=0))
    Lp.p=NULL;
  else
  {
    #ifdef HAVE_PLURAL
    if ( rIsPluralRing(currRing) )
    {
      if(pHasNotCF(p, si.p))
      {
        if(ncRingType(currRing) == nc_lie)
        {
          // generalized prod-crit for lie-type
          strat->cp++;
          Lp.p = nc_p_Bracket_qq(pCopy(p),si.p, currRing);
        }
        else
        if( ALLOW_PROD_CRIT(strat) )
        {
          // product criterion for homogeneous case in SCA
          strat->cp++;
          Lp.p = NULL;
        }
        else
        {
          Lp.p = // nc_CreateSpoly(strat->S[i].p,p,currRing);
                nc_CreateShortSpoly(si.p, p, currRing);

          assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
          pNext(Lp.p) = strat->tail; // !!!
        }
      }
      else
      {
        Lp.p = // nc_CreateSpoly(strat->S[i].p,p,currRing);
              nc_CreateShortSpoly(si.p, p, currRing);

        assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
        pNext(Lp.p) = strat->tail; // !!!
      }
    }
    else
    #endif
    {
      assume(!rIsPluralRing(currRing));
      Lp.p = ksCreateShortSpoly(si.p, p, strat->tailRing);
    }
  }
  // store from which element this pair comes from for further tests
  //Lp.from = strat->S.size();
  if(sigCmp==currRing->OrdSgn)
  {
    // pSig > sSig
    pDelete (&sSigMult);
    Lp.sig    = pSigMult;
    Lp.sevSig = ~pSigMultNegSev;
  }
  else
  {
    // pSig < sSig
    pDelete (&pSigMult);
    Lp.sig    = sSigMult;
    Lp.sevSig = ~sSigMultNegSev;
  }
  if (Lp.p == NULL)
  {
    if (Lp.lcm!=NULL) pLmFree(Lp.lcm);
    int pos = posInSyz(strat, Lp.sig);
    enterSyz(Lp, strat, pos);
  }
  else
  {
    // testing by rewCrit3 = Arris Rewritten Criterion (for F5 nothing happens!)
    if (strat->rewCrit3(Lp.sig,~Lp.sevSig,Lp.p,strat,strat->S.cend()))
    {
      pLmFree(Lp.lcm);
      pDelete(&Lp.sig);
      pDelete (&m1);
      pDelete (&m2);
      return;
    }
    // in any case Lp is checked up to the next strat->P which is added
    // to S right after this critical pair creation.
    // NOTE: this even holds if the 2nd generator gives the bigger signature
    //       moreover, this improves rewCriterion,
    //       i.e. strat->checked > strat->from if and only if the 2nd generator
    //       gives the bigger signature.
    // Capture an iterator at end() — the next enter_bba (which happens
    // immediately after this pair creation) will append at this position
    // and the iterator will correctly point at that newly-appended
    // element at resume time. See sBasisSet iterator-invalidation
    // contract (rule 1: append does not invalidate prior iterators).
    Lp.checked = strat->S.cend();
    // at this point it is clear that the pair will be added to L, since it has
    // passed all tests up to now

  // adds buchberger's first criterion
    if (pLmCmp(m2,pHead(p)) == 0)
    {
      Lp.prod_crit = TRUE; // Product Criterion
#if 0
      int pos = posInSyz(strat, Lp.sig);
      enterSyz(Lp, strat, pos);
      pDelete (&m1);
      pDelete (&m2);
      return;
#endif
    }
    pDelete (&m1);
    pDelete (&m2);
#if DEBUGF5
    PrintS("SIGNATURE OF PAIR:  ");
    pWrite(Lp.sig);
#endif
    /*- the pair (S[i],p) enters B -*/
    Lp.p1 = si.p;
    Lp.p2 = p;

    if (
        (!rIsPluralRing(currRing))
//      ||  (rIsPluralRing(currRing) && (ncRingType(currRing) != nc_lie))
       )
    {
      assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
      pNext(Lp.p) = strat->tail; // !!!
    }

    if (atR >= 0)
    {
      Lp.i_r1 = si.s_2_r;
      Lp.i_r2 = atR;
      // Creation-time CONSISTENCY check.
      // atR is an R-SLOT INDEX; the authoritative check is R[atR]->p.
      // (Historically this compared against T[atR].p, but T-array
      // positions shift under posInT-driven sorted insertion while
      // R slots are stable.  See progress 20Apr2026-1030.md.)
      if (atR < strat->T.size())
      {
        TObject *r_entry = strat->R[atR];
        poly R_atR_p = (r_entry != NULL) ? r_entry->p : NULL;
        poly T_atR_p = strat->T[atR].p;
        if (R_atR_p != NULL && R_atR_p != p)
          kt_debug_tag(
            "enterOnePair:R[atR]->p!=p (BORN_R_INCONSISTENT, REAL)",
            (void*)p, atR, si.s_2_r);
        if (T_atR_p != p && R_atR_p == p)
          kt_debug_tag(
            "enterOnePair:T[atR]!=p (cosmetic; R consistent)",
            (void*)p, atR, si.s_2_r);
      }
      else
      {
        // atR out of bounds at creation time — pair is created with
        // i_r2 pointing beyond current T array.  This is a race with
        // T not yet grown to atR+1.
        kt_debug_tag("enterOnePair:atR>=T.size (SKIPPED_CHECK)",
                     (void*)p, atR, (int)strat->T.size());
      }
    }
    else
    {
      Lp.i_r1 = -1;
      Lp.i_r2 = -1;
    }
    strat->initEcartPair(&Lp,si.p,p,si.ecart,ecart);

    if (TEST_OPT_INTSTRATEGY)
    {
      if (!rIsPluralRing(currRing)
      && !rField_is_Ring(currRing)
      && (Lp.p->coef!=NULL))
        nDelete(&(Lp.p->coef));
    }

    strat_B(strat).push(Lp);
  }
}


#ifdef DEBUGF5
static void enterOnePairSigRing (const SElement &si, poly p, poly pSig, int from, int ecart, int isFromQ, kStrategy strat, int atR = -1)
#else
static void enterOnePairSigRing (const SElement &si, poly p, poly pSig, int, int ecart, int isFromQ, kStrategy strat, int atR = -1)
#endif
{
  #if ALL_VS_JUST
  //Over rings, if we construct the strong pair, do not add the spair
  if(rField_is_Ring(currRing))
  {
    number s,t,d;
    d = n_ExtGcd(pGetCoeff(p), pGetCoeff(si.p), &s, &t, currRing->cf);

    if (!nIsZero(s) && !nIsZero(t))  // evtl. durch divBy tests ersetzen
    {
      nDelete(&d);
      nDelete(&s);
      nDelete(&t);
      return;
    }
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
  }
  #endif

  poly m1 = NULL,m2 = NULL; // we need the multipliers for the s-polynomial to compute
              // the corresponding signatures for criteria checks
  LObject  Lp;
  poly pSigMult = p_Copy(pSig,currRing);
  poly sSigMult = p_Copy(si.sig,currRing);
  unsigned long pSigMultNegSev,sSigMultNegSev;
  Lp.i_r = -1;

#ifdef KDEBUG
  Lp.ecart=0; Lp.length=0;
#endif
  /*- computes the lcm(s[i],p) -*/
  Lp.lcm = pInit();
  k_GetLeadTerms(p,si.p,currRing,m1,m2,currRing);
#ifndef HAVE_RATGRING
  pLcm(p,si.p,Lp.lcm);
#elif defined(HAVE_RATGRING)
  if (rIsRatGRing(currRing))
    pLcmRat(p,si.p,Lp.lcm, currRing->real_var_start); // int rat_shift
  else
    pLcm(p,si.p,Lp.lcm);
#endif
  pSetm(Lp.lcm);
  Lp.sev_lcm = p_GetShortExpVector(Lp.lcm, currRing);

  // set coeffs of multipliers m1 and m2
  if(rField_is_Ring(currRing))
  {
    number s = nCopy(pGetCoeff(si.p));
    number t = nCopy(pGetCoeff(p));
    pSetCoeff0(Lp.lcm, n_Lcm(s, t, currRing->cf));
    ksCheckCoeff(&s, &t, currRing->cf);
    pSetCoeff0(m1,s);
    pSetCoeff0(m2,t);
  }
  else
  {
    pSetCoeff0(m1, nInit(1));
    pSetCoeff0(m2, nInit(1));
  }
#ifdef DEBUGF5
  Print("P1  ");
  pWrite(pHead(p));
  Print("P2  ");
  pWrite(pHead(si.p));
  Print("M1  ");
  pWrite(m1);
  Print("M2  ");
  pWrite(m2);
#endif

  // get multiplied signatures for testing
  pSigMult = pp_Mult_mm(pSigMult,m1,currRing);
  if(pSigMult != NULL)
    pSigMultNegSev = ~p_GetShortExpVector(pSigMult,currRing);
  sSigMult = pp_Mult_mm(sSigMult,m2,currRing);
  if(sSigMult != NULL)
    sSigMultNegSev = ~p_GetShortExpVector(sSigMult,currRing);
//#if 1
#ifdef DEBUGF5
  Print("----------------\n");
  pWrite(pSigMult);
  pWrite(sSigMult);
  Print("----------------\n");
  Lp.checked  = strat->S.cbegin();
#endif
  int sigCmp;
  if(pSigMult != NULL && sSigMult != NULL)
  {
    if(rField_is_Ring(currRing))
      sigCmp = p_LtCmpNoAbs(pSigMult,sSigMult,currRing);
    else
      sigCmp = p_LmCmp(pSigMult,sSigMult,currRing);
  }
  else
  {
    if(pSigMult == NULL)
    {
      if(sSigMult == NULL)
        sigCmp = 0;
      else
        sigCmp = -1;
    }
    else
      sigCmp = 1;
  }
//#if 1
#if DEBUGF5
  Print("IN PAIR GENERATION - COMPARING SIGS: %d\n",sigCmp);
  pWrite(pSigMult);
  pWrite(sSigMult);
#endif
  //In the ring case we already build the sig
  if(rField_is_Ring(currRing))
  {
    if(sigCmp == 0)
    {
      //sigdrop since we loose the signature
      strat->sigdrop = TRUE;
      //Try to reduce it as far as we can via redRing
      if(rField_is_Ring(currRing))
      {
        poly p1 = p_Copy(p,currRing);
        poly p2 = p_Copy(si.p,currRing);
        p1 = p_Mult_mm(p1,m1,currRing);
        p2 = p_Mult_mm(p2,m2,currRing);
        Lp.p = p_Sub(p1,p2,currRing);
        if(Lp.p != NULL)
          Lp.sev = p_GetShortExpVector(Lp.p,currRing);
      }
      int red_result = redRing(&Lp,strat);
      if(red_result == 0)
      {
        // Cancel the sigdrop
        p_Delete(&Lp.sig,currRing);Lp.sig = NULL;
        strat->sigdrop = FALSE;
        return;
      }
      else
      {
        strat->enterS(strat->P, strat, strat->T.size(), strat->S.end());
        #if 1
        strat->enterS(Lp, strat, strat->T.size()-1, strat->S.begin());
        #endif
        return;
      }
    }
    if(pSigMult != NULL && sSigMult != NULL && p_LmCmp(pSigMult,sSigMult,currRing) == 0)
    {
      //Same lm, have to subtract
      Lp.sig = p_Sub(pCopy(pSigMult),pCopy(sSigMult),currRing);
    }
    else
    {
      if(sigCmp == 1)
      {
        Lp.sig = pCopy(pSigMult);
      }
      if(sigCmp == -1)
      {
        Lp.sig = pNeg(pCopy(sSigMult));
      }
    }
    Lp.sevSig = p_GetShortExpVector(Lp.sig,currRing);
  }

  #if 0
  if(sigCmp==0)
  {
    // printf("!!!!   EQUAL SIGS   !!!!\n");
    // pSig = sSig, delete element due to Rewritten Criterion
    pDelete(&pSigMult);
    pDelete(&sSigMult);
    if (rField_is_Ring(currRing))
      pLmDelete(Lp.lcm);
    else
      pLmFree(Lp.lcm);
    pDelete (&m1);
    pDelete (&m2);
    return;
  }
  #endif
  // testing by syzCrit = F5 Criterion
  // testing by rewCrit1 = Rewritten Criterion
  // NOTE: Arri's Rewritten Criterion is tested below, we need Lp.p for it!
  if  ( strat->syzCrit(pSigMult,pSigMultNegSev,strat) ||
        strat->syzCrit(sSigMult,sSigMultNegSev,strat)
        // With this rewCrit activated i get a wrong deletion in sba_int_56.tst
        //|| strat->rewCrit1(sSigMult,sSigMultNegSev,Lp.lcm,strat,i+1)
      )
  {
    pDelete(&pSigMult);
    pDelete(&sSigMult);
    if (rField_is_Ring(currRing))
      pLmDelete(Lp.lcm);
    else
      pLmFree(Lp.lcm);
    pDelete (&m1);
    pDelete (&m2);
    return;
  }
  /*
  *the pair (S[i],p) enters B if the spoly != 0
  */
  /*-  compute the short s-polynomial -*/
  if (strat->fromT && !TEST_OPT_INTSTRATEGY)
    pNorm(p);

  if ((si.p==NULL) || (p==NULL))
    return;

  if ((strat->hasFromQ) && (isFromQ!=0) && (si.fromQ!=0))
    Lp.p=NULL;
  else
  {
    //Build p
    if(rField_is_Ring(currRing))
    {
      poly p1 = p_Copy(p,currRing);
      poly p2 = p_Copy(si.p,currRing);
      p1 = p_Mult_mm(p1,m1,currRing);
      p2 = p_Mult_mm(p2,m2,currRing);
      Lp.p = p_Sub(p1,p2,currRing);
      if(Lp.p != NULL)
        Lp.sev = p_GetShortExpVector(Lp.p,currRing);
    }
    else
    {
      #ifdef HAVE_PLURAL
      if ( rIsPluralRing(currRing) )
      {
        if(ncRingType(currRing) == nc_lie)
        {
          // generalized prod-crit for lie-type
          strat->cp++;
          Lp.p = nc_p_Bracket_qq(pCopy(p),si.p, currRing);
        }
        else
        if( ALLOW_PROD_CRIT(strat) )
        {
          // product criterion for homogeneous case in SCA
          strat->cp++;
          Lp.p = NULL;
        }
        else
        {
          Lp.p = // nc_CreateSpoly(strat->S[i].p,p,currRing);
                nc_CreateShortSpoly(si.p, p, currRing);

          assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
          pNext(Lp.p) = strat->tail; // !!!
        }
      }
      else
      #endif
      {
        assume(!rIsPluralRing(currRing));
        Lp.p = ksCreateShortSpoly(si.p, p, strat->tailRing);
      }
    }
  }
  // store from which element this pair comes from for further tests
  //Lp.from = strat->S.size();
  if(rField_is_Ring(currRing))
  {
    //Put the sig to be > 0
    if(!nGreaterZero(pGetCoeff(Lp.sig)))
    {
      Lp.sig = pNeg(Lp.sig);
      Lp.p = pNeg(Lp.p);
    }
  }
  else
  {
    if(sigCmp==currRing->OrdSgn)
    {
      // pSig > sSig
      pDelete (&sSigMult);
      Lp.sig    = pSigMult;
      Lp.sevSig = ~pSigMultNegSev;
    }
    else
    {
      // pSig < sSig
      pDelete (&pSigMult);
      Lp.sig    = sSigMult;
      Lp.sevSig = ~sSigMultNegSev;
    }
  }
  if (Lp.p == NULL)
  {
    if (Lp.lcm!=NULL) pLmFree(Lp.lcm);
    int pos = posInSyz(strat, Lp.sig);
    enterSyz(Lp, strat, pos);
  }
  else
  {
    // testing by rewCrit3 = Arris Rewritten Criterion (for F5 nothing happens!)
    if (strat->rewCrit3(Lp.sig,~Lp.sevSig,Lp.p,strat,strat->S.cend()))
    {
      pLmFree(Lp.lcm);
      pDelete(&Lp.sig);
      pDelete (&m1);
      pDelete (&m2);
      return;
    }
    // in any case Lp is checked up to the next strat->P which is added
    // to S right after this critical pair creation.
    // NOTE: this even holds if the 2nd generator gives the bigger signature
    //       moreover, this improves rewCriterion,
    //       i.e. strat->checked > strat->from if and only if the 2nd generator
    //       gives the bigger signature.
    // Capture an iterator at end() — the next enter_bba (which happens
    // immediately after this pair creation) will append at this position
    // and the iterator will correctly point at that newly-appended
    // element at resume time. See sBasisSet iterator-invalidation
    // contract (rule 1: append does not invalidate prior iterators).
    Lp.checked = strat->S.cend();
    // at this point it is clear that the pair will be added to L, since it has
    // passed all tests up to now

  // adds buchberger's first criterion
    if (pLmCmp(m2,pHead(p)) == 0)
    {
      Lp.prod_crit = TRUE; // Product Criterion
#if 0
      int pos = posInSyz(strat, Lp.sig);
      enterSyz(Lp, strat, pos);
      pDelete (&m1);
      pDelete (&m2);
      return;
#endif
    }
    pDelete (&m1);
    pDelete (&m2);
#if DEBUGF5
    PrintS("SIGNATURE OF PAIR:  ");
    pWrite(Lp.sig);
#endif
    /*- the pair (S[i],p) enters B -*/
    Lp.p1 = si.p;
    Lp.p2 = p;

    if (
        (!rIsPluralRing(currRing))
//      ||  (rIsPluralRing(currRing) && (ncRingType(currRing) != nc_lie))
      && !rField_is_Ring(currRing)
       )
    {
      assume(pNext(Lp.p)==NULL); // TODO: this may be violated whenever ext.prod.crit. for Lie alg. is used
      pNext(Lp.p) = strat->tail; // !!!
    }

    if (atR >= 0)
    {
      Lp.i_r1 = si.s_2_r;
      Lp.i_r2 = atR;
      // Creation-time CONSISTENCY check.
      // atR is an R-SLOT INDEX; the authoritative check is R[atR]->p.
      // (Historically this compared against T[atR].p, but T-array
      // positions shift under posInT-driven sorted insertion while
      // R slots are stable.  See progress 20Apr2026-1030.md.)
      if (atR < strat->T.size())
      {
        TObject *r_entry = strat->R[atR];
        poly R_atR_p = (r_entry != NULL) ? r_entry->p : NULL;
        poly T_atR_p = strat->T[atR].p;
        if (R_atR_p != NULL && R_atR_p != p)
          kt_debug_tag(
            "enterOnePair:R[atR]->p!=p (BORN_R_INCONSISTENT, REAL)",
            (void*)p, atR, si.s_2_r);
        if (T_atR_p != p && R_atR_p == p)
          kt_debug_tag(
            "enterOnePair:T[atR]!=p (cosmetic; R consistent)",
            (void*)p, atR, si.s_2_r);
      }
      else
      {
        // atR out of bounds at creation time — pair is created with
        // i_r2 pointing beyond current T array.  This is a race with
        // T not yet grown to atR+1.
        kt_debug_tag("enterOnePair:atR>=T.size (SKIPPED_CHECK)",
                     (void*)p, atR, (int)strat->T.size());
      }
    }
    else
    {
      Lp.i_r1 = -1;
      Lp.i_r2 = -1;
    }
    strat->initEcartPair(&Lp,si.p,p,si.ecart,ecart);

    if (TEST_OPT_INTSTRATEGY)
    {
      if (!rIsPluralRing(currRing)
      && !rField_is_Ring(currRing)
      && (Lp.p->coef!=NULL))
        nDelete(&(Lp.p->coef));
    }
    // Check for sigdrop
    if(rField_is_Ring(currRing) && pLtCmp(Lp.sig,pSig) == -1)
    {
      strat->sigdrop = TRUE;
      // Completely reduce it
      int red_result = redRing(&Lp,strat);
      if(red_result == 0)
      {
        // Reduced to 0
        strat->sigdrop = FALSE;
        p_Delete(&Lp.sig,currRing);Lp.sig = NULL;
        return;
      }
      else
      {
        strat->enterS(strat->P, strat, strat->T.size(), strat->S.end());
        // 0 - add just the original poly causing the sigdrop, 1 - add also this
        #if 1
        strat->enterS(Lp, strat, strat->T.size(), strat->S.begin());
        #endif
        return;
      }
    }
    strat->L.push(Lp);
  }
}

/*2
* put the pair (s[i],p) into the set L, ecart=ecart(p)
* in the case that s forms a SB of (s)
*/
void enterOnePairSpecial (const SElement &si,poly p,int ecart,kStrategy strat, int atR = -1)
{
  //PrintS("try ");wrp(strat->S[i].p);PrintS(" and ");wrp(p);PrintLn();
  if(pHasNotCF(p,si.p))
  {
    //PrintS("prod-crit\n");
    if(ALLOW_PROD_CRIT(strat))
    {
      //PrintS("prod-crit\n");
      strat->cp++;
      return;
    }
  }

  LObject  Lp;
  Lp.i_r = -1;

  Lp.lcm = p_Lcm(p,si.p,currRing);
  Lp.sev_lcm = p_GetShortExpVector(Lp.lcm, currRing);
  /*-  compute the short s-polynomial -*/

  #ifdef HAVE_PLURAL
  if (rIsPluralRing(currRing))
  {
    Lp.p = nc_CreateShortSpoly(si.p,p, currRing); // ??? strat->tailRing?
  }
  else
  #endif
    Lp.p = ksCreateShortSpoly(si.p,p,strat->tailRing);

  if (Lp.p == NULL)
  {
     //PrintS("short spoly==NULL\n");
     pLmFree(Lp.lcm);
  }
  else
  {
    /*- the pair (S[i],p) enters L -*/
    Lp.p1 = si.p;
    Lp.p2 = p;
    if (atR >= 0)
    {
      Lp.i_r1 = si.s_2_r;
      Lp.i_r2 = atR;
    }
    else
    {
      Lp.i_r1 = -1;
      Lp.i_r2 = -1;
    }
    assume(pNext(Lp.p) == NULL);
    pNext(Lp.p) = strat->tail;
    strat->initEcartPair(&Lp,si.p,p,si.ecart,ecart);
    if (TEST_OPT_INTSTRATEGY)
    {
      if (!rIsPluralRing(currRing)
      && !rField_is_Ring(currRing)
      && (Lp.p->coef!=NULL))
        nDelete(&(Lp.p->coef));
    }
    strat->L.push(Lp);
  }
}

/*2
* merge set B into L
*/
void kMergeBintoL(kStrategy strat)
{
  while (!strat_B(strat).empty()) {
    auto Lobj = strat_B(strat).top();
    strat_B(strat).pop();
    if (g_reduce_this_dispatch) {
      char *lcm_lm = kt_lm_str(Lobj.lcm);
      char *p_lm = kt_lm_str(Lobj.p);
      kt_reduce_logf(
        "LINSERT tid=%d atT=%d site=kMergeBintoL i_r1=%d i_r2=%d "
        "lcm_lm=%s p_lm=%s ecart=%d FDeg=%ld length=%d\n",
        kt_debug_tid, (int)strat->T.size(),
        Lobj.i_r1, Lobj.i_r2,
        lcm_lm ? lcm_lm : "NULL", p_lm ? p_lm : "NULL",
        (int)Lobj.ecart, (long)Lobj.GetpFDeg(), (int)Lobj.length);
      if (lcm_lm) omFree(lcm_lm);
      if (p_lm) omFree(p_lm);
    }
    // Task 325 event log: LINSERT.  pointers are the two bases; lcm
    // ends up as a registered poly via kevlog_register_poly.
    if (g_event_log_enabled) {
      kevlog_emit(EVT_LINSERT, (uint16_t)kt_debug_tid,
                  (uint16_t)strat->T.size(), 0,
                  (uint32_t)Lobj.i_r1, (uint32_t)Lobj.i_r2,
                  (uint32_t)Lobj.ecart, 0,
                  (const void *)Lobj.p1, (const void *)Lobj.p2);
      kevlog_register_poly(Lobj.lcm);
    }
    strat->L.push(Lobj);
  }
  strat_B(strat).clear();  // reset flat_ array to prevent unbounded growth
}

/* merge set B into L, and return a vector of iterators pointing to the new
 * elements in L, guaranteed to be in the same order they appear in L
 *
 * The ordering is done to mimic previous versions of Singular so as
 * to ensure that regression tests pass.  I know of no other reason to
 * sort these iterators.
 */

std::vector<LSet::iterator> kMergeBintoL_and_return_iterators(kStrategy strat)
{
  std::vector<LSet::iterator> iterators;
  iterators.reserve(strat_B(strat).size());
  while (!strat_B(strat).empty()) {
    auto Lobj = strat_B(strat).top();
    strat_B(strat).pop();
    if (g_reduce_this_dispatch) {
      char *lcm_lm = kt_lm_str(Lobj.lcm);
      char *p_lm = kt_lm_str(Lobj.p);
      kt_reduce_logf(
        "LINSERT tid=%d atT=%d site=kMergeBintoL_ret i_r1=%d i_r2=%d "
        "lcm_lm=%s p_lm=%s ecart=%d FDeg=%ld length=%d\n",
        kt_debug_tid, (int)strat->T.size(),
        Lobj.i_r1, Lobj.i_r2,
        lcm_lm ? lcm_lm : "NULL", p_lm ? p_lm : "NULL",
        (int)Lobj.ecart, (long)Lobj.GetpFDeg(), (int)Lobj.length);
      if (lcm_lm) omFree(lcm_lm);
      if (p_lm) omFree(p_lm);
    }
    if (g_event_log_enabled) {
      kevlog_emit(EVT_LINSERT, (uint16_t)kt_debug_tid,
                  (uint16_t)strat->T.size(), 1, // flags=1 => _ret variant
                  (uint32_t)Lobj.i_r1, (uint32_t)Lobj.i_r2,
                  (uint32_t)Lobj.ecart, 0,
                  (const void *)Lobj.p1, (const void *)Lobj.p2);
      kevlog_register_poly(Lobj.lcm);
    }
    iterators.push_back(strat->L.push(Lobj));
  }
  // Sort iterators to match the ordering of their objects in L
  std::sort(iterators.begin(), iterators.end(),
    [&strat](LSet::iterator a, LSet::iterator b) {
      return strat->L.key_comp()(*a, *b);
    });
  strat_B(strat).clear();  // reset flat_ array to prevent unbounded growth
  return iterators;
}

/*2
*the pairset B of pairs of type (s[i],p) is complete now. It will be updated
*using the chain-criterion in B and L and enters B to L
*/
void chainCritNormal (poly p,int ecart,kStrategy strat)
{
  int j;
  unsigned long sev_p = p_GetShortExpVector(p, currRing);

  if (g_paircrit_this_dispatch) {
    char *p_lm = paircrit_lm_str(p);
    kt_paircrit_logf(
      "CHAIN_ENTER tid=%d atT=%d p_lm=%s ecart=%d B.size=%d S.size=%d L.size=%d "
      "Gebauer=%d fromT=%d sugarCrit=%d local_hits=%d\n",
      kt_debug_tid, (int)strat->S.size(),
      p_lm ? p_lm : "NULL", ecart,
      (int)strat_B(strat).size(), (int)strat->S.size(),
      (int)strat->L.size(),
      strat->Gebauer ? 1 : 0, strat->fromT ? 1 : 0,
      strat->sugarCrit ? 1 : 0,
      (t_local_pairtest_hits != NULL)
        ? (int)t_local_pairtest_hits->size() : -1);
    if (p_lm) omFree(p_lm);
  }

  // sev_flat_ is now maintained incrementally by LSet's insert/erase/pop/
  // reorder/copy/move methods.  No rebuild needed here.

  /*
  *pairtest[i] is TRUE if spoly(S[i],p) == 0.
  *In this case all elements in B such
  *that their lcm is divisible by the leading term of S[i] can be canceled
  *
  * Parallel phase-1 mode: iterate the thread-local pairtest_hits
  * vector (populated by this drainer's enterOnePair calls) so we don't
  * cross-contaminate with peer drainers' pairtest hits.  Serial mode
  * (t_local_pairtest_hits == NULL): use the historical S-scan + clear.
  */
  if (t_local_pairtest_hits != NULL)
  {
#ifdef HAVE_SHIFTBBA
    if (rIsLPRing(currRing))
    {
      for (SElement *sit : *t_local_pairtest_hits)
      {
        for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
        {
          if (pLPDivisibleBy(sit->p, it->lcm))
          {
            if (g_paircrit_this_dispatch) {
              char *k_lm = paircrit_lm_str(sit->p);
              char *v_lm = paircrit_lm_str(it->lcm);
              kt_paircrit_logf(
                "CHAINKILL tid=%d atT=%d src=local_hits_LP killer_lm=%s "
                "killer_arrival=%lu victim_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                k_lm ? k_lm : "NULL", (unsigned long)sit->arrival_id,
                v_lm ? v_lm : "NULL");
              if (k_lm) omFree(k_lm);
              if (v_lm) omFree(v_lm);
            }
            it = strat_B(strat).erase(it);
            strat->c3++;
          }
          else ++it;
        }
      }
    }
    else
#endif
    {
      for (SElement *sit : *t_local_pairtest_hits)
      {
        for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
        {
          if (!(sit->sev & ~it->sev_lcm)
              && pDivisibleBy(sit->p, it->lcm))
          {
            if (g_paircrit_this_dispatch) {
              char *k_lm = paircrit_lm_str(sit->p);
              char *v_lm = paircrit_lm_str(it->lcm);
              kt_paircrit_logf(
                "CHAINKILL tid=%d atT=%d src=local_hits killer_lm=%s "
                "killer_arrival=%lu victim_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                k_lm ? k_lm : "NULL", (unsigned long)sit->arrival_id,
                v_lm ? v_lm : "NULL");
              if (k_lm) omFree(k_lm);
              if (v_lm) omFree(v_lm);
            }
            it = strat_B(strat).erase(it);
            strat->c3++;
          }
          else ++it;
        }
      }
    }
    t_local_pairtest_hits->clear();
  }
  else if (strat->S.has_pairtest())
  {
#ifdef HAVE_SHIFTBBA
    // only difference is pLPDivisibleBy instead of pDivisibleBy
    if (rIsLPRing(currRing))
    {
      for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
      {
        if (selement_pairtest_load(*sit))
        {
          for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
          {
            if (pLPDivisibleBy(sit->p,it->lcm))
            {
              if (g_paircrit_this_dispatch) {
                char *k_lm = paircrit_lm_str(sit->p);
                char *v_lm = paircrit_lm_str(it->lcm);
                kt_paircrit_logf(
                  "CHAINKILL tid=%d atT=%d src=S_pairtest_LP killer_idx=%ld "
                  "killer_lm=%s victim_lcm=%s\n",
                  kt_debug_tid, (int)strat->S.size(),
                  (long)sit.index(),
                  k_lm ? k_lm : "NULL",
                  v_lm ? v_lm : "NULL");
                if (k_lm) omFree(k_lm);
                if (v_lm) omFree(v_lm);
              }
              it = strat_B(strat).erase(it);
              strat->c3++;
            }
            else
              ++it;
          }
        }
      }
    }
    else
#endif
    {
      /*- i.e. there is an i with pairtest[i]==TRUE -*/
      for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
      {
        if (selement_pairtest_load(*sit))
        {
          for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
          {
            if (!(sit->sev & ~it->sev_lcm)
            && pDivisibleBy(sit->p,it->lcm))
            {
              if (g_paircrit_this_dispatch) {
                char *k_lm = paircrit_lm_str(sit->p);
                char *v_lm = paircrit_lm_str(it->lcm);
                kt_paircrit_logf(
                  "CHAINKILL tid=%d atT=%d src=S_pairtest killer_idx=%ld "
                  "killer_lm=%s victim_lcm=%s\n",
                  kt_debug_tid, (int)strat->S.size(),
                  (long)sit.index(),
                  k_lm ? k_lm : "NULL",
                  v_lm ? v_lm : "NULL");
                if (k_lm) omFree(k_lm);
                if (v_lm) omFree(v_lm);
              }
              it = strat_B(strat).erase(it);
              strat->c3++;
            }
            else
              ++it;
          }
        }
      }
    }
    strat->S.clear_pairtest();
  }
  if (strat->Gebauer || strat->fromT)
  {
    if (strat->sugarCrit)
    {
    /*
    *suppose L[j] == (s,r) and p/lcm(s,r)
    *and lcm(s,r)#lcm(s,p) and lcm(s,r)#lcm(r,p)
    *and in case the sugar is o.k. then L[j] can be canceled
    */
      {
        for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
        {
          if (sugarDivisibleBy(ecart,it->ecart)
          && ((it->p == strat->tail) || (rHasGlobalOrdering(currRing)))
          && pCompareChain(p,it->p1,it->p2,it->lcm)
          && (it->p == strat->tail))
          {
            if (g_paircrit_this_dispatch) {
              char *p_lm = paircrit_lm_str(p);
              char *v_lm = paircrit_lm_str(it->lcm);
              kt_paircrit_logf(
                "LKILL tid=%d atT=%d src=pCompareChain_sugar_GM killer_lm=%s "
                "victim_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                p_lm ? p_lm : "NULL", v_lm ? v_lm : "NULL");
              if (p_lm) omFree(p_lm);
              if (v_lm) omFree(v_lm);
            }
            it = strat->L.erase(it);
            strat->c3++;
          }
          else
            ++it;
        }
      }
      /*
      *this is GEBAUER-MOELLER:
      *in B all elements with the same lcm except the "best"
      *(i.e. the last one in B with this property) will be canceled
      */
      {
        const unsigned long* sev = strat_B(strat).sev_flat_data();
        const size_t n = strat_B(strat).sev_flat_size();
        for (size_t i = 0; i < n; i++)
        {
          if (sev[i] == 0) continue;
          for (size_t j = i + 1; j < n; j++)
          {
            if (sev[j] == 0) continue;
            LObject* a = strat_B(strat).flat_ptr(i);
            LObject* b = strat_B(strat).flat_ptr(j);
            if (a == NULL || b == NULL) continue;
            if (pLmEqual(a->lcm,b->lcm))
            {
              strat->c3++;
              if (a->ecart < b->ecart)
              {
                if (g_paircrit_this_dispatch) {
                  char *v_lm = paircrit_lm_str(b->lcm);
                  char *k_lm = paircrit_lm_str(a->lcm);
                  kt_paircrit_logf(
                    "GMKILL tid=%d atT=%d src=sugar_GM_ecart victim_lcm=%s "
                    "killer_lcm=%s victim_ecart=%d killer_ecart=%d\n",
                    kt_debug_tid, (int)strat->S.size(),
                    v_lm ? v_lm : "NULL", k_lm ? k_lm : "NULL",
                    b->ecart, a->ecart);
                  if (v_lm) omFree(v_lm);
                  if (k_lm) omFree(k_lm);
                }
                strat_B(strat).erase(strat_B(strat).uiter_at(j));
              }
              else if (a->ecart > b->ecart)
              {
                if (g_paircrit_this_dispatch) {
                  char *v_lm = paircrit_lm_str(a->lcm);
                  char *k_lm = paircrit_lm_str(b->lcm);
                  kt_paircrit_logf(
                    "GMKILL tid=%d atT=%d src=sugar_GM_ecart victim_lcm=%s "
                    "killer_lcm=%s victim_ecart=%d killer_ecart=%d\n",
                    kt_debug_tid, (int)strat->S.size(),
                    v_lm ? v_lm : "NULL", k_lm ? k_lm : "NULL",
                    a->ecart, b->ecart);
                  if (v_lm) omFree(v_lm);
                  if (k_lm) omFree(k_lm);
                }
                strat_B(strat).erase(strat_B(strat).uiter_at(i));
                break;  // i is gone, move to next i
              }
              else
              {
                // Equal ecart: use key_comp() tiebreaker for determinism
                if (strat_B(strat).key_comp()(*a, *b))
                {
                  if (g_paircrit_this_dispatch) {
                    char *v_lm = paircrit_lm_str(b->lcm);
                    char *k_lm = paircrit_lm_str(a->lcm);
                    kt_paircrit_logf(
                      "GMKILL tid=%d atT=%d src=sugar_GM_tiebreak victim_lcm=%s "
                      "killer_lcm=%s\n",
                      kt_debug_tid, (int)strat->S.size(),
                      v_lm ? v_lm : "NULL", k_lm ? k_lm : "NULL");
                    if (v_lm) omFree(v_lm);
                    if (k_lm) omFree(k_lm);
                  }
                  strat_B(strat).erase(strat_B(strat).uiter_at(j));
                }
                else
                {
                  if (g_paircrit_this_dispatch) {
                    char *v_lm = paircrit_lm_str(a->lcm);
                    char *k_lm = paircrit_lm_str(b->lcm);
                    kt_paircrit_logf(
                      "GMKILL tid=%d atT=%d src=sugar_GM_tiebreak victim_lcm=%s "
                      "killer_lcm=%s\n",
                      kt_debug_tid, (int)strat->S.size(),
                      v_lm ? v_lm : "NULL", k_lm ? k_lm : "NULL");
                    if (v_lm) omFree(v_lm);
                    if (k_lm) omFree(k_lm);
                  }
                  strat_B(strat).erase(strat_B(strat).uiter_at(i));
                  break;
                }
              }
            }
          }
        }
      }
    }
    else /*sugarCrit*/
    {
      /*
      *suppose L[j] == (s,r) and p/lcm(s,r)
      *and lcm(s,r)#lcm(s,p) and lcm(s,r)#lcm(r,p)
      *and in case the sugar is o.k. then L[j] can be canceled
      */
      {
        for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
        {
          if (pCompareChain(p,it->p1,it->p2,it->lcm)
          && ((pNext(it->p) == strat->tail)||(rHasGlobalOrdering(currRing))))
          {
            if (g_paircrit_this_dispatch) {
              char *p_lm = paircrit_lm_str(p);
              char *v_lm = paircrit_lm_str(it->lcm);
              kt_paircrit_logf(
                "LKILL tid=%d atT=%d src=pCompareChain_GM killer_lm=%s "
                "victim_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                p_lm ? p_lm : "NULL", v_lm ? v_lm : "NULL");
              if (p_lm) omFree(p_lm);
              if (v_lm) omFree(v_lm);
            }
            it = strat->L.erase(it);
            strat->c3++;
          }
          else
            ++it;
        }
      }
      /*
      *this is GEBAUER-MOELLER:
      *in B all elements with the same lcm except the "best"
      *(i.e. the last one in B with this property) will be canceled
      */
      {
        const unsigned long* sev = strat_B(strat).sev_flat_data();
        const size_t n = strat_B(strat).sev_flat_size();
        for (size_t i = 0; i < n; i++)
        {
          if (sev[i] == 0) continue;
          for (size_t j = i + 1; j < n; j++)
          {
            if (sev[j] == 0) continue;
            LObject* a = strat_B(strat).flat_ptr(i);
            LObject* b = strat_B(strat).flat_ptr(j);
            if (a == NULL || b == NULL) continue;
            if (pLmEqual(a->lcm,b->lcm))
            {
              strat->c3++;
              // Erase the worse element; keep the one that sorts first
              if (strat_B(strat).key_comp()(*a, *b))
              {
                if (g_paircrit_this_dispatch) {
                  char *v_lm = paircrit_lm_str(b->lcm);
                  char *k_lm = paircrit_lm_str(a->lcm);
                  kt_paircrit_logf(
                    "GMKILL tid=%d atT=%d src=GM_tiebreak victim_lcm=%s "
                    "killer_lcm=%s\n",
                    kt_debug_tid, (int)strat->S.size(),
                    v_lm ? v_lm : "NULL", k_lm ? k_lm : "NULL");
                  if (v_lm) omFree(v_lm);
                  if (k_lm) omFree(k_lm);
                }
                strat_B(strat).erase(strat_B(strat).uiter_at(j));
              }
              else
              {
                if (g_paircrit_this_dispatch) {
                  char *v_lm = paircrit_lm_str(a->lcm);
                  char *k_lm = paircrit_lm_str(b->lcm);
                  kt_paircrit_logf(
                    "GMKILL tid=%d atT=%d src=GM_tiebreak victim_lcm=%s "
                    "killer_lcm=%s\n",
                    kt_debug_tid, (int)strat->S.size(),
                    v_lm ? v_lm : "NULL", k_lm ? k_lm : "NULL");
                  if (v_lm) omFree(v_lm);
                  if (k_lm) omFree(k_lm);
                }
                strat_B(strat).erase(strat_B(strat).uiter_at(i));
                break;  // i is gone
              }
            }
          }
        }
      }
    }
    /*
    *the elements of B enter L
    */
    kMergeBintoL(strat);
  }
  else
  {
    {
      for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
      {
        #ifdef HAVE_SHIFTBBA
        if ((it->p1!=NULL) &&
        pCompareChain(p,it->p1,it->p2,it->lcm))
        #else
        if (pCompareChain(p,it->p1,it->p2,it->lcm))
        #endif
        {
          if ((pNext(it->p) == strat->tail)||(rHasGlobalOrdering(currRing)))
          {
            if (g_paircrit_this_dispatch) {
              char *p_lm = paircrit_lm_str(p);
              char *v_lm = paircrit_lm_str(it->lcm);
              char *vp1_lm = paircrit_lm_str(it->p1);
              char *vp2_lm = paircrit_lm_str(it->p2);
              kt_paircrit_logf(
                "LKILL tid=%d atT=%d src=pCompareChain_nonGebauer killer_lm=%s "
                "victim_lcm=%s victim_p1_lm=%s victim_p2_lm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                p_lm ? p_lm : "NULL", v_lm ? v_lm : "NULL",
                vp1_lm ? vp1_lm : "NULL", vp2_lm ? vp2_lm : "NULL");
              if (p_lm) omFree(p_lm);
              if (v_lm) omFree(v_lm);
              if (vp1_lm) omFree(vp1_lm);
              if (vp2_lm) omFree(vp2_lm);
            }
            it = strat->L.erase(it);
            strat->c3++;
            continue;
          }
        }
        ++it;
      }
    }
    /*
    *this is our MODIFICATION of GEBAUER-MOELLER:
    *Merge B into L, then deduplicate the B-origin elements using
    *pair_index for O(1) triangle checks.  kMergeBintoL_and_return_iterators
    *gives us iterators to the new elements in L-order, avoiding a full
    *scan of L.
    */
    std::vector<LSet::iterator> bvec = kMergeBintoL_and_return_iterators(strat);
    /* Deduplicate: for each pair of B-origin elements with the same lcm,
     * do the triangle check via isInPairsetL (O(1) pair_index lookup).
     * Elements removed from consideration are nulled out (set to end())
     * in bvec rather than erased, to avoid memmove. */
    LSet::iterator endL = strat->L.end();
    for (size_t ji = 0; ji < bvec.size(); ji++)
    {
      if (bvec[ji] == endL) continue;
      for (size_t ii = ji + 1; ii < bvec.size(); ii++)
      {
        if (bvec[ii] != endL && pLmEqual(bvec[ji]->lcm, bvec[ii]->lcm))
        {
          /*B[i] could be canceled but we search for a better one to cancel*/
          strat->c3++;
          auto lt = bvec[ii] + 1;
          if (isInPairsetL(lt,bvec[ji]->p1,bvec[ii]->p1,strat)
          && (pNext(lt->p) == strat->tail)
          && (!pLmEqual(bvec[ii]->p,lt->p))
          && pDivisibleBy(p,lt->lcm))
          {
            if (g_paircrit_this_dispatch) {
              char *p_lm = paircrit_lm_str(p);
              char *v_lm = paircrit_lm_str(lt->lcm);
              char *bji_lm = paircrit_lm_str(bvec[ji]->lcm);
              kt_paircrit_logf(
                "LKILL tid=%d atT=%d src=bvec_triangle_L_erase killer_lm=%s "
                "victim_lcm=%s bvec[ji]_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                p_lm ? p_lm : "NULL", v_lm ? v_lm : "NULL",
                bji_lm ? bji_lm : "NULL");
              if (p_lm) omFree(p_lm);
              if (v_lm) omFree(v_lm);
              if (bji_lm) omFree(bji_lm);
            }
            /*
            *"NOT equal(...)" because in case of "equal" the element L[l]
            *is "older" and has to be from theoretical point of view behind
            *L[i], but we do not want to reorder L
            */
            strat->L.erase(lt);
            /*
            *L[l] will be canceled, we cannot cancel L[i] later on,
            *so we null it out in bvec to prevent re-processing
            */
            bvec[ii] = endL;
          }
          else
          {
            if (g_paircrit_this_dispatch) {
              char *p_lm = paircrit_lm_str(p);
              char *v_lm = paircrit_lm_str(bvec[ii]->lcm);
              char *bji_lm = paircrit_lm_str(bvec[ji]->lcm);
              kt_paircrit_logf(
                "LKILL tid=%d atT=%d src=bvec_triangle_bvec_erase killer_lm=%s "
                "victim_lcm=%s bvec[ji]_lcm=%s\n",
                kt_debug_tid, (int)strat->S.size(),
                p_lm ? p_lm : "NULL", v_lm ? v_lm : "NULL",
                bji_lm ? bji_lm : "NULL");
              if (p_lm) omFree(p_lm);
              if (v_lm) omFree(v_lm);
              if (bji_lm) omFree(bji_lm);
            }
            strat->L.erase(bvec[ii]);
            bvec[ii] = endL;
          }
        }
      }
    }
  }
}
/*2
*the pairset B of pairs of type (s[i],p) is complete now. It will be updated
*without the chain-criterion in B and L and enters B to L
*/
void chainCritOpt_1 (poly,int,kStrategy strat)
{
  if (strat->S.has_pairtest())
    strat->S.clear_pairtest();
  /*
  *the elements of B enter L
  */
  kMergeBintoL(strat);
}
/*2
*the pairset B of pairs of type (s[i],p) is complete now. It will be updated
*using the chain-criterion in B and L and enters B to L
*/
void chainCritSig (poly p,int /*ecart*/,kStrategy strat)
{
  /*
  *Merge B into L, then deduplicate B-origin elements with the same LCM.
  *Uses kMergeBintoL_and_return_iterators + pair_index for O(1) triangle
  *checks, iterating only over B-origin elements via bvec.
  *
  *Only B-origin elements have p2==p; no tail-marking needed since bvec
  *nulling (bvec[ii]=endL) prevents re-processing.
  */
  std::vector<LSet::iterator> bvec = kMergeBintoL_and_return_iterators(strat);
  LSet::iterator endL = strat->L.end();
  for (size_t ji = 0; ji < bvec.size(); ji++)
  {
    if (bvec[ji] == endL) continue;
    for (size_t ii = ji + 1; ii < bvec.size(); ii++)
    {
      if (bvec[ii] == endL) continue;
      /* sev equality pre-filter: if sev_lcm values differ, LCMs can't
       * be equal, so skip the expensive pLmEqual call */
      if ((bvec[ji]->sev_lcm == bvec[ii]->sev_lcm)
          && pLmEqual(bvec[ji]->lcm, bvec[ii]->lcm))
      {
        /*bvec[ii] could be canceled but we search for a better one to cancel*/
        strat->c3++;
        auto lt = bvec[ii] + 1;
        if (isInPairsetL(lt,bvec[ji]->p1,bvec[ii]->p1,strat)
        && (pNext(lt->p) == strat->tail)
        && (!pLmEqual(bvec[ii]->p,lt->p))
        && pDivisibleBy(p,lt->lcm))
        {
          /*
          *"NOT equal(...)" because in case of "equal" the element L[l]
          *is "older" and has to be from theoretical point of view behind
          *L[i], but we do not want to reorder L
          */
          strat->L.erase(lt);
          /*
          *L[l] will be canceled, we cannot cancel L[i] later on,
          *so we null it out in bvec to prevent re-processing
          */
          bvec[ii] = endL;
        }
        else
        {
          strat->L.erase(bvec[ii]);
          bvec[ii] = endL;
        }
      }
    }
  }
}
#ifdef HAVE_RATGRING
void chainCritPart (poly p,int ecart,kStrategy strat)
{
  int j;
  unsigned long sev_p = p_GetShortExpVector(p, currRing);

  /*
  *pairtest[i] is TRUE if spoly(S[i],p) == 0.
  *In this case all elements in B such
  *that their lcm is divisible by the leading term of S[i] can be canceled
  *
  * Parallel phase-1: iterate thread-local hits (if set).
  */
  if (t_local_pairtest_hits != NULL)
  {
    for (SElement *sit : *t_local_pairtest_hits)
    {
      for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
      {
        if (_p_LmDivisibleByPart(sit->p,currRing,
           it->lcm,currRing,
           currRing->real_var_start,currRing->real_var_end))
        {
          it = strat_B(strat).erase(it);
          strat->c3++;
        }
        else ++it;
      }
    }
    t_local_pairtest_hits->clear();
  }
  else if (strat->S.has_pairtest())
  {
    /*- i.e. there is an i with pairtest[i]==TRUE -*/
    for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
    {
      if (selement_pairtest_load(*sit))
      {
        for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
        {
          if (_p_LmDivisibleByPart(sit->p,currRing,
             it->lcm,currRing,
             currRing->real_var_start,currRing->real_var_end))
          {
            if(TEST_OPT_DEBUG)
            {
               Print("chain-crit-part: S[]=");
               p_wrp(sit->p,currRing);
               Print(" divide B[i].lcm=");
               p_wrp(it->lcm,currRing);
               PrintLn();
            }
            it = strat_B(strat).erase(it);
            strat->c3++;
          }
          else
            ++it;
        }
      }
    }
    strat->S.clear_pairtest();
  }
  if (strat->Gebauer || strat->fromT)
  {
    if (strat->sugarCrit)
    {
    /*
    *suppose L[j] == (s,r) and p/lcm(s,r)
    *and lcm(s,r)#lcm(s,p) and lcm(s,r)#lcm(r,p)
    *and in case the sugar is o.k. then L[j] can be canceled
    */
      {
        for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
        {
          if (sugarDivisibleBy(ecart,it->ecart)
          && ((pNext(it->p) == strat->tail) || (rHasGlobalOrdering(currRing)))
          && pCompareChainPart(p,it->p1,it->p2,it->lcm)
          && (it->p == strat->tail))
          {
            if(TEST_OPT_DEBUG)
            {
               PrintS("chain-crit-part: pCompareChainPart p=");
               p_wrp(p,currRing);
               Print(" delete L");
               p_wrp(it->lcm,currRing);
               PrintLn();
            }
            it = strat->L.erase(it);
            strat->c3++;
          }
          else
            ++it;
        }
      }
      /*
      *this is GEBAUER-MOELLER:
      *in B all elements with the same lcm except the "best"
      *(i.e. the last one in B with this property) will be canceled
      */
      {
        const unsigned long* sev = strat_B(strat).sev_flat_data();
        const size_t n = strat_B(strat).sev_flat_size();
        for (size_t i = 0; i < n; i++)
        {
          if (sev[i] == 0) continue;
          for (size_t j = i + 1; j < n; j++)
          {
            if (sev[j] == 0) continue;
            LObject* a = strat_B(strat).flat_ptr(i);
            LObject* b = strat_B(strat).flat_ptr(j);
            if (a == NULL || b == NULL) continue;
            if (pLmEqual(a->lcm,b->lcm))
            {
              strat->c3++;
              if (a->ecart < b->ecart)
              {
                if(TEST_OPT_DEBUG)
                {
                  Print("chain-crit-part: sugar B[j].lcm=");
                  p_wrp(a->lcm,currRing);
                  Print(" delete B[i]");
                  p_wrp(b->lcm,currRing);
                  PrintLn();
                }
                strat_B(strat).erase(strat_B(strat).uiter_at(j));
              }
              else if (a->ecart > b->ecart)
              {
                if(TEST_OPT_DEBUG)
                {
                  Print("chain-crit-part: sugar B[i].lcm=");
                  p_wrp(b->lcm,currRing);
                  Print(" delete B[j]");
                  p_wrp(a->lcm,currRing);
                  PrintLn();
                }
                strat_B(strat).erase(strat_B(strat).uiter_at(i));
                break;  // i is gone, move to next i
              }
              else
              {
                // Equal ecart: use key_comp() tiebreaker for determinism
                if (strat_B(strat).key_comp()(*a, *b))
                  strat_B(strat).erase(strat_B(strat).uiter_at(j));
                else
                {
                  strat_B(strat).erase(strat_B(strat).uiter_at(i));
                  break;
                }
              }
            }
          }
        }
      }
    }
    else /*sugarCrit*/
    {
      /*
      *suppose L[j] == (s,r) and p/lcm(s,r)
      *and lcm(s,r)#lcm(s,p) and lcm(s,r)#lcm(r,p)
      *and in case the sugar is o.k. then L[j] can be canceled
      */
      {
        for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
        {
          if (pCompareChainPart(p,it->p1,it->p2,it->lcm)
          && ((pNext(it->p) == strat->tail)||(rHasGlobalOrdering(currRing))))
          {
            if(TEST_OPT_DEBUG)
            {
              PrintS("chain-crit-part: sugar:pCompareChainPart p=");
              p_wrp(p,currRing);
              Print(" delete L[j]");
              p_wrp(it->lcm,currRing);
              PrintLn();
            }
            it = strat->L.erase(it);
            strat->c3++;
          }
          else
            ++it;
        }
      }
      /*
      *this is GEBAUER-MOELLER:
      *in B all elements with the same lcm except the "best"
      *(i.e. the last one in B with this property) will be canceled
      */
      {
        const unsigned long* sev = strat_B(strat).sev_flat_data();
        const size_t n = strat_B(strat).sev_flat_size();
        for (size_t i = 0; i < n; i++)
        {
          if (sev[i] == 0) continue;
          for (size_t j = i + 1; j < n; j++)
          {
            if (sev[j] == 0) continue;
            LObject* a = strat_B(strat).flat_ptr(i);
            LObject* b = strat_B(strat).flat_ptr(j);
            if (a == NULL || b == NULL) continue;
            if (pLmEqual(a->lcm,b->lcm))
            {
              if(TEST_OPT_DEBUG)
              {
                Print("chain-crit-part: equal lcm B[j].lcm=");
                p_wrp(a->lcm,currRing);
                Print(" delete B[i]\n");
              }
              strat->c3++;
              // Erase the worse element; keep the one that sorts first
              if (strat_B(strat).key_comp()(*a, *b))
                strat_B(strat).erase(strat_B(strat).uiter_at(j));
              else
              {
                strat_B(strat).erase(strat_B(strat).uiter_at(i));
                break;  // i is gone
              }
            }
          }
        }
      }
    }
    /*
    *the elements of B enter L
    */
    kMergeBintoL(strat);
  }
  else
  {
    {
      for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
      {
        if (pCompareChainPart(p,it->p1,it->p2,it->lcm)
        && ((pNext(it->p) == strat->tail)||(rHasGlobalOrdering(currRing))))
        {
          if(TEST_OPT_DEBUG)
          {
            PrintS("chain-crit-part: pCompareChainPart p=");
            p_wrp(p,currRing);
            Print(" delete L[j]");
            p_wrp(it->lcm,currRing);
            PrintLn();
          }
          it = strat->L.erase(it);
          strat->c3++;
        }
        else
          ++it;
      }
    }
    /*
    *this is our MODIFICATION of GEBAUER-MOELLER:
    *Merge B into L, then deduplicate B-origin elements using
    *pair_index for O(1) triangle checks via bvec.
    *Uses _p_LmDivisibleByPart instead of pDivisibleBy.
    */
    std::vector<LSet::iterator> bvec = kMergeBintoL_and_return_iterators(strat);
    LSet::iterator endL = strat->L.end();
    for (size_t ji = 0; ji < bvec.size(); ji++)
    {
      if (bvec[ji] == endL) continue;
      for (size_t ii = ji + 1; ii < bvec.size(); ii++)
      {
        if (bvec[ii] == endL) continue;
        if ((bvec[ji]->sev_lcm == bvec[ii]->sev_lcm)
            && pLmEqual(bvec[ji]->lcm, bvec[ii]->lcm))
        {
          /*bvec[ii] could be canceled but we search for a better one to cancel*/
          strat->c3++;
          auto lt = bvec[ii] + 1;
          if (isInPairsetL(lt,bvec[ji]->p1,bvec[ii]->p1,strat)
          && (pNext(lt->p) == strat->tail)
          && (!pLmEqual(bvec[ii]->p,lt->p))
          && _p_LmDivisibleByPart(p,currRing,
                         lt->lcm,currRing,
                         currRing->real_var_start, currRing->real_var_end))
          {
            /*
            *"NOT equal(...)" because in case of "equal" the element L[l]
            *is "older" and has to be from theoretical point of view behind
            *L[i], but we do not want to reorder L
            */
            if(TEST_OPT_DEBUG)
            {
              PrintS("chain-crit-part: divisible_by p=");
              p_wrp(p,currRing);
              Print(" delete L[l]");
              p_wrp(lt->lcm,currRing);
              PrintLn();
            }
            strat->L.erase(lt);
            /*
            *L[l] will be canceled, we cannot cancel L[i] later on,
            *so we null it out in bvec to prevent re-processing
            */
            bvec[ii] = endL;
          }
          else
          {
            if(TEST_OPT_DEBUG)
            {
              PrintS("chain-crit-part: divisible_by(2) p=");
              p_wrp(p,currRing);
              Print(" delete L[i]");
              p_wrp(bvec[ii]->lcm,currRing);
              PrintLn();
            }
            strat->L.erase(bvec[ii]);
            bvec[ii] = endL;
          }
        }
      }
    }
  }
}
#endif

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L
*/
void initenterpairs (poly h,int k,int ecart,int isFromQ,kStrategy strat, int atR/* = -1*/)
{

  if ((strat->syzComp==0)
  || (pGetComp(h)<=strat->syzComp))
  {
    int j;
    BOOLEAN new_pair=FALSE;

    if (pGetComp(h)==0)
    {
      /* for Q!=NULL: build pairs (f,q),(f1,f2), but not (q1,q2)*/
      if ((isFromQ)&&(strat->hasFromQ))
      {
        for (auto sit=strat->S.begin(); sit!=strat->S.end(); ++sit)
        {
          if (!arrival_id_ok(*sit)) continue; // skip h and peer-drain survivors
          if (!sit->fromQ)
          {
            new_pair=TRUE;
            strat->enterOnePair(*sit,h,ecart,isFromQ,strat, atR);
          //Print("j:%d, L.size():%d\n",sit.index(),(int)strat->L.size());
          }
        }
      }
      else
      {
        new_pair=TRUE;
        for (auto sit=strat->S.begin(); sit!=strat->S.end(); ++sit)
        {
          if (!arrival_id_ok(*sit)) continue; // skip h and peer-drain survivors
          strat->enterOnePair(*sit,h,ecart,isFromQ,strat, atR);
          //Print("j:%d, L.size():%d\n",sit.index(),(int)strat->L.size());
        }
      }
    }
    else
    {
      for (auto sit=strat->S.begin(); sit!=strat->S.end(); ++sit)
      {
        if (!arrival_id_ok(*sit)) continue; // skip h and peer-drain survivors
        if ((pGetComp(h)==pGetComp(sit->p))
        || (pGetComp(sit->p)==0))
        {
          new_pair=TRUE;
          strat->enterOnePair(*sit,h,ecart,isFromQ,strat, atR);
        //Print("j:%d, Ll:%d\n",sit.index(),strat->Ll);
        }
      }
    }
    if (new_pair)
    {
    #ifdef HAVE_RATGRING
      if (currRing->real_var_start>0)
        chainCritPart(h,ecart,strat);
      else
    #endif
      strat->chainCrit(h,ecart,strat);
    }
    kMergeBintoL(strat);
  }
}

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L
*using signatures <= only for signature-based standard basis algorithms
*/

void initenterpairsSig (poly h,poly hSig,int hFrom,int k,int ecart,int isFromQ,kStrategy strat, int atR = -1)
{

  if ((strat->syzComp==0)
  || (pGetComp(h)<=strat->syzComp))
  {
    int j;
    BOOLEAN new_pair=FALSE;

    if (pGetComp(h)==0)
    {
      /* for Q!=NULL: build pairs (f,q),(f1,f2), but not (q1,q2)*/
      if ((isFromQ)&&(strat->hasFromQ))
      {
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if (!sjt->fromQ)
          {
            new_pair=TRUE;
            enterOnePairSig(*sjt,sjt,h,hSig,hFrom,ecart,isFromQ,strat, atR);
          //Print("j:%d, L.size():%d\n",j,(int)strat->L.size());
          }
        }
      }
      else
      {
        new_pair=TRUE;
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          enterOnePairSig(*sjt,sjt,h,hSig,hFrom,ecart,isFromQ,strat, atR);
          //Print("j:%d, L.size():%d\n",j,(int)strat->L.size());
        }
      }
    }
    else
    {
      for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
      {
        if ((pGetComp(h)==pGetComp(sjt->p))
        || (pGetComp(sjt->p)==0))
        {
          new_pair=TRUE;
          enterOnePairSig(*sjt,sjt,h,hSig,hFrom,ecart,isFromQ,strat, atR);
        //Print("j:%d, Ll:%d\n",j,strat->Ll);
        }
      }
    }

    if (new_pair)
    {
#ifdef HAVE_RATGRING
      if (currRing->real_var_start>0)
        chainCritPart(h,ecart,strat);
      else
#endif
      strat->chainCrit(h,ecart,strat);
    }
  }
}

void initenterpairsSigRing (poly h,poly hSig,int hFrom,int k,int ecart,int isFromQ,kStrategy strat, int atR = -1)
{

  if ((strat->syzComp==0)
  || (pGetComp(h)<=strat->syzComp))
  {
    int j;

    if (pGetComp(h)==0)
    {
      /* for Q!=NULL: build pairs (f,q),(f1,f2), but not (q1,q2)*/
      if ((isFromQ)&&(strat->hasFromQ))
      {
        for (auto sjt = strat->S.begin(); sjt != strat->S.end() && !strat->sigdrop; ++sjt)
        {
          if (!sjt->fromQ)
          {
            enterOnePairSigRing(*sjt,h,hSig,hFrom,ecart,isFromQ,strat, atR);
          //Print("j:%d, L.size():%d\n",j,(int)strat->L.size());
          }
        }
      }
      else
      {
        for (auto sjt = strat->S.begin(); sjt != strat->S.end() && !strat->sigdrop; ++sjt)
        {
          enterOnePairSigRing(*sjt,h,hSig,hFrom,ecart,isFromQ,strat, atR);
          //Print("j:%d, L.size():%d\n",j,(int)strat->L.size());
        }
      }
    }
    else
    {
      for (auto sjt = strat->S.begin(); sjt != strat->S.end() && !strat->sigdrop; ++sjt)
      {
        if ((pGetComp(h)==pGetComp(sjt->p))
        || (pGetComp(sjt->p)==0))
        {
          enterOnePairSigRing(*sjt,h,hSig,hFrom,ecart,isFromQ,strat, atR);
        //Print("j:%d, Ll:%d\n",j,strat->Ll);
        }
      }
    }

#if 0
    if (new_pair)
    {
#ifdef HAVE_RATGRING
      if (currRing->real_var_start>0)
        chainCritPart(h,ecart,strat);
      else
#endif
      strat->chainCrit(h,ecart,strat);
    }
#endif
  }
}

/*2
*the pairset B of pairs of type (s[i],p) is complete now. It will be updated
*using the chain-criterion in B and L and enters B to L
*/
void chainCritRing (poly p,int, kStrategy strat)
{
  int j;
  /*
  *pairtest[i] is TRUE if spoly(S[i],p) == 0.
  *In this case all elements in B such
  *that their lcm is divisible by the leading term of S[i] can be canceled
  *
  * Parallel phase-1: iterate thread-local hits.
  */
  if (t_local_pairtest_hits != NULL)
  {
    for (SElement *sit : *t_local_pairtest_hits)
    {
      for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
      {
        if (pDivisibleBy(sit->p, it->lcm)
            && n_DivBy(pGetCoeff(it->lcm), pGetCoeff(sit->p), currRing->cf))
        {
          it = strat_B(strat).erase(it);
          strat->c3++;
        }
        else ++it;
      }
    }
    t_local_pairtest_hits->clear();
  }
  else if (strat->S.has_pairtest())
  {
    {
      /*- i.e. there is an i with pairtest[i]==TRUE -*/
      for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
      {
        if (selement_pairtest_load(*sit))
        {
          for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); )
          {
            if (pDivisibleBy(sit->p,it->lcm) && n_DivBy(pGetCoeff(it->lcm), pGetCoeff(sit->p),currRing->cf))
            {
#ifdef KDEBUG
              if (TEST_OPT_DEBUG)
              {
                PrintS("--- chain criterion func chainCritRing type 1\n");
                PrintS("S[].p:");
                wrp(sit->p);
                PrintS("  strat_B(strat)[i].lcm:");
                wrp(it->lcm);PrintLn();
                pWrite(it->p);
                pWrite(it->p1);
                pWrite(it->p2);
                wrp(it->lcm);
                PrintLn();
              }
#endif
              it = strat_B(strat).erase(it);
              strat->c3++;
            }
            else
              ++it;
          }
        }
      }
    }
    strat->S.clear_pairtest();
  }
  assume(!(strat->Gebauer || strat->fromT));
  {
    unsigned long sev_p = p_GetShortExpVector(p, currRing);
    for (auto it = strat->L.ufbegin_lcm(sev_p); it != strat->L.ufend_lcm(); )
    {
      if ((it->lcm != NULL) && n_DivBy(pGetCoeff(it->lcm), pGetCoeff(p), currRing->cf)
      && pCompareChain(p,it->p1,it->p2,it->lcm)
      && ((pNext(it->p) == strat->tail)||(rHasGlobalOrdering(currRing))))
      {
        it = strat->L.erase(it);
        strat->c3++;
#ifdef KDEBUG
        if (TEST_OPT_DEBUG)
        {
          PrintS("--- chain criterion func chainCritRing type 2\n");
          PrintS("  p:");
          wrp(p);
          PrintLn();
        }
#endif
      }
      else
        ++it;
    }
  }
  /*
  *this is our MODIFICATION of GEBAUER-MOELLER:
  *Merge B into L, then deduplicate the B-origin elements using
  *pair_index for O(1) triangle checks.  kMergeBintoL_and_return_iterators
  *gives us iterators to the new elements in L-order, avoiding a full
  *scan of L.
  *
  *The outer loop iterates worst-to-first (ji = high index = worse position).
  *For each pair with the same LCM, the inner element ii (better position)
  *gets canceled, matching the original Ring semantics where the worse-
  *positioned element survives.  The n_DivBy coefficient check is preserved.
  */
  std::vector<LSet::iterator> bvec = kMergeBintoL_and_return_iterators(strat);
  LSet::iterator endL = strat->L.end();
  /* Iterate worst-to-first: bvec is sorted best(0) to worst(size-1),
   * so ji starts at the end and works backward */
  for (int ji = (int)bvec.size() - 1; ji >= 0; ji--)
  {
    if (bvec[ji] == endL) continue;
    for (int ii = ji - 1; ii >= 0; ii--)
    {
      if (bvec[ii] == endL) continue;
      /* bvec[ji] is worse (outer), bvec[ii] is better (inner) */
      if (n_DivBy(pGetCoeff(bvec[ji]->lcm), pGetCoeff(bvec[ii]->lcm), currRing->cf)
          && pLmEqual(bvec[ji]->lcm, bvec[ii]->lcm))
      {
        /*bvec[ii] could be canceled but we search for a better one to cancel*/
        strat->c3++;
#ifdef KDEBUG
        if (TEST_OPT_DEBUG)
        {
          PrintS("--- chain criterion func chainCritRing type 3\n");
          PrintS("strat->L[j].lcm:");
          wrp(bvec[ji]->lcm);
          PrintS("  strat->L[i].lcm:");
          wrp(bvec[ii]->lcm);
          PrintLn();
        }
#endif
        auto lt = bvec[ii] + 1;
        if (isInPairsetL(lt,bvec[ji]->p1,bvec[ii]->p1,strat)
        && (pNext(lt->p) == strat->tail)
        && (!pLmEqual(bvec[ii]->p,lt->p))
        && pDivisibleBy(p,lt->lcm))
        {
          /*
          *"NOT equal(...)" because in case of "equal" the element L[l]
          *is "older" and has to be from theoretical point of view behind
          *L[i], but we do not want to reorder L
          */
          strat->L.erase(lt);
          /*
          *L[l] will be canceled, we cannot cancel L[i] later on,
          *so we null it out in bvec to prevent re-processing
          */
          bvec[ii] = endL;
        }
        else
        {
          strat->L.erase(bvec[ii]);
          bvec[ii] = endL;
        }
      }
    }
  }
}

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L
*/
void initenterstrongPairs (poly h,int k,int ecart,int isFromQ,kStrategy strat, int atR = -1)
{
  if (!nIsOne(pGetCoeff(h)))
  {
    int j;
    BOOLEAN new_pair=FALSE;

    if (pGetComp(h)==0)
    {
      /* for Q!=NULL: build pairs (f,q),(f1,f2), but not (q1,q2)*/
      if ((isFromQ)&&(strat->hasFromQ))
      {
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if (!arrival_id_ok(*sjt)) continue; // parallel phase-1 filter
          if (!sjt->fromQ)
          {
            new_pair=TRUE;
            enterOneStrongPoly(*sjt,h,ecart,isFromQ,strat, atR, FALSE);
          }
        }
      }
      else
      {
        new_pair=TRUE;
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if (!arrival_id_ok(*sjt)) continue; // parallel phase-1 filter
          enterOneStrongPoly(*sjt,h,ecart,isFromQ,strat, atR, FALSE);
        }
      }
    }
    else
    {
      for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
      {
        if (!arrival_id_ok(*sjt)) continue; // parallel phase-1 filter
        if ((pGetComp(h)==pGetComp(sjt->p))
        || (pGetComp(sjt->p)==0))
        {
          new_pair=TRUE;
          enterOneStrongPoly(*sjt,h,ecart,isFromQ,strat, atR, FALSE);
        }
      }
    }
    if (new_pair)
    {
    #ifdef HAVE_RATGRING
      if (currRing->real_var_start>0)
        chainCritPart(h,ecart,strat);
      else
    #endif
      strat->chainCrit(h,ecart,strat);
    }
    kMergeBintoL(strat);
  }
}

static void initenterstrongPairsSig (poly h,poly hSig, int k,int ecart,int isFromQ,kStrategy strat, int atR = -1)
{
  const int iCompH = pGetComp(h);
  if (!nIsOne(pGetCoeff(h)))
  {
    int j;

    for (auto sjt = strat->S.begin(); sjt != strat->S.end() && !strat->sigdrop; ++sjt)
    {
      // Print("j:%d, Ll:%d\n",j,strat->Ll);
//      if (((unsigned long) pGetCoeff(h) % (unsigned long) pGetCoeff(strat->S[j].p) != 0) &&
//         ((unsigned long) pGetCoeff(strat->S[j].p) % (unsigned long) pGetCoeff(h) != 0))
      if (((iCompH == pGetComp(sjt->p))
      || (0 == pGetComp(sjt->p)))
      && ((iCompH<=strat->syzComp)||(strat->syzComp==0)))
      {
        enterOneStrongPolySig(*sjt,h,hSig,ecart,isFromQ,strat, atR);
      }
    }
  }
}

/*2
* Generates spoly(0, h) if applicable. Assumes ring has zero divisors
*/
void enterExtendedSpoly(poly h,kStrategy strat)
{
  if (nIsOne(pGetCoeff(h))) return;
  number gcd;
  number zero=n_Init(0,currRing->cf);
  bool go = false;
  if (n_DivBy(zero, pGetCoeff(h), currRing->cf))
  {
    gcd = n_Ann(pGetCoeff(h),currRing->cf);
    go = true;
  }
  else
    gcd = n_Gcd(zero, pGetCoeff(h), strat->tailRing->cf);
  if (go || !nIsOne(gcd))
  {
    poly p = h->next;
    if (!go)
    {
      number tmp = gcd;
      gcd = n_Ann(gcd,currRing->cf);
      nDelete(&tmp);
    }
    p_Test(p,strat->tailRing);
    p = __pp_Mult_nn(p, gcd, strat->tailRing);

    if (p != NULL)
    {
      if (TEST_OPT_PROT)
      {
        PrintS("Z");
      }
#ifdef KDEBUG
      if (TEST_OPT_DEBUG)
      {
        PrintS("--- create zero spoly: ");
        p_wrp(h,currRing,strat->tailRing);
        PrintS(" ---> ");
      }
#endif
      poly tmp = pInit();
      pSetCoeff0(tmp, pGetCoeff(p));
      for (int i = 1; i <= rVar(currRing); i++)
      {
        pSetExp(tmp, i, p_GetExp(p, i, strat->tailRing));
      }
      if (rRing_has_Comp(currRing) && rRing_has_Comp(strat->tailRing))
      {
        p_SetComp(tmp, __p_GetComp(p, strat->tailRing), currRing);
      }
      p_Setm(tmp, currRing);
      p = p_LmFreeAndNext(p, strat->tailRing);
      pNext(tmp) = p;
      LObject Lp;
      Lp.Init();
      Lp.p = tmp;
      Lp.tailRing = strat->tailRing;
      if (Lp.p!=NULL)
      {
        strat->initEcart(&Lp);
        Lp.sev = pGetShortExpVector(Lp.p);
        if (strat->tailRing != currRing)
        {
          Lp.t_p = k_LmInit_currRing_2_tailRing(Lp.p, strat->tailRing);
        }
#ifdef KDEBUG
        if (TEST_OPT_DEBUG)
        {
          p_wrp(tmp,currRing,strat->tailRing);
          PrintLn();
        }
#endif
        strat->L.push(Lp);
      }
    }
  }
  nDelete(&zero);
  nDelete(&gcd);
}

void enterExtendedSpolySig(poly h,poly hSig,kStrategy strat)
{
  if (nIsOne(pGetCoeff(h))) return;
  number gcd;
  number zero=n_Init(0,currRing->cf);
  bool go = false;
  if (n_DivBy(zero, pGetCoeff(h), currRing->cf))
  {
    gcd = n_Ann(pGetCoeff(h),currRing->cf);
    go = true;
  }
  else
    gcd = n_Gcd(zero, pGetCoeff(h), strat->tailRing->cf);
  if (go || !nIsOne(gcd))
  {
    poly p = h->next;
    if (!go)
    {
      number tmp = gcd;
      gcd = n_Ann(gcd,currRing->cf);
      nDelete(&tmp);
    }
    p_Test(p,strat->tailRing);
    p = __pp_Mult_nn(p, gcd, strat->tailRing);

    if (p != NULL)
    {
      if (TEST_OPT_PROT)
      {
        PrintS("Z");
      }
#ifdef KDEBUG
      if (TEST_OPT_DEBUG)
      {
        PrintS("--- create zero spoly: ");
        p_wrp(h,currRing,strat->tailRing);
        PrintS(" ---> ");
      }
#endif
      poly tmp = pInit();
      pSetCoeff0(tmp, pGetCoeff(p));
      for (int i = 1; i <= rVar(currRing); i++)
      {
        pSetExp(tmp, i, p_GetExp(p, i, strat->tailRing));
      }
      if (rRing_has_Comp(currRing) && rRing_has_Comp(strat->tailRing))
      {
        p_SetComp(tmp, __p_GetComp(p, strat->tailRing), currRing);
      }
      p_Setm(tmp, currRing);
      p = p_LmFreeAndNext(p, strat->tailRing);
      pNext(tmp) = p;
      LObject Lp;
      Lp.Init();
      Lp.p = tmp;
      //printf("\nOld\n");pWrite(h);pWrite(hSig);
      #if EXT_POLY_NEW
      Lp.sig = __pp_Mult_nn(hSig, gcd, currRing);
      if(Lp.sig == NULL || nIsZero(pGetCoeff(Lp.sig)))
      {
        strat->sigdrop = TRUE;
        //Try to reduce it as far as we can via redRing
        int red_result = redRing(&Lp,strat);
        if(red_result == 0)
        {
          // Cancel the sigdrop
          p_Delete(&Lp.sig,currRing);Lp.sig = NULL;
          strat->sigdrop = FALSE;
        }
        else
        {
          strat->enterS(strat->P, strat, strat->T.size(), strat->S.end());
          #if 1
          strat->enterS(Lp, strat, strat->T.size()-1, strat->S.begin());
          #endif
        }
        nDelete(&zero);
        nDelete(&gcd);
        return;
      }
      #else
      Lp.sig = pOne();
      if(!strat->L.empty())
        p_SetComp(Lp.sig,pGetComp(strat->L.end()->sig)+1,currRing);
      else
        p_SetComp(Lp.sig,pGetComp(hSig)+1,currRing);
      #endif
      Lp.tailRing = strat->tailRing;
      if (Lp.p!=NULL)
      {
        strat->initEcart(&Lp);
        Lp.sev = pGetShortExpVector(Lp.p);
        if (strat->tailRing != currRing)
        {
          Lp.t_p = k_LmInit_currRing_2_tailRing(Lp.p, strat->tailRing);
        }
#ifdef KDEBUG
        if (TEST_OPT_DEBUG)
        {
          p_wrp(tmp,currRing,strat->tailRing);
          PrintLn();
        }
#endif
  //pWrite(h);pWrite(hSig);pWrite(Lp.p);pWrite(Lp.sig);printf("\n------------------\n");getchar();
        strat->L.push(Lp);
      }
    }
  }
  nDelete(&gcd);
  nDelete(&zero);
}

void clearSbatch (poly h,int k,sBasisSet::iterator pos,kStrategy strat)
{
  if ( (!strat->fromT)
  && ((strat->syzComp==0)
    ||(pGetComp(h)<=strat->syzComp)
  ))
  {
    // In serial mode, enterS has not yet placed h in S when clearSbatch is
    // called, so k = size() - 1 holds.  In parallel phase-1 mode, enterS
    // runs in phase 0 (before clearSbatch), so S already contains h; the
    // arrival_id filter skips h (and concurrent peer survivors) so the
    // iteration bound is conservative.
    unsigned long h_sev = pGetShortExpVector(h);
    for (auto it = pos; it != strat->S.end(); ++it)
    {
      if (!arrival_id_ok(*it)) continue; // skip h and peer-drain survivors
      strat->S.clear_if_divisible(h, h_sev, it, strat);
    }
  }
}

/*2
* Generates a sufficient set of spolys (maybe just a finite generating
* set of the syzygys)
*/
void superenterpairs (poly h,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR)
{
  assume (rField_is_Ring(currRing));
#if HAVE_SHIFTBBA
  assume(!rIsLPRing(currRing)); /* LP should use enterpairsShift */
#endif
  // enter also zero divisor * poly, if this is non zero and of smaller degree
  if (!(rField_is_Domain(currRing))) enterExtendedSpoly(h, strat);
  initenterstrongPairs(h, k, ecart, 0, strat, atR);
  initenterpairs(h, k, ecart, 0, strat, atR);
  clearSbatch(h, k, pos, strat);
}

void superenterpairsSig (poly h,poly hSig,int hFrom,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR)
{
  assume (rField_is_Ring(currRing));
  // enter also zero divisor * poly, if this is non zero and of smaller degree
  if (!(rField_is_Domain(currRing))) enterExtendedSpolySig(h, hSig, strat);
  if(strat->sigdrop) return;
  initenterpairsSigRing(h, hSig, hFrom, k, ecart, 0, strat, atR);
  if(strat->sigdrop) return;
  initenterstrongPairsSig(h, hSig, k, ecart, 0, strat, atR);
  if(strat->sigdrop) return;
  clearSbatch(h, k, pos, strat);
}

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L(via initenterpairs)
*superfluous elements in S will be deleted
*/
void enterpairs (poly h,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR)
{
  assume (!rField_is_Ring(currRing));
  initenterpairs(h,k,ecart,0,strat, atR);
  if ( (!strat->fromT)
  && ((strat->syzComp==0)
    ||(pGetComp(h)<=strat->syzComp)))
  {
    // See clearSbatch comment: in parallel phase-1 mode, h is already in S
    // at `pos`; arrival_id filter skips it.  k may not equal
    // strat->S.size() - 1 in parallel mode.
    unsigned long h_sev = pGetShortExpVector(h);
    for (auto it = pos; it != strat->S.end(); ++it)
    {
      if (!arrival_id_ok(*it)) continue; // skip h and peer-drain survivors
      strat->S.clear_if_divisible(h, h_sev, it, strat);
    }
  }
}

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L(via initenterpairs)
*superfluous elements in S will be deleted
*this is a special variant of signature-based algorithms including the
*signatures for criteria checks
*/
void enterpairsSig (poly h,poly hSig,int hFrom,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR)
{
  assume (!rField_is_Ring(currRing));
  initenterpairsSig(h,hSig,hFrom,k,ecart,0,strat, atR);
  if ( (!strat->fromT)
  && ((strat->syzComp==0)
    ||(pGetComp(h)<=strat->syzComp)))
  {
    // sba() is serial (not reached from the parallel drain); arrival_id
    // filter is a no-op here but kept for uniformity.
    unsigned long h_sev = pGetShortExpVector(h);
    for (auto it = pos; it != strat->S.end(); ++it)
    {
      if (!arrival_id_ok(*it)) continue;
      strat->S.clear_if_divisible(h, h_sev, it, strat);
    }
  }
}

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L(via initenterpairs)
*superfluous elements in S will be deleted
*/
void enterpairsSpecial (poly h,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR = -1)
{
  const int iCompH = pGetComp(h);

  if (rField_is_Ring(currRing))
  {
    for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
    {
      const int iCompSj = pGetComp(sjt->p);
      if ((iCompH==iCompSj)
          //|| (0==iCompH) // can only happen,if iCompSj==0
          || (0==iCompSj))
      {
        enterOnePairRing(*sjt,h,ecart,FALSE,strat, atR);
      }
    }
    kMergeBintoL(strat);
  }
  else
  {
    for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
    {
      const int iCompSj = pGetComp(sjt->p);
      if ((iCompH==iCompSj)
          //|| (0==iCompH) // can only happen,if iCompSj==0
          || (0==iCompSj))
      {
        enterOnePairSpecial(*sjt,h,ecart,strat, atR);
      }
    }
  }

  if (strat->noClearS) return;

  {
    // Caller contract: k == strat->S.size()-1 on entry.
    assume(k == strat->S.size() - 1);
    unsigned long h_sev = pGetShortExpVector(h);
    for (auto it = pos; it != strat->S.end(); ++it)
      strat->S.clear_if_divisible(h, h_sev, it, strat);
  }
}

/*2
*reorders  s with respect to posInS,
*suc is the first changed index or zero
*/

/*2
* sBasisSet::reorder — reorder S with respect to find_pos.
* Replaces the old free reorderS function.
*/
void sBasisSet::reorder(int *suc, kStrategy strat)
{
  int i,j,at;

  // Compact first: reorder walks by physical index (elem(i)) and
  // shifts whole SElements. Tombstones in the array would corrupt the
  // shift semantics. After compact(), physical count == live_count,
  // and *suc (a live-index set by a prior reorder call or 0 initially)
  // is a valid physical index into the compacted sequence.
  compact();

  int new_suc = size();
  i = *suc;
  if (i < 0) i = 0;

  for (; i < size(); i++)
  {
    // Search within the first i elements for the new position of elem(i).
    at = find_pos(elem(i).p, elem(i).ecart, iterator(this, i)).index();
    if (at != i)
    {
      if (new_suc > at) new_suc = at;
      {
        // Shift entire SElements so all fields stay in sync.
        SElement tmp = elem(i);
        for (j=i; j>=at+1; j--)
          elem(j) = elem(j-1);
        elem(at) = tmp;
      }
    }
  }
  if (new_suc < size()) *suc = new_suc;
  else                   *suc = -1;
}



/*2
* sBasisSet::find_pos — binary search for sorted insertion position.
* Replaces the old free posInS function. Public form scans all of S.
*
* Tombstone handling: find_pos does index-based binary search, which
* cannot cope with tombstoned entries (their .p may have been set to
* NULL by callers before the erase, e.g. updateS/redBba). Compact the
* array first if any tombstones exist. This is amortized O(n) per
* find_pos call only when tombstones accumulated since the last
* compact; insert-heavy phases pay no compaction cost.
*/
sBasisSet::iterator sBasisSet::find_pos(const poly p, int ecart_p)
{
  if (deleted_count_ > 0) compact();
  return find_pos(p, ecart_p, end());
}

/* Bounded form: search within [0 .. end_bound). Only used internally
 * by find_divisor_search_bound (which must ensure the set is already
 * compacted). */
sBasisSet::iterator sBasisSet::find_pos(const poly p, int ecart_p, iterator end_bound)
{
  if (empty()) return begin();
  int length = end_bound.index() - 1;
  if (length < 0) return begin();
  int i;
  int an = 0;
  int en = length;
  int cmp_int = currRing->OrdSgn;
  if ((rHasMixedOrdering(currRing))
#ifdef HAVE_PLURAL
  && (currRing->real_var_start==0)
#endif
#if 0
  || ((strat->ak>0) && ((currRing->order[0]==ringorder_c)||((currRing->order[0]==ringorder_C))))
#endif
  )
  {
    int o=p_Deg(p,currRing);
    int oo=p_Deg(elem(length).p,currRing);

    if ((oo<o)
    || ((o==oo) && (pLmCmp(elem(length).p,p)!= cmp_int)))
      return iterator(this, length+1);

    loop
    {
      if (an >= en-1)
      {
        if ((p_Deg(elem(an).p,currRing)>=o) && (pLmCmp(elem(an).p,p) == cmp_int))
        {
          return iterator(this, an);
        }
        return iterator(this, en);
      }
      i=(an+en) / 2;
      if ((p_Deg(elem(i).p,currRing)>=o) && (pLmCmp(elem(i).p,p) == cmp_int)) en=i;
      else                              an=i;
    }
  }
  else
  {
    if (rField_is_Ring(currRing))
    {
      if (pLmCmp(elem(length).p,p)== -cmp_int)
        return iterator(this, length+1);
      int cmp;
      loop
      {
        if (an >= en-1)
        {
          cmp = pLmCmp(elem(an).p,p);
          if (cmp == cmp_int)  return iterator(this, an);
          if (cmp == -cmp_int) return iterator(this, en);
          if (n_DivBy(pGetCoeff(p), pGetCoeff(elem(an).p), currRing->cf)) return iterator(this, en);
          return iterator(this, an);
        }
        i = (an+en) / 2;
        cmp = pLmCmp(elem(i).p,p);
        if (cmp == cmp_int)         en = i;
        else if (cmp == -cmp_int)   an = i;
        else
        {
          if (n_DivBy(pGetCoeff(p), pGetCoeff(elem(i).p), currRing->cf)) an = i;
          else en = i;
        }
      }
    }
    else
    if (pLmCmp(elem(length).p,p)== -cmp_int)
      return iterator(this, length+1);

    loop
    {
      if (an >= en-1)
      {
        if (pLmCmp(elem(an).p,p) == cmp_int) return iterator(this, an);
        if (pLmCmp(elem(an).p,p) == -cmp_int) return iterator(this, en);
        if ((cmp_int!=1)
        && ((elem(an).ecart)>ecart_p))
          return iterator(this, an);
        return iterator(this, en);
      }
      i=(an+en) / 2;
      if (pLmCmp(elem(i).p,p) == cmp_int) en=i;
      else if (pLmCmp(elem(i).p,p) == -cmp_int) an=i;
      else
      {
        if ((cmp_int!=1)
        &&((elem(i).ecart)<ecart_p))
          en=i;
        else
          an=i;
      }
    }
  }
}

/* sBasisSet::find_divisor_search_bound — encapsulates the search-range
 * narrowing used by kFindDivisibleByInS / kFindDivisibleByInS_noCF.
 *
 * For Ring/component/lex orderings, narrowing is not safe and we return
 * max unchanged. Otherwise, we use the bounded find_pos with ecart_p=0
 * to compute the position p would occupy if inserted, then clamp to
 * max. The caller uses the returned iterator as an exclusive upper
 * bound on the divisor scan.
 */
sBasisSet::iterator sBasisSet::find_divisor_search_bound(poly p, iterator max, kStrategy strat)
{
  if (rField_is_Ring(currRing) || (strat->ak > 0) || currRing->pLexOrder)
    return max;
  // Bounded find_pos returns the insertion point within [0, max). The
  // old code used an inclusive "ende = fp + 1" clamp; translated to an
  // exclusive iterator-end, that's fp + 2 clamped to max.
  // Compact first: the bounded find_pos does not compact itself so that
  // callers passing a pre-computed max iterator aren't silently affected.
  if (deleted_count_ > 0) {
    int max_idx = max.index();
    compact();
    // After compact, max_idx may now be past physical count. Re-clamp.
    if (max_idx > count) max_idx = count;
    max = iterator(this, max_idx);
  }
  iterator ende = find_pos(p, 0, max);
  // Advance two raw slots: one for the "+1" inclusive -> exclusive shift
  // (we want to include the element at fp), one more for the historic
  // "+1" overshoot. Cap at max.
  if (ende.index() + 2 <= max.index())
    ende = iterator(this, ende.index() + 2);
  else
    ende = max;
  return ende;
}

/*
* sBasisSet::find_pos_monfirst — binary search for monfirst insertion position.
* Replaces the old free posInSMonFirst function. Scans all of S.
*/
sBasisSet::iterator sBasisSet::find_pos_monfirst(const poly p)
{
  if (deleted_count_ > 0) compact();
  if (empty()) return begin();
  int length = size() - 1;
  if (pNext(p) == NULL)
  {
    // p is a monomial — insert among the leading monomials
    int mon = 0;
    for (int i = 0; i <= length; i++)
    {
      if (elem(i).p != NULL && pNext(elem(i).p) == NULL)
        mon++;
    }
    int o = p_Deg(p, currRing);
    int op = p_Deg(elem(mon).p, currRing);

    if ((op < o)
    || ((op == o) && (pLtCmp(elem(mon).p, p) == -1)))
      return iterator(this, length + 1);
    int i;
    int an = 0;
    int en = mon;
    loop
    {
      if (an >= en - 1)
      {
        op = p_Deg(elem(an).p, currRing);
        if ((op < o)
        || ((op == o) && (pLtCmp(elem(an).p, p) == -1)))
          return iterator(this, en);
        return iterator(this, an);
      }
      i = (an + en) / 2;
      op = p_Deg(elem(i).p, currRing);
      if ((op < o)
      || ((op == o) && (pLtCmp(elem(i).p, p) == -1)))
        an = i;
      else
        en = i;
    }
  }
  else /*if(pNext(p) != NULL)*/
  {
    int o = p_Deg(p, currRing);
    int op = p_Deg(elem(length).p, currRing);

    if ((op < o)
    || ((op == o) && (pLtCmp(elem(length).p, p) == -1)))
      return iterator(this, length + 1);
    int i;
    int an = 0;
    for (i = 0; i <= length; i++)
      if (elem(i).p != NULL && pNext(elem(i).p) == NULL)
        an++;
    int en = length;
    loop
    {
      if (an >= en - 1)
      {
        op = p_Deg(elem(an).p, currRing);
        if ((op < o)
        || ((op == o) && (pLtCmp(elem(an).p, p) == -1)))
          return iterator(this, en);
        return iterator(this, an);
      }
      i = (an + en) / 2;
      op = p_Deg(elem(i).p, currRing);
      if ((op < o)
      || ((op == o) && (pLtCmp(elem(i).p, p) == -1)))
        an = i;
      else
        en = i;
    }
  }
}

/*
* sBasisSet::insert — insert at position determined by ordering mode.
*/
sBasisSet::iterator sBasisSet::insert(const SElement& val, kStrategy strat)
{
  ensure_capacity(count + 1);
  switch (order_) {
  case SORDER_APPEND:
    BlockArray<SElement>::push_back(val);
    live_count_++;
    return iterator(this, count - 1);
  case SORDER_MONFIRST: {
    iterator pos = find_pos_monfirst(val.p);
    return insert_at(pos.index(), val);
  }
  case SORDER_STANDARD:
  default: {
    iterator pos = find_pos(val.p, val.ecart);
    return insert_at(pos.index(), val);
  }
  }
}


// sorts by degree and pLtCmp in the block between start,end;
// but puts pure monomials at the beginning
int posInIdealMonFirst (const ideal F, const poly p,int start,int end)
{
  if(end < 0 || end >= IDELEMS(F))
    end = IDELEMS(F);
  if (end<0) return 0;
  if(pNext(p) == NULL) return start;
  polyset set=F->m;
  int o = p_Deg(p,currRing);
  int op;
  int i;
  int an = start;
  for(i=start;i<end;i++)
    if(set[i] != NULL && pNext(set[i]) == NULL)
      an++;
  if(an == end-1)
    return end;
  int en= end;
  loop
  {
    if(an>=en)
      return en;
    if (an == en-1)
    {
      op = p_Deg(set[an],currRing);
      if ((op < o)
      || ((op == o) && (pLtCmp(set[an],p) == -1)))
        return en;
      return an;
    }
    i=(an+en) / 2;
    op = p_Deg(set[i],currRing);
    if ((op < o)
    || ((op == o) && (pLtCmp(set[i],p) == -1)))
      an=i;
    else
      en=i;
  }
}


/*2
* looks up the position of p in set
* the position is the last one
*/
int posInT0 (const BlockArray<TObject> &,const int length,LObject &)
{
  return (length+1);
}


/*2
* looks up the position of p in T
* set[0] is the smallest with respect to the ordering-procedure
* pComp
*/
int posInT1 (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  if (pLmCmp(set[length].p,p.p)!= currRing->OrdSgn) return length+1;

  int i;
  int an = 0;
  int en= length;
  int cmp_int=currRing->OrdSgn;

  loop
  {
    if (an >= en-1)
    {
      if (pLmCmp(set[an].p,p.p) == cmp_int) return an;
      return en;
    }
    i=(an+en) / 2;
    if (pLmCmp(set[i].p,p.p) == cmp_int) en=i;
    else                                 an=i;
  }
}

/*2
* looks up the position of p in T
* set[0] is the smallest with respect to the ordering-procedure
* length
*/
int posInT2 (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;
  p.GetpLength();
  if (set[length].length<p.length) return length+1;

  int i;
  int an = 0;
  int en= length;

  loop
  {
    if (an >= en-1)
    {
      if (set[an].length>p.length) return an;
      return en;
    }
    i=(an+en) / 2;
    if (set[i].length>p.length) en=i;
    else                        an=i;
  }
}

/*2
* looks up the position of p in T
* set[0] is the smallest with respect to the ordering-procedure
* totaldegree,pComp
*/
int posInT11 (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int o = p.GetpFDeg();
  int op = set[length].GetpFDeg();
  int cmp_int=currRing->OrdSgn;

  if ((op < o)
  || ((op == o) && (pLmCmp(set[length].p,p.p) != cmp_int)))
    return length+1;

  int i;
  int an = 0;
  int en= length;

  loop
  {
    if (an >= en-1)
    {
      op= set[an].GetpFDeg();
      if ((op > o)
      || (( op == o) && (pLmCmp(set[an].p,p.p) == cmp_int)))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg();
    if (( op > o)
    || (( op == o) && (pLmCmp(set[i].p,p.p) == cmp_int)))
      en=i;
    else
      an=i;
  }
}

int posInT11Ring (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int o = p.GetpFDeg();
  int op = set[length].GetpFDeg();

  if ((op < o)
  || ((op == o) && (pLtCmpOrdSgnDiffP(set[length].p,p.p))))
    return length+1;

  int i;
  int an = 0;
  int en= length;

  loop
  {
    if (an >= en-1)
    {
      op= set[an].GetpFDeg();
      if ((op > o)
      || (( op == o) && (pLtCmpOrdSgnEqP(set[an].p,p.p))))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg();
    if (( op > o)
    || (( op == o) && (pLtCmpOrdSgnEqP(set[i].p,p.p))))
      en=i;
    else
      an=i;
  }
}

/*2
* looks up the position of p in T
* set[0] is the smallest with respect to the ordering-procedure
* totaldegree,pComp
*/
int posInT110 (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;
  p.GetpLength();

  int o = p.GetpFDeg();
  int op = set[length].GetpFDeg();
  int cmp_int=currRing->OrdSgn;

  if (( op < o)
  || (( op == o) && (set[length].length<p.length))
  || (( op == o) && (set[length].length == p.length)
     && (pLmCmp(set[length].p,p.p) != cmp_int)))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg();
      if (( op > o)
      || (( op == o) && (set[an].length > p.length))
      || (( op == o) && (set[an].length == p.length)
         && (pLmCmp(set[an].p,p.p) == cmp_int)))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg();
    if (( op > o)
    || (( op == o) && (set[i].length > p.length))
    || (( op == o) && (set[i].length == p.length)
       && (pLmCmp(set[i].p,p.p) == cmp_int)))
      en=i;
    else
      an=i;
  }
}

int posInT110Ring (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;
  p.GetpLength();

  int o = p.GetpFDeg();
  int op = set[length].GetpFDeg();

  if (( op < o)
  || (( op == o) && (set[length].length<p.length))
  || (( op == o) && (set[length].length == p.length)
     && (pLtCmpOrdSgnDiffP(set[length].p,p.p))))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg();
      if (( op > o)
      || (( op == o) && (set[an].length > p.length))
      || (( op == o) && (set[an].length == p.length)
         && (pLtCmpOrdSgnEqP(set[an].p,p.p))))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg();
    if (( op > o)
    || (( op == o) && (set[i].length > p.length))
    || (( op == o) && (set[i].length == p.length)
       && (pLtCmpOrdSgnEqP(set[i].p,p.p))))
      en=i;
    else
      an=i;
  }
}

/*2
* looks up the position of p in set
* set[0] is the smallest with respect to the ordering-procedure
* pFDeg
*/
int posInT13 (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int o = p.GetpFDeg();

  if (set[length].GetpFDeg() <= o)
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      if (set[an].GetpFDeg() > o)
        return an;
      return en;
    }
    i=(an+en) / 2;
    if (set[i].GetpFDeg() > o)
      en=i;
    else
      an=i;
  }
}

// determines the position based on: 1.) Ecart 2.) pLength
int posInT_EcartpLength(const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;
  int ol = p.GetpLength();
  int op=p.ecart;
  int oo=set[length].ecart;

  if ((oo < op) || ((oo==op) && (set[length].length <= ol)))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      int oo=set[an].ecart;
      if((oo > op)
         || ((oo==op) && (set[an].pLength > ol)))
        return an;
      return en;
    }
    i=(an+en) / 2;
    int oo=set[i].ecart;
    if ((oo > op)
        || ((oo == op) && (set[i].pLength > ol)))
      en=i;
    else
      an=i;
  }
}

/*2
* looks up the position of p in set
* set[0] is the smallest with respect to the ordering-procedure
* maximaldegree, pComp
*/
int posInT15 (const BlockArray<TObject> &set,const int length,LObject &p)
/*{
 *int j=0;
 * int o;
 *
 * o = p.GetpFDeg()+p.ecart;
 * loop
 * {
 *   if ((set[j].GetpFDeg()+set[j].ecart > o)
 *   || ((set[j].GetpFDeg()+set[j].ecart == o)
 *     && (pLmCmp(set[j].p,p.p) == currRing->OrdSgn)))
 *   {
 *     return j;
 *   }
 *   j++;
 *   if (j > length) return j;
 * }
 *}
 */
{
  if (length==-1) return 0;

  int o = p.GetpFDeg() + p.ecart;
  int op = set[length].GetpFDeg()+set[length].ecart;
  int cmp_int=currRing->OrdSgn;

  if ((op < o)
  || ((op == o)
     && (pLmCmp(set[length].p,p.p) != cmp_int)))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg()+set[an].ecart;
      if (( op > o)
      || (( op  == o) && (pLmCmp(set[an].p,p.p) == cmp_int)))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg()+set[i].ecart;
    if (( op > o)
    || (( op == o) && (pLmCmp(set[i].p,p.p) == cmp_int)))
      en=i;
    else
      an=i;
  }
}

int posInT15Ring (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int o = p.GetpFDeg() + p.ecart;
  int op = set[length].GetpFDeg()+set[length].ecart;

  if ((op < o)
  || ((op == o)
     && (pLtCmpOrdSgnDiffP(set[length].p,p.p))))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg()+set[an].ecart;
      if (( op > o)
      || (( op  == o) && (pLtCmpOrdSgnEqP(set[an].p,p.p))))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg()+set[i].ecart;
    if (( op > o)
    || (( op == o) && (pLtCmpOrdSgnEqP(set[i].p,p.p))))
      en=i;
    else
      an=i;
  }
}

/*2
* looks up the position of p in set
* set[0] is the smallest with respect to the ordering-procedure
* pFDeg+ecart, ecart, pComp
*/
int posInT17 (const BlockArray<TObject> &set,const int length,LObject &p)
/*
*{
* int j=0;
* int  o;
*
*  o = p.GetpFDeg()+p.ecart;
*  loop
*  {
*    if ((pFDeg(set[j].p)+set[j].ecart > o)
*    || (((pFDeg(set[j].p)+set[j].ecart == o)
*      && (set[j].ecart < p.ecart)))
*    || ((pFDeg(set[j].p)+set[j].ecart == o)
*      && (set[j].ecart==p.ecart)
*      && (pLmCmp(set[j].p,p.p)==currRing->OrdSgn)))
*      return j;
*    j++;
*    if (j > length) return j;
*  }
* }
*/
{
  if (length==-1) return 0;

  int o = p.GetpFDeg() + p.ecart;
  int op = set[length].GetpFDeg()+set[length].ecart;
  int cmp_int=currRing->OrdSgn;

  if ((op < o)
  || (( op == o) && (set[length].ecart > p.ecart))
  || (( op == o) && (set[length].ecart==p.ecart)
     && (pLmCmp(set[length].p,p.p) != cmp_int)))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg()+set[an].ecart;
      if (( op > o)
      || (( op == o) && (set[an].ecart < p.ecart))
      || (( op  == o) && (set[an].ecart==p.ecart)
         && (pLmCmp(set[an].p,p.p) == cmp_int)))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg()+set[i].ecart;
    if ((op > o)
    || (( op == o) && (set[i].ecart < p.ecart))
    || (( op == o) && (set[i].ecart == p.ecart)
       && (pLmCmp(set[i].p,p.p) == cmp_int)))
      en=i;
    else
      an=i;
  }
}

int posInT17Ring (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int o = p.GetpFDeg() + p.ecart;
  int op = set[length].GetpFDeg()+set[length].ecart;

  if ((op < o)
  || (( op == o) && (set[length].ecart > p.ecart))
  || (( op == o) && (set[length].ecart==p.ecart)
     && (pLtCmpOrdSgnDiffP(set[length].p,p.p))))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg()+set[an].ecart;
      if (( op > o)
      || (( op == o) && (set[an].ecart < p.ecart))
      || (( op  == o) && (set[an].ecart==p.ecart)
         && (pLtCmpOrdSgnEqP(set[an].p,p.p))))
        return an;
      return en;
    }
    i=(an+en) / 2;
    op = set[i].GetpFDeg()+set[i].ecart;
    if ((op > o)
    || (( op == o) && (set[i].ecart < p.ecart))
    || (( op == o) && (set[i].ecart == p.ecart)
       && (pLtCmpOrdSgnEqP(set[i].p,p.p))))
      en=i;
    else
      an=i;
  }
}

/*2
* looks up the position of p in set
* set[0] is the smallest with respect to the ordering-procedure
* pGetComp, pFDeg+ecart, ecart, pComp
*/
int posInT17_c (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int cc = (-1+2*currRing->order[0]==ringorder_c);
  /* cc==1 for (c,..), cc==-1 for (C,..) */
  int o = p.GetpFDeg() + p.ecart;
  int c = pGetComp(p.p)*cc;
  int cmp_int=currRing->OrdSgn;

  if (pGetComp(set[length].p)*cc < c)
    return length+1;
  if (pGetComp(set[length].p)*cc == c)
  {
    int op = set[length].GetpFDeg()+set[length].ecart;
    if ((op < o)
    || ((op == o) && (set[length].ecart > p.ecart))
    || ((op == o) && (set[length].ecart==p.ecart)
       && (pLmCmp(set[length].p,p.p) != cmp_int)))
      return length+1;
  }

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      if (pGetComp(set[an].p)*cc < c)
        return en;
      if (pGetComp(set[an].p)*cc == c)
      {
        int op = set[an].GetpFDeg()+set[an].ecart;
        if ((op > o)
        || ((op == o) && (set[an].ecart < p.ecart))
        || ((op == o) && (set[an].ecart==p.ecart)
           && (pLmCmp(set[an].p,p.p) == cmp_int)))
          return an;
      }
      return en;
    }
    i=(an+en) / 2;
    if (pGetComp(set[i].p)*cc > c)
      en=i;
    else if (pGetComp(set[i].p)*cc == c)
    {
      int op = set[i].GetpFDeg()+set[i].ecart;
      if ((op > o)
      || ((op == o) && (set[i].ecart < p.ecart))
      || ((op == o) && (set[i].ecart == p.ecart)
         && (pLmCmp(set[i].p,p.p) == cmp_int)))
        en=i;
      else
        an=i;
    }
    else
      an=i;
  }
}

int posInT17_cRing (const BlockArray<TObject> &set,const int length,LObject &p)
{
  if (length==-1) return 0;

  int cc = (-1+2*currRing->order[0]==ringorder_c);
  /* cc==1 for (c,..), cc==-1 for (C,..) */
  int o = p.GetpFDeg() + p.ecart;
  int c = pGetComp(p.p)*cc;

  if (pGetComp(set[length].p)*cc < c)
    return length+1;
  if (pGetComp(set[length].p)*cc == c)
  {
    int op = set[length].GetpFDeg()+set[length].ecart;
    if ((op < o)
    || ((op == o) && (set[length].ecart > p.ecart))
    || ((op == o) && (set[length].ecart==p.ecart)
       && (pLtCmpOrdSgnDiffP(set[length].p,p.p))))
      return length+1;
  }

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      if (pGetComp(set[an].p)*cc < c)
        return en;
      if (pGetComp(set[an].p)*cc == c)
      {
        int op = set[an].GetpFDeg()+set[an].ecart;
        if ((op > o)
        || ((op == o) && (set[an].ecart < p.ecart))
        || ((op == o) && (set[an].ecart==p.ecart)
           && (pLtCmpOrdSgnEqP(set[an].p,p.p))))
          return an;
      }
      return en;
    }
    i=(an+en) / 2;
    if (pGetComp(set[i].p)*cc > c)
      en=i;
    else if (pGetComp(set[i].p)*cc == c)
    {
      int op = set[i].GetpFDeg()+set[i].ecart;
      if ((op > o)
      || ((op == o) && (set[i].ecart < p.ecart))
      || ((op == o) && (set[i].ecart == p.ecart)
         && (pLtCmpOrdSgnEqP(set[i].p,p.p))))
        en=i;
      else
        an=i;
    }
    else
      an=i;
  }
}

/*2
* looks up the position of p in set
* set[0] is the smallest with respect to
* ecart, pFDeg, length
*/
int posInT19 (const BlockArray<TObject> &set,const int length,LObject &p)
{
  p.GetpLength();
  if (length==-1) return 0;

  int o = p.ecart;
  int op=p.GetpFDeg();

  if (set[length].ecart < o)
    return length+1;
  if (set[length].ecart == o)
  {
    int oo=set[length].GetpFDeg();
    if ((oo < op) || ((oo==op) && (set[length].length < p.length)))
      return length+1;
  }

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      if (set[an].ecart > o)
        return an;
      if (set[an].ecart == o)
      {
        int oo=set[an].GetpFDeg();
        if((oo > op)
        || ((oo==op) && (set[an].length > p.length)))
          return an;
      }
      return en;
    }
    i=(an+en) / 2;
    if (set[i].ecart > o)
      en=i;
    else if (set[i].ecart == o)
    {
      int oo=set[i].GetpFDeg();
      if ((oo > op)
      || ((oo == op) && (set[i].length > p.length)))
        en=i;
      else
       an=i;
    }
    else
      an=i;
  }
}

/* Ordering procedure:
 *    - total degree (reversed)
 *    - p1 == NULL before p1 != NULL
 *    - leading monomial
 */

int compareLSpecial (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg();
  auto dr = rhs.GetpFDeg();
  if (dl < dr) return -1;
  if (dl > dr) return 1;

  if ((lhs.p1 == NULL) && (rhs.p1 != NULL)) return 1;
  if ((lhs.p1 != NULL) && (rhs.p1 == NULL)) return -1;

  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - leading monomial
 */

int compareL0 (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - leading monomial (ring version)
 */

int compareL0Ring (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  return (pLtCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - signature
 */

int compareLSig (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  return (pLtCmp(lhs.sig,rhs.sig) * currRing->OrdSgn);
}

/* UNUSED Ordering procedure:
 *    - signature
 *    - total degree
 *    - leading monomial
 */

int compareLSigRing (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  assume(currRing->OrdSgn == 1 && rField_is_Ring(currRing));
  auto cmp = pLtCmp(lhs.sig,rhs.sig);
  if (cmp != 1) return cmp;
  if (lhs.FDeg > rhs.FDeg) return -1;
  if (lhs.FDeg < rhs.FDeg) return 1;
  return pLtCmp(lhs.p, rhs.p);
}

// for sba, sorting syzygies
int posInSyz (const kStrategy strat, poly sig)
{
  if (strat->syzl==0) return 0;
  int cmp_int=currRing->OrdSgn;
  if (pLtCmp(strat->syz[strat->syzl-1],sig) != cmp_int)
    return strat->syzl;
  int i;
  int an = 0;
  int en= strat->syzl-1;
  loop
  {
    if (an >= en-1)
    {
      if (pLtCmp(strat->syz[an],sig) != cmp_int) return en;
      return an;
    }
    i=(an+en) / 2;
    if (pLtCmp(strat->syz[i],sig) != cmp_int) an=i;
    else                                      en=i;
    /*aend. fuer lazy == in !=- machen */
  }
}

/* Ordering procedure:
 *    - all elements equal
 *
 * Only used in F5C
 */

int compareLF5C (const LObject &, const LObject &, const kStrategy)
{
  return 0;
}

/* Ordering procedure:
 *    - total degree (reversed)
 *    - leading monomial
 */

int compareL11 (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg();
  auto dr = rhs.GetpFDeg();
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - total degree (reversed)
 *    - leading monomial (ring version)
 */

int compareL11Ring (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg();
  auto dr = rhs.GetpFDeg();
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  return (pLtCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - total degree (reversed)
 *    - coefficient (larger first)
 */

int compareL11Ringls (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  if (lhs.FDeg < rhs.FDeg) return -1;
  if (lhs.FDeg > rhs.FDeg) return 1;

  number lcl = pGetCoeff(lhs.p);
  number lcr = pGetCoeff(rhs.p);

  // Ensure positive coefficients
  if (!nGreaterZero(lcl))
    lcl = nInpNeg(nCopy(lcl));
  if (!nGreaterZero(lcr))
    lcr = nInpNeg(nCopy(lcr));

  if (nGreater(lcl, lcr)) return 1;
  if (nGreater(lcr, lcl)) return -1;
  return 0;
}

/*2 Position for rings L: Here I am
* looks up the position of polynomial p in set
* e is the ecart of p
* set[length] is the smallest element in set with respect
* to the ordering-procedure totaldegree,pComp
*/
inline int getIndexRng(long coeff)
{
  if (coeff == 0) return -1;
  long tmp = coeff;
  int ind = 0;
  while (tmp % 2 == 0)
  {
    tmp = tmp / 2;
    ind++;
  }
  return ind;
}

/*{
  if (length < 0) return 0;

  int o = p->GetpFDeg();
  int op = set[length].GetpFDeg();

  int inde = getIndexRng((unsigned long) pGetCoeff(set[length].p));
  int indp = getIndexRng((unsigned long) pGetCoeff(p->p));
  int inda;
  int indi;

  if ((inda > indp) || ((inda == inde) && ((op > o) || ((op == o) && (pLmCmp(set[length].p,p->p) != -currRing->OrdSgn)))))
    return length + 1;
  int i;
  int an = 0;
  inda = getIndexRng((unsigned long) pGetCoeff(set[an].p));
  int en = length;
  loop
  {
    if (an >= en-1)
    {
      op = set[an].GetpFDeg();
      if ((indp > inda) || ((indp == inda) && ((op > o) || ((op == o) && (pLmCmp(set[an].p,p->p) != -currRing->OrdSgn)))))
        return en;
      return an;
    }
    i = (an + en) / 2;
    indi = getIndexRng((unsigned long) pGetCoeff(set[i].p));
    op = set[i].GetpFDeg();
    if ((indi > indp) || ((indi == indp) && ((op > o) || ((op == o) && (pLmCmp(set[i].p,p->p) != -currRing->OrdSgn)))))
    // if ((op > o) || ((op == o) && (pLmCmp(set[i].p,p->p) != -currRing->OrdSgn)))
    {
      an = i;
      inda = getIndexRng((unsigned long) pGetCoeff(set[an].p));
    }
    else
      en = i;
  }
} */

/* Ordering procedure:
 *    - total degree (reversed)
 *    - length (reversed)
 *    - leading monomial
 */

int compareL110 (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg();
  auto dr = rhs.GetpFDeg();
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  if (lhs.length < rhs.length) return -1;
  if (lhs.length > rhs.length) return 1;
  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - total degree (reversed)
 *    - length (reversed)
 *    - leading monomial (ring version)
 */

int compareL110Ring (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg();
  auto dr = rhs.GetpFDeg();
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  if (lhs.length < rhs.length) return -1;
  if (lhs.length > rhs.length) return 1;
  return (pLtCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - total degree (reversed)
 */

int compareL13 (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg();
  auto dr = rhs.GetpFDeg();
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  return 0;
}

/* Ordering procedure:
 *    - maximal degree (reversed)
 *    - leading monomial
 */

int compareL15 (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg() + lhs.ecart;
  auto dr = rhs.GetpFDeg() + rhs.ecart;
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - maximal degree (reversed)
 *    - leading monomial (ring version)
 */

int compareL15Ring (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg() + lhs.ecart;
  auto dr = rhs.GetpFDeg() + rhs.ecart;
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  return (pLtCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - maximal degree (reversed)
 *    - ecart (reversed)
 *    - leading monomial
 */

int compareL17 (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg() + lhs.ecart;
  auto dr = rhs.GetpFDeg() + rhs.ecart;
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  if (lhs.ecart < rhs.ecart) return -1;
  if (lhs.ecart > rhs.ecart) return 1;
  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - maximal degree (reversed)
 *    - ecart (reversed)
 *    - leading monomial (ring version)
 */

int compareL17Ring (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  auto dl = lhs.GetpFDeg() + lhs.ecart;
  auto dr = rhs.GetpFDeg() + rhs.ecart;
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  if (lhs.ecart < rhs.ecart) return -1;
  if (lhs.ecart > rhs.ecart) return 1;
  return (pLtCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - component (reversed)
 *    - maximal degree (reversed)
 *    - ecart (reversed)
 *    - leading monomial
 */

int compareL17_c (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  int cc = (-1+2*currRing->order[0]==ringorder_c);
  long cl = pGetComp(lhs.p)*cc;
  long cr = pGetComp(rhs.p)*cc;
  if (cl < cr) return -1;
  if (cl > cr) return 1;

  auto dl = lhs.GetpFDeg() + lhs.ecart;
  auto dr = rhs.GetpFDeg() + rhs.ecart;
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  if (lhs.ecart < rhs.ecart) return -1;
  if (lhs.ecart > rhs.ecart) return 1;
  return (pLmCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/* Ordering procedure:
 *    - component (reversed)
 *    - maximal degree (reversed)
 *    - ecart (reversed)
 *    - leading monomial (ring version)
 */

int compareL17_cRing (const LObject &lhs, const LObject &rhs, const kStrategy)
{
  int cc = (-1+2*currRing->order[0]==ringorder_c);
  long cl = pGetComp(lhs.p)*cc;
  long cr = pGetComp(rhs.p)*cc;
  if (cl < cr) return -1;
  if (cl > cr) return 1;

  auto dl = lhs.GetpFDeg() + lhs.ecart;
  auto dr = rhs.GetpFDeg() + rhs.ecart;
  if (dl < dr) return -1;
  if (dl > dr) return 1;
  if (lhs.ecart < rhs.ecart) return -1;
  if (lhs.ecart > rhs.ecart) return 1;
  return (pLtCmp(lhs.p,rhs.p) * currRing->OrdSgn);
}

/*
 * SYZYGY CRITERION for signature-based standard basis algorithms
 */
BOOLEAN syzCriterion(poly sig, unsigned long not_sevSig, kStrategy strat)
{
//#if 1
#ifdef DEBUGF5
  PrintS("syzygy criterion checks:  ");
  pWrite(sig);
#endif
  for (int k=0; k<strat->syzl; k++)
  {
    //printf("-%d",k);
//#if 1
#ifdef DEBUGF5
    Print("checking with: %d / %d --  \n",k,strat->syzl);
    pWrite(pHead(strat->syz[k]));
#endif
    if (p_LmShortDivisibleBy(strat->syz[k], strat->sevSyz[k], sig, not_sevSig, currRing)
    && (!rField_is_Ring(currRing) ||
    (n_DivBy(pGetCoeff(sig), pGetCoeff(strat->syz[k]),currRing->cf) && pLtCmp(sig,strat->syz[k]) == 1)))
    {
//#if 1
#ifdef DEBUGF5
      PrintS("DELETE!\n");
#endif
      strat->nrsyzcrit++;
      //printf("- T -\n\n");
      return TRUE;
    }
  }
  //printf("- F -\n\n");
  return FALSE;
}

/*
 * SYZYGY CRITERION for signature-based standard basis algorithms
 */
BOOLEAN syzCriterionInc(poly sig, unsigned long not_sevSig, kStrategy strat)
{
//#if 1
  if(sig == NULL)
    return FALSE;
#ifdef DEBUGF5
  PrintS("--- syzygy criterion checks:  ");
  pWrite(sig);
#endif
  int comp = (int)__p_GetComp(sig, currRing);
  int min, max;
  if (comp<=1)
    return FALSE;
  else
  {
    min = strat->syzIdx[comp-2];
    //printf("SYZIDX %d/%d\n",strat->syzIdx[comp-2],comp-2);
    //printf("SYZIDX %d/%d\n",strat->syzIdx[comp-1],comp-1);
    //printf("SYZIDX %d/%d\n",strat->syzIdx[comp],comp);
    if (comp == strat->currIdx)
    {
      max = strat->syzl;
    }
    else
    {
      max = strat->syzIdx[comp-1];
    }
    for (int k=min; k<max; k++)
    {
#ifdef F5DEBUG
      Print("COMP %d/%d - MIN %d - MAX %d - SYZL %ld\n",comp,strat->currIdx,min,max,strat->syzl);
      Print("checking with: %d --  ",k);
      pWrite(pHead(strat->syz[k]));
#endif
      if (p_LmShortDivisibleBy(strat->syz[k], strat->sevSyz[k], sig, not_sevSig, currRing)
      && (!rField_is_Ring(currRing) ||
      (n_DivBy(pGetCoeff(sig), pGetCoeff(strat->syz[k]),currRing->cf) && pLtCmp(sig,strat->syz[k]) == 1)))
      {
        strat->nrsyzcrit++;
        return TRUE;
      }
    }
    return FALSE;
  }
}

/*
 * REWRITTEN CRITERION for signature-based standard basis algorithms
 */
BOOLEAN faugereRewCriterion(poly sig, unsigned long not_sevSig, poly /*lm*/, kStrategy strat, sBasisSet::const_iterator start)
{
  //printf("Faugere Rewritten Criterion\n");
  if(rField_is_Ring(currRing))
    return FALSE;
//#if 1
#ifdef DEBUGF5
  PrintS("rewritten criterion checks:  ");
  pWrite(sig);
#endif
  // Scan from start forward (previously backward from S.size()-1 to
  // start; direction is irrelevant for correctness since we return
  // TRUE on any match).
  for (auto sit_k = start; sit_k != strat->S.cend(); ++sit_k)
  {
//#if 1
#ifdef DEBUGF5
    PrintS("checking with:  ");
    pWrite(sit_k->sig);
    pWrite(pHead(sit_k->p));
#endif
    if (p_LmShortDivisibleBy(sit_k->sig, sit_k->sevSig, sig, not_sevSig, currRing))
    {
//#if 1
#ifdef DEBUGF5
      PrintS("DELETE!\n");
#endif
      strat->nrrewcrit++;
      return TRUE;
    }
  }
#ifdef DEBUGF5
  PrintS("ALL ELEMENTS OF S\n----------------------------------------\n");
  for (auto sit_kk = strat->S.begin(); sit_kk != strat->S.end(); ++sit_kk)
  {
    pWrite(pHead(sit_kk->p));
  }
  PrintS("------------------------------\n");
#endif
  return FALSE;
}

/*
 * REWRITTEN CRITERION for signature-based standard basis algorithms
 ***************************************************************************
 * TODO:This should become the version of Arri/Perry resp. Bjarke/Stillman *
 ***************************************************************************
 */

// real implementation of arri's rewritten criterion, only called once in
// kstd2.cc, right before starting reduction
// IDEA:  Arri says that it is enough to consider 1 polynomial for each unique
//        signature appearing during the computations. Thus we first of all go
//        through strat->L and delete all other pairs of the same signature,
//        keeping only the one with least possible leading monomial. After this
//        we check if we really need to compute this critical pair at all: There
//        can be elements already in strat->S whose signatures divide the
//        signature of the critical pair in question and whose multiplied
//        leading monomials are smaller than the leading monomial of the
//        critical pair. In this situation we can discard the critical pair
//        completely.
BOOLEAN arriRewCriterion(poly /*sig*/, unsigned long /*not_sevSig*/, poly /*lm*/, kStrategy strat, sBasisSet::const_iterator start)
{
  if(rField_is_Ring(currRing))
    return FALSE;
  poly p1 = pOne();
  poly p2 = pOne();
  // Old loop was `ii > start` (strict), i.e. skipped the element at
  // `start` itself.  With the iterator API the "+1 collapse" means
  // callers now pass the iterator AT the first element to check, so
  // we use `!=` (inclusive).  rewCrit2's sole caller collapsed its +1;
  // rewCrit3 callers still pass S.cend() (no change in scan range).
  for (auto sit_ii = start; sit_ii != strat->S.cend(); ++sit_ii)
  {
    if (p_LmShortDivisibleBy(sit_ii->sig, sit_ii->sevSig, strat->P.sig, ~strat->P.sevSig, currRing))
    {
      p_ExpVectorSum(p1,strat->P.sig,sit_ii->p,currRing);
      p_ExpVectorSum(p2,sit_ii->sig,strat->P.p,currRing);
      if (!(pLmCmp(p1,p2) == 1))
      {
        pDelete(&p1);
        pDelete(&p2);
        return TRUE;
      }
    }
  }
  pDelete(&p1);
  pDelete(&p2);
  return FALSE;
}

BOOLEAN arriRewCriterionPre(poly sig, unsigned long not_sevSig, poly lm, kStrategy strat, sBasisSet::const_iterator /*start*/)
{
  //Over Rings, there are still some changes to do: considering coeffs
  if(rField_is_Ring(currRing))
    return FALSE;
  auto found = strat_B(strat).uend();
  for (auto it = strat_B(strat).ubegin(); it != strat_B(strat).uend(); ++it)
  {
    if (pLmEqual(it->sig,sig))
    {
      found = it;
      break;
    }
  }
  if (found != strat_B(strat).uend())
  {
    if (pLmCmp(lm,found->GetLmCurrRing()) == -1)
    {
      strat_B(strat).erase(found);
    }
    else
    {
      return TRUE;
    }
  }
  poly p1 = pOne();
  poly p2 = pOne();
  for (auto sit_ii = strat->S.begin(); sit_ii != strat->S.end(); ++sit_ii)
  {
    if (p_LmShortDivisibleBy(sit_ii->sig, sit_ii->sevSig, sig, not_sevSig, currRing))
    {
      p_ExpVectorSum(p1,sig,sit_ii->p,currRing);
      p_ExpVectorSum(p2,sit_ii->sig,lm,currRing);
      if (!(pLmCmp(p1,p2) == 1))
      {
        pDelete(&p1);
        pDelete(&p2);
        return TRUE;
      }
    }
  }
  pDelete(&p1);
  pDelete(&p2);
  return FALSE;
}

/***************************************************************
 *
 * Tail reductions
 *
 ***************************************************************/
TObject* kFindDivisibleByInS_T(kStrategy strat, sBasisSet::const_iterator end, LObject* L, TObject *T, long ecart)
{
  const unsigned long not_sev = ~L->sev;
  poly p;
  ring r;
  L->GetLm(p, r);

  assume(~not_sev == p_GetShortExpVector(p, r));

  // Walk S as iterator; `end` is the exclusive upper-bound iterator.
  auto at_end = [&](sBasisSet::iterator it) {
    return it == strat->S.end() || it.index() >= end.index();
  };

  if (r == currRing)
  {
    auto sit = strat->S.begin();
    if(!rField_is_Ring(r))
    {
      for (; !at_end(sit); ++sit)
      {
  #if defined(PDEBUG) || defined(PDIV_DEBUG)
        if (sit->p != NULL && p_LmShortDivisibleBy(sit->p, sit->sev, p, not_sev, r) &&
            (ecart== LONG_MAX || ecart >= sit->ecart))
  #else
        if (!(sit->sev & not_sev) &&
            (ecart== LONG_MAX || ecart >= sit->ecart) &&
            p_LmDivisibleBy(sit->p, p, r))
  #endif
        {
          break;
        }
      }
    }
    else
    {
      for (; !at_end(sit); ++sit)
      {
  #if defined(PDEBUG) || defined(PDIV_DEBUG)
        if (sit->p != NULL
        && p_LmShortDivisibleBy(sit->p, sit->sev, p, not_sev, r)
        && (ecart== LONG_MAX || ecart >= sit->ecart)
        && n_DivBy(pGetCoeff(p), pGetCoeff(sit->p), r->cf))
  #else
        if (!(sit->sev & not_sev)
        && (ecart== LONG_MAX || ecart >= sit->ecart)
        && p_LmDivisibleBy(sit->p, p, r)
        && n_DivBy(pGetCoeff(p), pGetCoeff(sit->p), r->cf))
  #endif
        {
          break; // found
        }
      }
    }
    if (at_end(sit)) return NULL;
    // if called from NF, T objects do not exist:
    if (strat->T.empty() || sit->s_2_r == -1)
    {
      T->Set(sit->p, r, strat->tailRing);
      assume(T->GetpLength()==pLength(T->p != __null ? T->p : T->t_p));
      return T;
    }
    else
    {
      return strat->S.S_2_T(sit, strat);
    }
  }
  else
  {
    TObject* t;
    auto sit = strat->S.begin();
    if(!rField_is_Ring(r))
    {
      for (; !at_end(sit); ++sit)
      {
        assume(sit->s_2_r != -1);
  #if defined(PDEBUG) || defined(PDIV_DEBUG)
        t = strat->S.S_2_T(sit, strat);
        assume(t != NULL && t->t_p != NULL && t->tailRing == r);
        if (p_LmShortDivisibleBy(t->t_p, sit->sev, p, not_sev, r)
        && (ecart== LONG_MAX || ecart >= sit->ecart))
        {
          kt_debug_tag("kFindInS_T:6591", (void*)t, t->i_r, (int)pLength(t->t_p));
          t->pLength=pLength(t->t_p);
          return t;
        }
  #else
        if (! (sit->sev & not_sev)
        && (ecart== LONG_MAX || ecart >= sit->ecart))
        {
          if (sit->s_2_r == -1)
          {
            T->Set(sit->p, r, strat->tailRing);
            if (p_LmDivisibleBy(T->t_p != NULL ? T->t_p : T->p, p, r))
            {
              T->pLength=pLength(T->t_p != NULL ? T->t_p : T->p);
              return T;
            }
            continue;
          }
          t = strat->S.S_2_T(sit, strat);
          assume(t != NULL && t->t_p != NULL && t->tailRing == r && t->p == sit->p);
          if (p_LmDivisibleBy(t->t_p, p, r))
          {
            kt_debug_tag("kFindInS_T:6612", (void*)t, t->i_r, (int)pLength(t->t_p));
            t->pLength=pLength(t->t_p);
            return t;
          }
        }
  #endif
      }
      return NULL;
    }
    else
    {
      for (; !at_end(sit); ++sit)
      {
        assume(sit->s_2_r != -1);
  #if defined(PDEBUG) || defined(PDIV_DEBUG)
        t = strat->S.S_2_T(sit, strat);
        assume(t != NULL && t->t_p != NULL && t->tailRing == r);
        if (p_LmShortDivisibleBy(t->t_p, sit->sev, p, not_sev, r)
        && (ecart== LONG_MAX || ecart >= sit->ecart)
        && n_DivBy(pGetCoeff(p), pGetCoeff(t->t_p), r->cf))
        {
          kt_debug_tag("kFindInS_T:6632", (void*)t, t->i_r, (int)pLength(t->t_p));
          t->pLength=pLength(t->t_p);
          return t;
        }
  #else
        if (! (sit->sev & not_sev)
        && (ecart== LONG_MAX || ecart >= sit->ecart))
        {
          if (sit->s_2_r == -1)
          {
            T->Set(sit->p, r, strat->tailRing);
            if (p_LmDivisibleBy(T->t_p != NULL ? T->t_p : T->p, p, r))
            {
              T->pLength=pLength(T->t_p != NULL ? T->t_p : T->p);
              return T;
            }
            continue;
          }
          t = strat->S.S_2_T(sit, strat);
          assume(t != NULL && t->t_p != NULL && t->tailRing == r && t->p == sit->p);
          if (p_LmDivisibleBy(t->t_p, p, r)
          && n_DivBy(pGetCoeff(p), pGetCoeff(t->t_p), r->cf))
          {
            kt_debug_tag("kFindInS_T:6654", (void*)t, t->i_r, (int)pLength(t->t_p));
            t->pLength=pLength(t->t_p);
            return t;
          }
        }
  #endif
      }
      return NULL;
    }
  }
}

poly redtail (LObject* L, sBasisSet::const_iterator end, kStrategy strat)
{
  poly h, hn;
  strat->redTailChange=FALSE;

  L->GetP();
  poly p = L->p;
  if (strat->noTailReduction || pNext(p) == NULL)
    return p;

  LObject Ln(strat->tailRing);
  TObject* With;
  // placeholder in case strat->T.empty()
  TObject  With_s(strat->tailRing);
  h = p;
  hn = pNext(h);
  long op = strat->tailRing->pFDeg(hn, strat->tailRing);
  long e;
  int l;
  BOOLEAN save_HE=strat->kAllAxis;
  strat->kAllAxis |=
    ((Kstd1_deg>0) && (op<=Kstd1_deg)) || TEST_OPT_INFREDTAIL;

  while(hn != NULL)
  {
    op = strat->tailRing->pFDeg(hn, strat->tailRing);
    if ((Kstd1_deg>0)&&(op>Kstd1_deg)) goto all_done;
    e = strat->tailRing->pLDeg(hn, &l, strat->tailRing) - op;
    loop
    {
      Ln.Set(hn, strat->tailRing);
      Ln.sev = p_GetShortExpVector(hn, strat->tailRing);
      if (strat->kAllAxis)
        With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s);
      else
        With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s, e);
      if (With == NULL) break;
      // Suspect #1 from the static-audit starting list: writes
      // With->pLength=0 on a pointer returned by
      // kFindDivisibleByInS_T, which could be a T entry shared with
      // other threads.  Tag with With->i_r so we can correlate with
      // a later ksCreateSpoly STALE_l2 event on the same index.
      kt_debug_tag("redtail:With->pLength=0", (void*)With,
                   With->i_r, (int)With->pLength);
      With->length=0;
      With->pLength=0;
      strat->redTailChange=TRUE;
      if (ksReducePolyTail(L, With, h, strat->kNoetherTail()))
      {
        // reducing the tail would violate the exp bound
        if (kStratChangeTailRing(strat, L))
        {
          strat->kAllAxis = save_HE;
          return redtail(L, end, strat);
        }
        else
          return NULL;
      }
      hn = pNext(h);
      if (hn == NULL) goto all_done;
      op = strat->tailRing->pFDeg(hn, strat->tailRing);
      if ((Kstd1_deg>0)&&(op>Kstd1_deg)) goto all_done;
      e = strat->tailRing->pLDeg(hn, &l, strat->tailRing) - op;
    }
    h = hn;
    hn = pNext(h);
  }

  all_done:
  if (strat->redTailChange)
  {
    L->pLength = 0;
  }
  strat->kAllAxis = save_HE;
  return p;
}

poly redtail (poly p, sBasisSet::const_iterator end, kStrategy strat)
{
  LObject L(p, currRing);
  return redtail(&L, end, strat);
}

// `end` is the exclusive upper-bound iterator: reduction considers
// S-elements in [strat->S.begin(), end).
poly redtailBba (LObject* L, sBasisSet::const_iterator end, kStrategy strat, BOOLEAN withT, BOOLEAN normalize)
{
  strat->redTailChange=FALSE;
  if (strat->noTailReduction) return L->GetLmCurrRing();
  poly h, p;
  p = h = L->GetLmTailRing();
  if ((h==NULL) || (pNext(h)==NULL))
    return L->GetLmCurrRing();

  TObject* With;
  // placeholder in case strat->T.empty()
  TObject  With_s(strat->tailRing);

  LObject Ln(pNext(h), strat->tailRing);
  Ln.GetpLength();

  pNext(h) = NULL;
  if (L->p != NULL)
  {
    pNext(L->p) = NULL;
    if (L->t_p != NULL) pNext(L->t_p) = NULL;
  }
  L->pLength = 1;

  Ln.PrepareRed(strat->use_buckets);

  int cnt=REDTAIL_CANONICALIZE;
  while(!Ln.IsNull())
  {
    loop
    {
      if (TEST_OPT_IDLIFT)
      {
        if (Ln.p!=NULL)
        {
          if ((int)__p_GetComp(Ln.p,currRing)> strat->syzComp) break;
        }
        else
        {
          if ((int)__p_GetComp(Ln.t_p,strat->tailRing)> strat->syzComp) break;
        }
      }
      Ln.SetShortExpVector();
      if (withT)
      {
        int j;
        j = kFindDivisibleByInT(strat, &Ln);
        if (j < 0) break;
        // Task 511 t-iterator-migrate-hot-path: route the T[j] access
        // through the skipping iterator.  kFindDivisibleByInT returns
        // the physical index of a slot it already observed as published
        // (via iterator walk or tobject_published_load gate), so
        // iterator_at_T lands directly on T[j] without skipping.
        auto t_it = iterator_at_T(strat->T, j);
        With = &(*t_it);
        assume(With->GetpLength()==pLength(With->p != __null ? With->p : With->t_p));
      }
      else
      {
        With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s);
        if (With == NULL) break;
        assume(With->GetpLength()==pLength(With->p != __null ? With->p : With->t_p));
      }
      cnt--;
      if (cnt==0)
      {
        cnt=REDTAIL_CANONICALIZE;
        /*poly tmp=*/Ln.CanonicalizeP();
        if (normalize)
        {
          Ln.Normalize();
          //pNormalize(tmp);
          //if (TEST_OPT_PROT) { PrintS("n"); mflush(); }
        }
      }
      if (normalize && (!TEST_OPT_INTSTRATEGY) && (!nIsOne(pGetCoeff(With->p))))
      {
        With->pNorm();
      }
      strat->redTailChange=TRUE;
      if (ksReducePolyTail(L, With, &Ln))
      {
        // reducing the tail would violate the exp bound
        //  set a flag and hope for a retry (in bba)
        strat->completeReduce_retry=TRUE;
        if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
        do
        {
          pNext(h) = Ln.LmExtractAndIter();
          pIter(h);
          L->pLength++;
        } while (!Ln.IsNull());
        goto all_done;
      }
      if (Ln.IsNull()) goto all_done;
      if (! withT) With_s.Init(currRing);
    }
    pNext(h) = Ln.LmExtractAndIter();
    pIter(h);
    pNormalize(h);
    L->pLength++;
  }

  all_done:
  Ln.Delete();
  if (L->p != NULL) pNext(L->p) = pNext(p);

  if (strat->redTailChange)
  {
    L->length = 0;
    L->pLength = 0;
  }

  //if (TEST_OPT_PROT) { PrintS("N"); mflush(); }
  //L->Normalize(); // HANNES: should have a test
  kTest_L(L,strat);
  return L->GetLmCurrRing();
}

poly redtailBbaBound (LObject* L, sBasisSet::const_iterator end, kStrategy strat, int bound, BOOLEAN withT, BOOLEAN normalize)
{
  strat->redTailChange=FALSE;
  if (strat->noTailReduction) return L->GetLmCurrRing();
  poly h, p;
  p = h = L->GetLmTailRing();
  if ((h==NULL) || (pNext(h)==NULL))
    return L->GetLmCurrRing();

  TObject* With;
  // placeholder in case strat->T.empty()
  TObject  With_s(strat->tailRing);

  LObject Ln(pNext(h), strat->tailRing);
  Ln.pLength = L->GetpLength() - 1;

  pNext(h) = NULL;
  if (L->p != NULL) pNext(L->p) = NULL;
  L->pLength = 1;

  Ln.PrepareRed(strat->use_buckets);

  int cnt=REDTAIL_CANONICALIZE;
  while(!Ln.IsNull())
  {
    loop
    {
      if (TEST_OPT_IDLIFT)
      {
        if (Ln.p!=NULL)
        {
          if ((int)__p_GetComp(Ln.p,currRing)> strat->syzComp) break;
        }
        else
        {
          if ((int)__p_GetComp(Ln.t_p,strat->tailRing)> strat->syzComp) break;
        }
      }
      Ln.SetShortExpVector();
      if (withT)
      {
        int j;
        j = kFindDivisibleByInT(strat, &Ln);
        if (j < 0) break;
        // Task 511 t-iterator-migrate-hot-path: route T[j] access through
        // the skipping iterator (j is already known published by
        // kFindDivisibleByInT).
        auto t_it = iterator_at_T(strat->T, j);
        With = &(*t_it);
      }
      else
      {
        With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s);
        if (With == NULL) break;
      }
      cnt--;
      if (cnt==0)
      {
        cnt=REDTAIL_CANONICALIZE;
        /*poly tmp=*/Ln.CanonicalizeP();
        if (normalize)
        {
          Ln.Normalize();
          //pNormalize(tmp);
          //if (TEST_OPT_PROT) { PrintS("n"); mflush(); }
        }
      }
      if (normalize && (!TEST_OPT_INTSTRATEGY) && (!nIsOne(pGetCoeff(With->p))))
      {
        With->pNorm();
      }
      strat->redTailChange=TRUE;
      if (ksReducePolyTail(L, With, &Ln))
      {
        // reducing the tail would violate the exp bound
        //  set a flag and hope for a retry (in bba)
        strat->completeReduce_retry=TRUE;
        if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
        do
        {
          pNext(h) = Ln.LmExtractAndIter();
          pIter(h);
          L->pLength++;
        } while (!Ln.IsNull());
        goto all_done;
      }
      if(!Ln.IsNull())
      {
        Ln.GetP();
        Ln.p = pJet(Ln.p,bound);
      }
      if (Ln.IsNull())
      {
        goto all_done;
      }
      if (! withT) With_s.Init(currRing);
    }
    pNext(h) = Ln.LmExtractAndIter();
    pIter(h);
    pNormalize(h);
    L->pLength++;
  }

  all_done:
  Ln.Delete();
  if (L->p != NULL) pNext(L->p) = pNext(p);

  if (strat->redTailChange)
  {
    L->length = 0;
    L->pLength = 0;
  }

  //if (TEST_OPT_PROT) { PrintS("N"); mflush(); }
  //L->Normalize(); // HANNES: should have a test
  kTest_L(L,strat);
  return L->GetLmCurrRing();
}

void redtailBbaAlsoLC_Z (LObject* L, kStrategy strat )
// normalize=FALSE, withT=FALSE, coeff=Z
{
  strat->redTailChange=FALSE;

  poly h, p;
  p = h = L->GetLmTailRing();
  if ((h==NULL) || (pNext(h)==NULL))
    return;

  TObject* With;
  LObject Ln(pNext(h), strat->tailRing);
  Ln.GetpLength();

  pNext(h) = NULL;
  if (L->p != NULL)
  {
    pNext(L->p) = NULL;
    if (L->t_p != NULL) pNext(L->t_p) = NULL;
  }
  L->pLength = 1;

  Ln.PrepareRed(strat->use_buckets);

  int cnt=REDTAIL_CANONICALIZE;

  while(!Ln.IsNull())
  {
    loop
    {
      if (TEST_OPT_IDLIFT)
      {
        if (Ln.p!=NULL)
        {
          if ((int)__p_GetComp(Ln.p,currRing)> strat->syzComp) break;
        }
        else
        {
          if ((int)__p_GetComp(Ln.t_p,strat->tailRing)> strat->syzComp) break;
        }
      }
      Ln.SetShortExpVector();
      int j;
      j = kFindDivisibleByInT(strat, &Ln);
      if (j < 0)
      {
        j = kFindDivisibleByInT_Z(strat, &Ln);
        if (j < 0)
        {
          break;
        }
        else
        {
          /* reduction not cancelling a tail term, but reducing its coefficient */
          With = &(strat->T[j]);
          assume(With->GetpLength()==pLength(With->p != __null ? With->p : With->t_p));
          cnt--;
          if (cnt==0)
          {
            cnt=REDTAIL_CANONICALIZE;
            /*poly tmp=*/Ln.CanonicalizeP();
          }
          strat->redTailChange=TRUE;
          /* reduction cancelling a tail term */
          if (ksReducePolyTailLC_Z(L, With, &Ln))
          {
            // reducing the tail would violate the exp bound
            //  set a flag and hope for a retry (in bba)
            strat->completeReduce_retry=TRUE;
            if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
            do
            {
              pNext(h) = Ln.LmExtractAndIter();
              pIter(h);
              L->pLength++;
            } while (!Ln.IsNull());
            goto all_done;
          }
          /* we have to break since we did not cancel the term, but only decreased
           * its coefficient. */
          break;
        }
      } else {
        With = &(strat->T[j]);
        assume(With->GetpLength()==pLength(With->p != __null ? With->p : With->t_p));
        cnt--;
        if (cnt==0)
        {
          cnt=REDTAIL_CANONICALIZE;
          /*poly tmp=*/Ln.CanonicalizeP();
        }
        strat->redTailChange=TRUE;
        /* reduction cancelling a tail term */
        if (ksReducePolyTail_Z(L, With, &Ln))
        {
          // reducing the tail would violate the exp bound
          //  set a flag and hope for a retry (in bba)
          strat->completeReduce_retry=TRUE;
          if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
          do
          {
            pNext(h) = Ln.LmExtractAndIter();
            pIter(h);
            L->pLength++;
          } while (!Ln.IsNull());
          goto all_done;
        }
      }
      if (Ln.IsNull()) goto all_done;
    }
    pNext(h) = Ln.LmExtractAndIter();
    pIter(h);
    L->pLength++;
  }

  all_done:
  Ln.Delete();
  if (L->p != NULL) pNext(L->p) = pNext(p);

  if (strat->redTailChange)
  {
    L->length = 0;
    L->pLength = 0;
  }

  kTest_L(L, strat);
  return;
}

poly redtailBba_Z (LObject* L, sBasisSet::const_iterator end, kStrategy strat )
// normalize=FALSE, withT=FALSE, coeff=Z
{
  strat->redTailChange=FALSE;
  if (strat->noTailReduction) return L->GetLmCurrRing();
  poly h, p;
  p = h = L->GetLmTailRing();
  if ((h==NULL) || (pNext(h)==NULL))
    return L->GetLmCurrRing();

  TObject* With;
  // placeholder in case strat->T.empty()
  TObject  With_s(strat->tailRing);

  LObject Ln(pNext(h), strat->tailRing);
  Ln.pLength = L->GetpLength() - 1;

  pNext(h) = NULL;
  if (L->p != NULL) pNext(L->p) = NULL;
  L->pLength = 1;

  Ln.PrepareRed(strat->use_buckets);

  int cnt=REDTAIL_CANONICALIZE;
  while(!Ln.IsNull())
  {
    loop
    {
      Ln.SetShortExpVector();
      With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s);
      if (With == NULL) break;
      cnt--;
      if (cnt==0)
      {
        cnt=REDTAIL_CANONICALIZE;
        /*poly tmp=*/Ln.CanonicalizeP();
      }
      // we are in Z, do not call pNorm
      strat->redTailChange=TRUE;
      // test divisibility of coefs:
      Ln.GetLmCurrRing();
      With->GetLmCurrRing();

      if (ksReducePolyTail_Z(L, With, &Ln))
      {
        // reducing the tail would violate the exp bound
        //  set a flag and hope for a retry (in bba)
        strat->completeReduce_retry=TRUE;
        if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
        do
        {
          pNext(h) = Ln.LmExtractAndIter();
          pIter(h);
          L->pLength++;
        } while (!Ln.IsNull());
        goto all_done;
      }
      if (Ln.IsNull()) goto all_done;
      With_s.Init(currRing);
    }
    pNext(h) = Ln.LmExtractAndIter();
    pIter(h);
    pNormalize(h);
    L->pLength++;
  }

  all_done:
  Ln.Delete();
  if (L->p != NULL) pNext(L->p) = pNext(p);

  if (strat->redTailChange)
  {
    L->length = 0;
  }

  //if (TEST_OPT_PROT) { PrintS("N"); mflush(); }
  //L->Normalize(); // HANNES: should have a test
  kTest_L(L,strat);
  return L->GetLmCurrRing();
}

poly redtailBba_NF (poly p, kStrategy strat )
{
  strat->redTailChange=FALSE;
  if (strat->noTailReduction) return p;
  if ((p==NULL) || (pNext(p)==NULL))
    return p;

  poly h=p;
  p=pNext(p);
  pNext(h)=NULL;
  while(p!=NULL)
  {
    p=redNF(p,1,strat);
    if (p!=NULL)
    {
      poly hh=p;
      p=pNext(p);
      pNext(hh)=NULL;
      h=p_Add_q(h,hh,currRing);
    }
  }
  return h;
}

poly redtailBba_Ring (LObject* L, sBasisSet::const_iterator end, kStrategy strat )
// normalize=FALSE, withT=FALSE, coeff=Ring
{
  strat->redTailChange=FALSE;
  if (strat->noTailReduction) return L->GetLmCurrRing();
  poly h, p;
  p = h = L->GetLmTailRing();
  if ((h==NULL) || (pNext(h)==NULL))
    return L->GetLmCurrRing();

  TObject* With;
  // placeholder in case strat->T.empty()
  TObject  With_s(strat->tailRing);

  LObject Ln(pNext(h), strat->tailRing);
  Ln.pLength = L->GetpLength() - 1;

  pNext(h) = NULL;
  if (L->p != NULL) pNext(L->p) = NULL;
  L->pLength = 1;

  Ln.PrepareRed(strat->use_buckets);

  int cnt=REDTAIL_CANONICALIZE;
  while(!Ln.IsNull())
  {
    loop
    {
      Ln.SetShortExpVector();
      With_s.Init(currRing);
      With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s);
      if (With == NULL) break;
      cnt--;
      if (cnt==0)
      {
        cnt=REDTAIL_CANONICALIZE;
        /*poly tmp=*/Ln.CanonicalizeP();
      }
      // we are in a ring, do not call pNorm
      // test divisibility of coefs:
      poly p_Ln=Ln.GetLmCurrRing();
      poly p_With=With->GetLmCurrRing();
      if (n_DivBy(pGetCoeff(p_Ln),pGetCoeff(p_With), currRing->cf))
      {
        strat->redTailChange=TRUE;

        if (ksReducePolyTail_Z(L, With, &Ln))
        {
          // reducing the tail would violate the exp bound
          //  set a flag and hope for a retry (in bba)
          strat->completeReduce_retry=TRUE;
          if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
          do
          {
            pNext(h) = Ln.LmExtractAndIter();
            pIter(h);
            L->pLength++;
          } while (!Ln.IsNull());
          goto all_done;
        }
      }
      else break; /*proceed to next monomial*/
      if (Ln.IsNull()) goto all_done;
    }
    pNext(h) = Ln.LmExtractAndIter();
    pIter(h);
    pNormalize(h);
    L->pLength++;
  }

  all_done:
  Ln.Delete();
  if (L->p != NULL) pNext(L->p) = pNext(p);

  if (strat->redTailChange)
  {
    L->length = 0;
  }

  //if (TEST_OPT_PROT) { PrintS("N"); mflush(); }
  //L->Normalize(); // HANNES: should have a test
  kTest_L(L,strat);
  return L->GetLmCurrRing();
}

/*2
*checks the change degree and write progress report
*/
void message (int i,int* olddeg,LSet::size_type* reduc,kStrategy strat, int red_result)
{
  if (i != *olddeg)
  {
    Print("%d",i);
    *olddeg = i;
  }
  if (TEST_OPT_OLDSTD)
  {
    if (strat->L.size() != *reduc+1)
    {
      if (strat->L.size() != *reduc)
        Print("(%zu)",strat->L.size());
      else
        PrintS("-");
      *reduc = strat->L.size()-1;
    }
    else
      PrintS(".");
    mflush();
  }
  else
  {
    if (red_result == 0)
      PrintS("-");
    else if (red_result < 0)
      PrintS(".");
    if ((red_result > 0) || ((strat->L.size() % 100)==0))
    {
      if (strat->L.size() != *reduc+1 && strat->L.size() > 1)
      {
        Print("(%zu)",strat->L.size());
        *reduc = strat->L.size()-1;
      }
    }
  }
}

/*2
*statistics
*/
void messageStat (int hilbcount,kStrategy strat)
{
  //PrintS("\nUsage/Allocation of temporary storage:\n");
  //Print("%d/%d polynomials in standard base\n",srmax,IDELEMS(Shdl));
  //Print("%d/%d polynomials in set L (for lazy alg.)",strat->L.size(),strat->L.capacity());
  Print("product criterion:%d chain criterion:%d\n",strat->cp,strat->c3);
  if (hilbcount!=0) Print("hilbert series criterion:%d\n",hilbcount);
  #ifdef HAVE_SHIFTBBA
  /* in usual case strat->cv is 0, it gets changed only in shift routines */
  if (strat->cv!=0) Print("shift V criterion:%d\n",strat->cv);
  #endif
}

void messageStatSBA (int hilbcount,kStrategy strat)
{
  //PrintS("\nUsage/Allocation of temporary storage:\n");
  //Print("%d/%d polynomials in standard base\n",srmax,IDELEMS(Shdl));
  //Print("%d/%d polynomials in set L (for lazy alg.)",strat->L.size(),strat->L.capacity());
  Print("syz criterion:%d rew criterion:%d\n",strat->nrsyzcrit,strat->nrrewcrit);
  //Print("product criterion:%d chain criterion:%d\n",strat->cp,strat->c3);
  if (hilbcount!=0) Print("hilbert series criterion:%d\n",hilbcount);
  #ifdef HAVE_SHIFTBBA
  /* in usual case strat->cv is 0, it gets changed only in shift routines */
  if (strat->cv!=0) Print("shift V criterion:%d\n",strat->cv);
  #endif
}

#ifdef KDEBUG
/*2
*debugging output: all internal sets, if changed
*for testing purpose only/has to be changed for later use
*/
void messageSets (kStrategy strat)
{
  int i;
  if (strat->news)
  {
    PrintS("set S");
    for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
    {
      Print("\n  %d:",sit.index());
      p_wrp(sit->p, currRing, strat->tailRing);
      if (strat->hasFromQ && sit->fromQ)
        Print(" (from Q)");
    }
    strat->news = FALSE;
  }
  if (strat->newt)
  {
    PrintS("\nset T");
    for (i=0; i < strat->T.size(); i++)
    {
      Print("\n  %d:",i);
      strat->T[i].wrp();
      if (strat->T[i].length==0) strat->T[i].length=pLength(strat->T[i].p);
      Print(" o:%ld e:%d l:%d",
        strat->T[i].pFDeg(),strat->T[i].ecart,strat->T[i].length);
    }
    strat->newt = FALSE;
  }
  PrintS("\nset L");
  i=strat->L.size()-1;
  for (auto& Lp: strat->L) {
    Print("\n%d:",i);
    p_wrp(Lp.p1, currRing, strat->tailRing);
    PrintS("  ");
    p_wrp(Lp.p2, currRing, strat->tailRing);
    PrintS(" lcm: ");p_wrp(Lp.lcm, currRing);
    PrintS("\n  p : ");
    Lp.wrp();
    Print("  o:%ld e:%d l:%d",
          Lp.pFDeg(),Lp.ecart,Lp.length);
    i--;
  }
  PrintLn();
}

#endif


/*2
*construct the set s from F
*/
void initS (ideal F, ideal Q, kStrategy strat)
{
  int   i,pos;

  if (Q!=NULL) i=((IDELEMS(F)+IDELEMS(Q)+(setmaxTinc-1))/setmaxTinc)*setmaxTinc;
  else         i=((IDELEMS(F)+(setmaxTinc-1))/setmaxTinc)*setmaxTinc;
  if (i<setmaxTinc) i=setmaxT;
  strat->S.ensure_capacity(i);
  strat->hasFromQ=FALSE;
  strat->Srank=F->rank;
  /*- put polys into S -*/
  if (Q!=NULL)
  {
    strat->hasFromQ=TRUE;
    for (i=0; i<IDELEMS(Q); i++)
    {
      if (Q->m[i]!=NULL)
      {
        LObject h;
        h.p = pCopy(Q->m[i]);
        if (TEST_OPT_INTSTRATEGY)
        {
          h.pCleardenom(); // also does remove Content
        }
        else
        {
          h.pNorm();
        }
        if (rHasLocalOrMixedOrdering(currRing))
        {
          deleteHC(&h, strat);
        }
        if (h.p!=NULL)
        {
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          h.fromQ = 1;
          strat->enterS(h, strat, -1, strat->S.end());
        }
      }
    }
  }
  for (i=0; i<IDELEMS(F); i++)
  {
    if (F->m[i]!=NULL)
    {
      LObject h;
      h.p = pCopy(F->m[i]);
      if (rHasLocalOrMixedOrdering(currRing))
      {
        cancelunit(&h);  /*- tries to cancel a unit -*/
        deleteHC(&h, strat);
      }
      if (h.p!=NULL)
      // do not rely on the input being a SB!
      {
        if (TEST_OPT_INTSTRATEGY)
        {
          h.pCleardenom(); // also does remove Content
        }
        else
        {
          h.pNorm();
        }
        strat->initEcart(&h);
        h.sev = pGetShortExpVector(h.p);
        strat->enterS(h, strat, -1, strat->S.end());
      }
    }
  }
  /*- test, if a unit is in F -*/
  if (((!strat->S.empty()))
       && n_IsUnit(pGetCoeff(strat->S.begin()->p),currRing->cf)
       && pIsConstant(strat->S.begin()->p))
  {
    while (strat->S.size() > 1)
    {
      auto last = strat->S.end(); --last;
      strat->S.erase(last);
    }
  }
}

void initSL (ideal F, ideal Q,kStrategy strat)
{
  int   i,pos;

  if (Q!=NULL)
  {
    i=((IDELEMS(Q)+(setmaxTinc-1))/setmaxTinc)*setmaxTinc;
    if (i<setmaxTinc) i=setmaxT;
  }
  else i=setmaxT;
  strat->S.ensure_capacity(i);
  strat->hasFromQ=FALSE;
  strat->Srank=F->rank;
  /*- put polys into S -*/
  if (Q!=NULL)
  {
    strat->hasFromQ=TRUE;
    for (i=0; i<IDELEMS(Q); i++)
    {
      if (Q->m[i]!=NULL)
      {
        LObject h;
        h.p = pCopy(Q->m[i]);
        if (rHasLocalOrMixedOrdering(currRing))
        {
          deleteHC(&h,strat);
          cancelunit(&h);  /*- tries to cancel a unit -*/
        }
        if (TEST_OPT_INTSTRATEGY)
        {
          h.pCleardenom(); // also does remove Content
        }
        else
        {
          h.pNorm();
        }
        if (h.p!=NULL)
        {
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          h.fromQ = 1;
          strat->enterS(h, strat, -1, strat->S.end());
        }
        if(errorreported) return;
      }
    }
  }
  for (i=0; i<IDELEMS(F); i++)
  {
    if (F->m[i]!=NULL)
    {
      LObject h;
      h.p = pCopy(F->m[i]);
      if (h.p!=NULL)
      {
        if (rHasLocalOrMixedOrdering(currRing))
        {
          cancelunit(&h);  /*- tries to cancel a unit -*/
          deleteHC(&h, strat);
        }
        if (h.p!=NULL)
        {
          if (TEST_OPT_INTSTRATEGY)
          {
            h.pCleardenom(); // also does remove Content
          }
          else
          {
            h.pNorm();
          }
          if(errorreported) return;
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          strat->L.push(h);
        }
      }
    }
  }
  /*- test, if a unit is in F -*/

  if (! strat->L.empty()
       && n_IsUnit(pGetCoeff(strat->L.top().p), currRing->cf)
       && pIsConstant(strat->L.top().p))
  {
    // pop the top LObject, but don't erase it because we're going to put it back
    auto unit = strat->L.top();
    strat->L.pop();
    while (! strat->L.empty()) strat->L.pop_and_erase();
    strat->L.push(unit);
  }
}

void initSLSba (ideal F, ideal Q,kStrategy strat)
{
  int   i,pos;
  if (Q!=NULL)
  {
    i=((IDELEMS(Q)+(setmaxTinc-1))/setmaxTinc)*setmaxTinc;
    if (i<setmaxTinc) i=setmaxT;
  }
  else i=setmaxT;
  strat->S.ensure_capacity(i);
  strat->hasFromQ  =   FALSE;
  strat->Srank   =   F->rank;
  if (strat->sbaOrder != 1)
  {
    strat->syz    = (poly *)omAlloc0(i*sizeof(poly));
    strat->sevSyz = initsevS(i);
    strat->syzmax = i;
    strat->syzl   = 0;
  }
  /*- put polys into S -*/
  if (Q!=NULL)
  {
    strat->hasFromQ=TRUE;
    for (i=0; i<IDELEMS(Q); i++)
    {
      if (Q->m[i]!=NULL)
      {
        LObject h;
        h.p = pCopy(Q->m[i]);
        if (rHasLocalOrMixedOrdering(currRing))
        {
          deleteHC(&h,strat);
        }
        if (TEST_OPT_INTSTRATEGY)
        {
          h.pCleardenom(); // also does remove Content
        }
        else
        {
          h.pNorm();
        }
        if (h.p!=NULL)
        {
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          h.fromQ = 1;
          strat->enterS(h, strat, -1, strat->S.end());
        }
      }
    }
  }
  for (i=0; i<IDELEMS(F); i++)
  {
    if (F->m[i]!=NULL)
    {
      LObject h;
      h.p = pCopy(F->m[i]);
      h.sig = pOne();
      //h.sig = pInit();
      //p_SetCoeff(h.sig,nInit(1),currRing);
      p_SetComp(h.sig,i+1,currRing);
      // if we are working with the Schreyer order we generate it
      // by multiplying the initial signatures with the leading monomial
      // of the corresponding initial polynomials generating the ideal
      // => we can keep the underlying monomial order and get a Schreyer
      //    order without any bigger overhead
      if (strat->sbaOrder == 0 || strat->sbaOrder == 3)
      {
        p_ExpVectorAdd (h.sig,F->m[i],currRing);
      }
      h.sevSig = pGetShortExpVector(h.sig);
#ifdef DEBUGF5
      pWrite(h.p);
      pWrite(h.sig);
#endif
      if (h.p!=NULL)
      {
        if (rHasLocalOrMixedOrdering(currRing))
        {
          cancelunit(&h);  /*- tries to cancel a unit -*/
          deleteHC(&h, strat);
        }
        if (h.p!=NULL)
        {
          if (TEST_OPT_INTSTRATEGY)
          {
            h.pCleardenom(); // also does remove Content
          }
          else
          {
            h.pNorm();
          }
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          // Initial L entries have no prior rewCrit2 scan history, so
          // resume-scan starts from the beginning of S.  Without this
          // explicit init, h.checked would be the null iterator from
          // sLObject::Init()'s memset, segfaulting in rewCrit2 at the
          // first iteration of the main sba() loop.
          h.checked = strat->S.cbegin();
          strat->L.push(h);
        }
      }
      /*
      if (strat->sbaOrder != 1)
      {
        for(j=0;j<i;j++)
        {
          strat->syz[ctr] = pCopy(F->m[j]);
          p_SetCompP(strat->syz[ctr],i+1,currRing);
          // add LM(F->m[i]) to the signature to get a Schreyer order
          // without changing the underlying polynomial ring at all
          p_ExpVectorAdd (strat->syz[ctr],F->m[i],currRing);
          // since p_Add_q() destroys all input
          // data we need to recreate help
          // each time
          poly help = pCopy(F->m[i]);
          p_SetCompP(help,j+1,currRing);
          pWrite(strat->syz[ctr]);
          pWrite(help);
          printf("%d\n",pLmCmp(strat->syz[ctr],help));
          strat->syz[ctr] = p_Add_q(strat->syz[ctr],help,currRing);
          printf("%d. SYZ  ",ctr);
          pWrite(strat->syz[ctr]);
          strat->sevSyz[ctr] = p_GetShortExpVector(strat->syz[ctr],currRing);
          ctr++;
        }
        strat->syzl = ps;
      }
      */
    }
  }
  /*- test, if a unit is in F -*/

  if (! strat->L.empty()
       && n_IsUnit(pGetCoeff(strat->L.top().p), currRing->cf)
       && pIsConstant(strat->L.top().p))
  {
    // pop the top LObject, but don't erase it because we're going to put it back
    auto unit = strat->L.top();
    strat->L.pop();
    while (! strat->L.empty()) strat->L.pop_and_erase();
    strat->L.push(unit);
  }
}

void initSyzRules (kStrategy strat)
{
  if( strat->S.begin()->p )
  {
    // Note: size() >= 2 is not explicitly checked here, matching the
    // pre-migration behaviour (iterator_at(1) was also unchecked).
    auto second = strat->S.begin(); ++second;
    if( second->p && !rField_is_Ring(currRing))
    {
      omFreeSize(strat->syzIdx,(strat->syzidxmax)*sizeof(int));
      omFreeSize(strat->sevSyz,(strat->syzmax)*sizeof(unsigned long));
      omFreeSize(strat->syz,(strat->syzmax)*sizeof(poly));
    }
    int j, diff, comp, comp_old, ps=0, ctr=0;
    /************************************************************
     * computing the length of the syzygy array needed
     ***********************************************************/
    {
      // Trailing-pair iterator walk: prev=i-1, cur=i. Count is cur.index().
      auto prev = strat->S.begin();
      auto cur  = prev; ++cur;
      for (; cur != strat->S.end(); ++prev, ++cur)
      {
        if (pGetComp(prev->sig) != pGetComp(cur->sig))
        {
          ps += cur.index();
        }
      }
    }
    ps += strat->S.size();
    //comp              = pGetComp (strat->P.sig);
    comp              = strat->currIdx;
    strat->syzIdx     = initec(comp);
    strat->sevSyz     = initsevS(ps);
    strat->syz        = (poly *)omAlloc(ps*sizeof(poly));
    strat->syzmax     = ps;
    strat->syzl       = 0;
    strat->syzidxmax  = comp;
#if defined(DEBUGF5) || defined(DEBUGF51)
    PrintS("------------- GENERATING SYZ RULES NEW ---------------\n");
#endif
    j = 0;
    /************************************************************
     * generating the leading terms of the principal syzygies
     ***********************************************************/
    {
      auto prev = strat->S.begin();
      auto sit_i = prev; ++sit_i;
      for (; sit_i != strat->S.end(); ++prev, ++sit_i)
      {
      /**********************************************************
       * principal syzygies start with component index 2
       * the array syzIdx starts with index 0
       * => the rules for a signature with component comp start
       *    at strat->syz[strat->syzIdx[comp-2]] !
       *********************************************************/
      if (pGetComp(prev->sig) != pGetComp(sit_i->sig))
      {
        comp      = pGetComp(sit_i->sig);
        comp_old  = pGetComp(prev->sig);
        diff      = comp - comp_old - 1;
        // diff should be zero, but sometimes also the initial generating
        // elements of the input ideal reduce to zero. then there is an
        // index-gap between the signatures. for these in-between signatures we
        // can safely set syzIdx[j] = 0 as no such element will be ever computed
        // in the following.
        // doing this, we keep the relation "j = comp - 2" alive, which makes
        // jumps way easier when checking criteria
        while (diff>0)
        {
          strat->syzIdx[j]  = 0;
          diff--;
          j++;
        }
        strat->syzIdx[j]  = ctr;
        j++;
        LObject Q;
        int pos;
        for (auto sit_k = strat->S.begin(); sit_k != sit_i; ++sit_k)
        {
          Q.sig          = pOne();
          if(rField_is_Ring(currRing))
            p_SetCoeff(Q.sig,nCopy(p_GetCoeff(sit_k->p,currRing)),currRing);
          p_ExpVectorCopy(Q.sig,sit_k->p,currRing);
          p_SetCompP (Q.sig, comp, currRing);
          poly q          = p_One(currRing);
          if(rField_is_Ring(currRing))
            p_SetCoeff(q,nCopy(p_GetCoeff(sit_i->p,currRing)),currRing);
          p_ExpVectorCopy(q,sit_i->p,currRing);
          q               = p_Neg (q, currRing);
          p_SetCompP (q, __p_GetComp(sit_k->sig, currRing), currRing);
          Q.sig = p_Add_q (Q.sig, q, currRing);
          Q.sevSig  = p_GetShortExpVector(Q.sig,currRing);
          pos = posInSyz(strat, Q.sig);
          enterSyz(Q, strat, pos);
          ctr++;
        }
      }
      }
    }
    /**************************************************************
    * add syzygies for upcoming first element of new iteration step
    **************************************************************/
    comp      = strat->currIdx;
    {
      auto sit_last = strat->S.end(); --sit_last;
      comp_old  = pGetComp(sit_last->sig);
    }
    diff      = comp - comp_old - 1;
    // diff should be zero, but sometimes also the initial generating
    // elements of the input ideal reduce to zero. then there is an
    // index-gap between the signatures. for these in-between signatures we
    // can safely set syzIdx[j] = 0 as no such element will be ever computed
    // in the following.
    // doing this, we keep the relation "j = comp - 2" alive, which makes
    // jumps way easier when checking criteria
    while (diff>0)
    {
      strat->syzIdx[j]  = 0;
      diff--;
      j++;
    }
    strat->syzIdx[j]  = ctr;
    LObject Q;
    int pos;
    for (auto sit_k = strat->S.begin(); sit_k != strat->S.end(); ++sit_k)
    {
      Q.sig          = pOne();
      if(rField_is_Ring(currRing))
        p_SetCoeff(Q.sig,nCopy(p_GetCoeff(sit_k->p,currRing)),currRing);
      p_ExpVectorCopy(Q.sig,sit_k->p,currRing);
      p_SetCompP (Q.sig, comp, currRing);
      poly q          = p_One(currRing);
      if(rField_is_Ring(currRing))
        p_SetCoeff(q,nCopy(p_GetCoeff(strat->L.top().p,currRing)),currRing);
      p_ExpVectorCopy(q,strat->L.top().p,currRing);
      q               = p_Neg (q, currRing);
      p_SetCompP (q, __p_GetComp(sit_k->sig, currRing), currRing);
      Q.sig = p_Add_q (Q.sig, q, currRing);
      Q.sevSig = p_GetShortExpVector(Q.sig,currRing);
      pos = posInSyz(strat, Q.sig);
      enterSyz(Q, strat, pos);
      ctr++;
    }
//#if 1
#ifdef DEBUGF5
    PrintS("Principal syzygies:\n");
    Print("syzl   %d\n",strat->syzl);
    Print("syzmax %d\n",strat->syzmax);
    Print("ps     %d\n",ps);
    PrintS("--------------------------------\n");
    for(i=0;i<=strat->syzl-1;i++)
    {
      Print("%d - ",i);
      pWrite(strat->syz[i]);
    }
    for(i=0;i<strat->currIdx;i++)
    {
      Print("%d - %d\n",i,strat->syzIdx[i]);
    }
    PrintS("--------------------------------\n");
#endif
  }
}

/*2
*construct the set s from F and {P}
*/
void initSSpecial (ideal F, ideal Q, ideal P,kStrategy strat)
{
  int   i,pos;

  if (Q!=NULL)
  {
    i=((IDELEMS(Q)+(setmaxTinc-1))/setmaxTinc)*setmaxTinc;
    if (i<setmaxTinc) i=setmaxT;
  }
  else i=setmaxT;
  i=((i+IDELEMS(F)+IDELEMS(P)+setmax-1)/setmax)*setmax;
  strat->S.ensure_capacity(i);
  strat->hasFromQ=FALSE;
  strat->Srank=F->rank;

  /*- put polys into S -*/
  if (Q!=NULL)
  {
    strat->hasFromQ=TRUE;
    for (i=0; i<IDELEMS(Q); i++)
    {
      if (Q->m[i]!=NULL)
      {
        LObject h;
        h.p = pCopy(Q->m[i]);
        //if (TEST_OPT_INTSTRATEGY)
        //{
        //  h.pCleardenom(); // also does remove Content
        //}
        //else
        //{
        //  h.pNorm();
        //}
        if (rHasLocalOrMixedOrdering(currRing))
        {
          deleteHC(&h,strat);
        }
        if (h.p!=NULL)
        {
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          h.fromQ = 1;
          strat->enterS(h, strat, strat->T.size(), strat->S.end());
          enterT(h, strat);
        }
      }
    }
  }
  /*- put polys into S -*/
  for (i=0; i<IDELEMS(F); i++)
  {
    if (F->m[i]!=NULL)
    {
      LObject h;
      h.p = pCopy(F->m[i]);
      if (rHasLocalOrMixedOrdering(currRing))
      {
        deleteHC(&h,strat);
      }
      else if (TEST_OPT_REDTAIL || TEST_OPT_REDSB)
      {
        h.p=redtailBba(h.p,strat->S.end(),strat);
      }
      if (h.p!=NULL)
      {
        strat->initEcart(&h);
        h.sev = pGetShortExpVector(h.p);
        strat->enterS(h, strat, strat->T.size(), strat->S.end());
        enterT(h,strat);
      }
    }
  }
  for (i=0; i<IDELEMS(P); i++)
  {
    if (P->m[i]!=NULL)
    {
      LObject h;
      h.p=pCopy(P->m[i]);
      if (TEST_OPT_INTSTRATEGY)
      {
        h.pCleardenom();
      }
      else
      {
        h.pNorm();
      }
      if((!strat->S.empty()))
      {
        if (rHasGlobalOrdering(currRing))
        {
          h.p=redBba(h.p,strat->S.end(),strat);
          if ((h.p!=NULL)&&(TEST_OPT_REDTAIL || TEST_OPT_REDSB))
          {
            h.p=redtailBba(h.p,strat->S.end(),strat);
          }
        }
        else
        {
          h.p=redMora(h.p,strat->S.end(),strat);
        }
        if(h.p!=NULL)
        {
          strat->initEcart(&h);
          if (TEST_OPT_INTSTRATEGY)
          {
            h.pCleardenom();
          }
          else
          {
            h.is_normalized = 0;
            h.pNorm();
          }
          h.sev = pGetShortExpVector(h.p);
          h.SetpFDeg();
          auto pos_it = strat->S.find_pos(h.p,h.ecart);
          enterpairsSpecial(h.p,strat->S.size()-1,h.ecart,pos_it,strat,strat->T.size());
          strat->enterS(h, strat, strat->T.size(), strat->S.end());
          enterT(h,strat);
        }
      }
      else
      {
        h.sev = pGetShortExpVector(h.p);
        strat->initEcart(&h);
        strat->enterS(h, strat, strat->T.size(), strat->S.end());
        enterT(h,strat);
      }
    }
  }
}
/*2
*construct the set s from F and {P}
*/

void initSSpecialSba (ideal F, ideal Q, ideal P,kStrategy strat)
{
  int   i,pos;

  if (Q!=NULL)
  {
    i=((IDELEMS(Q)+(setmaxTinc-1))/setmaxTinc)*setmaxTinc;
    if (i<setmaxTinc) i=setmaxT;
  }
  else i=setmaxT;
  i=((i+IDELEMS(F)+IDELEMS(P)+setmax-1)/setmax)*setmax;
  strat->S.ensure_capacity(i);
  strat->hasFromQ=FALSE;
  strat->Srank=F->rank;
  /*- put polys into S -*/
  if (Q!=NULL)
  {
    strat->hasFromQ=TRUE;
    for (i=0; i<IDELEMS(Q); i++)
    {
      if (Q->m[i]!=NULL)
      {
        LObject h;
        h.p = pCopy(Q->m[i]);
        //if (TEST_OPT_INTSTRATEGY)
        //{
        //  h.pCleardenom(); // also does remove Content
        //}
        //else
        //{
        //  h.pNorm();
        //}
        if (rHasLocalOrMixedOrdering(currRing))
        {
          deleteHC(&h,strat);
        }
        if (h.p!=NULL)
        {
          strat->initEcart(&h);
          h.sev = pGetShortExpVector(h.p);
          h.fromQ = 1;
          strat->enterS(h, strat, strat->T.size(), strat->S.end());
          enterT(h, strat);
        }
      }
    }
  }
  /*- put polys into S -*/
  for (i=0; i<IDELEMS(F); i++)
  {
    if (F->m[i]!=NULL)
    {
      LObject h;
      h.p = pCopy(F->m[i]);
      if (rHasLocalOrMixedOrdering(currRing))
      {
        deleteHC(&h,strat);
      }
      else if (TEST_OPT_REDTAIL || TEST_OPT_REDSB)
      {
        h.p=redtailBba(h.p,strat->S.end(),strat);
      }
      if (h.p!=NULL)
      {
        strat->initEcart(&h);
        h.sev = pGetShortExpVector(h.p);
        strat->enterS(h, strat, strat->T.size(), strat->S.end());
        enterT(h,strat);
      }
    }
  }
  for (i=0; i<IDELEMS(P); i++)
  {
    if (P->m[i]!=NULL)
    {
      LObject h;
      h.p=pCopy(P->m[i]);
      if (TEST_OPT_INTSTRATEGY)
      {
        h.pCleardenom();
      }
      else
      {
        h.pNorm();
      }
      if((!strat->S.empty()))
      {
        if (rHasGlobalOrdering(currRing))
        {
          h.p=redBba(h.p,strat->S.end(),strat);
          if ((h.p!=NULL)&&(TEST_OPT_REDTAIL || TEST_OPT_REDSB))
          {
            h.p=redtailBba(h.p,strat->S.end(),strat);
          }
        }
        else
        {
          h.p=redMora(h.p,strat->S.end(),strat);
        }
        if(h.p!=NULL)
        {
          strat->initEcart(&h);
          if (TEST_OPT_INTSTRATEGY)
          {
            h.pCleardenom();
          }
          else
          {
            h.is_normalized = 0;
            h.pNorm();
          }
          h.sev = pGetShortExpVector(h.p);
          h.SetpFDeg();
          auto pos_it = strat->S.find_pos(h.p,h.ecart);
          enterpairsSpecial(h.p,strat->S.size()-1,h.ecart,pos_it,strat,strat->T.size());
          strat->enterS(h, strat, strat->T.size(), strat->S.end());
          enterT(h,strat);
        }
      }
      else
      {
        h.sev = pGetShortExpVector(h.p);
        strat->initEcart(&h);
        strat->enterS(h, strat, strat->T.size(), strat->S.end());
        enterT(h,strat);
      }
    }
  }
}

/*2
* reduces h using the set S
* procedure used in cancelunit1
*/
static poly redBba1 (poly h,int maxIndex,kStrategy strat)
{
  unsigned long not_sev = ~ pGetShortExpVector(h);

  for (auto sit = strat->S.begin(); sit != strat->S.end() && sit.index() <= maxIndex; ++sit)
  {
    if (pLmShortDivisibleBy(sit->p, sit->sev, h, not_sev))
       return ksOldSpolyRedNew(sit->p, h, strat->kNoetherTail());
  }
  return h;
}

/*2
*tests if p.p=monomial*unit and cancels the unit
*/
void cancelunit1 (LObject* p,int *suc, int index,kStrategy strat )
{
  int k;
  poly r,h,h1,q;

  if (!pIsVector((*p).p) && ((*p).ecart != 0))
  {
    // Leading coef have to be a unit: no
    // example 2x+4x2 should be simplified to 2x*(1+2x)
    // and 2 is not a unit in Z
    //if ( !(n_IsUnit(pGetCoeff((*p).p), currRing->cf)) ) return;
    k = 0;
    h1 = r = pCopy((*p).p);
    h =pNext(r);
    loop
    {
      if (h==NULL)
      {
        pDelete(&r);
        pDelete(&(pNext((*p).p)));
        (*p).ecart = 0;
        (*p).length = 1;
        (*p).pLength = 1;
        (*suc)=0;
        return;
      }
      if (!pDivisibleBy(r,h))
      {
        q=redBba1(h,index ,strat);
        if (q != h)
        {
          k++;
          pDelete(&h);
          pNext(h1) = h = q;
        }
        else
        {
          pDelete(&r);
          return;
        }
      }
      else
      {
        h1 = h;
        pIter(h);
      }
      if (k > 10)
      {
        pDelete(&r);
        return;
      }
    }
  }
}

#if 0
/*2
* reduces h using the elements from Q in the set S
* procedure used in updateS
* must not be used for elements of Q or elements of an ideal !
*/
static poly redQ (poly h, int /*j*/, kStrategy strat)
{
  unsigned long not_sev = ~ pGetShortExpVector(h);
  auto sit = strat->S.begin();
  while (sit != strat->S.end() && pGetComp(sit->p) != 0) ++sit;
  auto start = sit;
  while (sit != strat->S.end())
  {
    if (pLmShortDivisibleBy(sit->p, sit->sev, h, not_sev))
    {
      h = ksOldSpolyRed(sit->p, h, strat->kNoetherTail());
      if (h==NULL) return NULL;
      sit = start;
      not_sev = ~ pGetShortExpVector(h);
    }
    else ++sit;
  }
  return h;
}
#endif

/*2
* reduces h using the set S
* procedure used in updateS
*/
static poly redBba (poly h,sBasisSet::const_iterator end,kStrategy strat)
{
  unsigned long not_sev = ~ pGetShortExpVector(h);

  sBasisSet::const_iterator sit = strat->S.cbegin();
  while (sit != end)
  {
    if (pLmShortDivisibleBy(sit->p, sit->sev, h, not_sev))
    {
      h = ksOldSpolyRed(sit->p, h, strat->kNoetherTail());
      if (h==NULL) return NULL;
      sit = strat->S.cbegin();
      not_sev = ~ pGetShortExpVector(h);
    }
    else ++sit;
  }
  return h;
}

/*2
* reduces h using the set S
*e is the ecart of h
*procedure used in updateS
*/
static poly redMora (poly h,sBasisSet::const_iterator end,kStrategy strat)
{
  int  e,l;
  unsigned long not_sev = ~ pGetShortExpVector(h);

  if (end != strat->S.cbegin())
  {
    e = currRing->pLDeg(h,&l,currRing)-p_FDeg(h,currRing);
    sBasisSet::const_iterator sit = strat->S.cbegin();
    do
    {
      if (pLmShortDivisibleBy(sit->p, sit->sev, h, not_sev)
      && ((e >= sit->ecart) || (strat->kNoether!=NULL)))
      {
#ifdef KDEBUG
        if (TEST_OPT_DEBUG)
        {
          PrintS("reduce ");wrp(h);Print(" with S[%d] (",sit.index());wrp(sit->p);
        }
#endif
        h = ksOldSpolyRed(sit->p, h, strat->kNoetherTail());
#ifdef KDEBUG
        if(TEST_OPT_DEBUG)
        {
          PrintS(")\nto "); wrp(h); PrintLn();
        }
#endif
        // pDelete(&h);
        if (h == NULL) return NULL;
        e = currRing->pLDeg(h,&l,currRing)-p_FDeg(h,currRing);
        sit = strat->S.cbegin();
        not_sev = ~ pGetShortExpVector(h);
      }
      else ++sit;
    }
    while (sit != end);
  }
  return h;
}

/*2
*updates S:
*the result is a set of polynomials which are in
*normalform with respect to S
*/
void updateS(BOOLEAN toT,kStrategy strat)
{
  LObject h;
  int suc=0;
  poly redSi=NULL;
  BOOLEAN change,any_change;
//  Print("nach initS: updateS start mit sl=%d\n",(strat->S.size()-1));
//  for (i=0; i<=(strat->S.size()-1); i++)
//  {
//    Print("s%d:",i);
//    if (strat->hasFromQ) Print("(Q:%d) ",strat->S[i].fromQ);
//    pWrite(strat->S[i].p);
//  }
//  Print("currRing->OrdSgn=%d\n", currRing->OrdSgn);
  any_change=FALSE;
  if (rHasGlobalOrdering(currRing))
  {
    while (suc != -1)
    {
      // Walk a live iterator starting at index suc+1. After reorder we
      // re-enter here from begin() since positions have shuffled.
      auto sit = strat->S.begin();
      for (int skip = suc+1; skip > 0 && sit != strat->S.end(); --skip) ++sit;
      while (sit != strat->S.end())
      {
        change=FALSE;
        if(rField_is_Ring(currRing))
            any_change = FALSE;
        if (((!strat->hasFromQ) || (sit->fromQ==0)) && (sit != strat->S.begin()))
        {
          redSi = pHead(sit->p);
          sit->p = redBba(sit->p,sit,strat);
          //if ((strat->ak!=0)&&(strat->S[i].p!=NULL))
          //  strat->S[i].p=redQ(strat->S[i].p,i+1,strat); /*reduce S[i] mod Q*/
          if (pCmp(redSi,sit->p)!=0)
          {
            change=TRUE;
            any_change=TRUE;
            #ifdef KDEBUG
            if (TEST_OPT_DEBUG)
            {
              PrintS("reduce:");
              wrp(redSi);PrintS(" to ");p_wrp(sit->p, currRing, strat->tailRing);PrintLn();
            }
            #endif
            if (TEST_OPT_PROT)
            {
              if (sit->p==NULL)
                PrintS("V");
              else
                PrintS("v");
              mflush();
            }
          }
          pLmDelete(&redSi);
          if (sit->p==NULL)
          {
            sit = strat->S.erase_and_next(sit);
            continue;
          }
          else if (change)
          {
            if (TEST_OPT_INTSTRATEGY)
            {
              if (TEST_OPT_CONTENTSB)
              {
                number n;
                p_Cleardenom_n(sit->p, currRing, n);// also does remove Content
                if (!nIsOne(n))
                {
                  denominator_list denom=(denominator_list)omAlloc(sizeof(denominator_list_s));
                  denom->n=nInvers(n);
                  denom->next=DENOMINATOR_LIST;
                  DENOMINATOR_LIST=denom;
                }
                nDelete(&n);
              }
              else
              {
                sit->p=p_Cleardenom(sit->p, currRing);// also does remove Content
              }
            }
            else
            {
              pNorm(sit->p);
            }
            sit->sev = pGetShortExpVector(sit->p);
          }
        }
        ++sit;
      }
      if (any_change) strat->S.reorder(&suc,strat);
      else break;
    }
    if (toT)
    {
      for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
      {
        if ((!strat->hasFromQ) || (sit->fromQ==0))
        {
          h.p = redtailBba(sit->p,sit,strat);
          if (TEST_OPT_INTSTRATEGY)
          {
            h.pCleardenom();// also does remove Content
          }
        }
        else
        {
          h.p = sit->p;
        }
        strat->initEcart(&h);
        if (strat->honey)
        {
          sit->ecart = h.ecart;
        }
        if (sit->sev == 0) {sit->sev = pGetShortExpVector(h.p);}
        else assume(sit->sev == pGetShortExpVector(h.p));
        h.sev = sit->sev;
        /*puts the elements of S also to T*/
        strat->initEcart(&h);
        /*if (toT) - already checked*/ enterT(h,strat);
        sit->s_2_r = strat->T.size()-1;
#ifdef HAVE_SHIFTBBA
        if (/*(toT) && */(currRing->isLPring))
          enterTShift(h, strat);
#endif
      }
    }
  }
  else
  {
    while (suc != -1)
    {
      auto sit = strat->S.begin();
      for (int skip = suc; skip > 0 && sit != strat->S.end(); --skip) ++sit;
      while (sit != strat->S.end())
      {
        change=FALSE;
        if (((!strat->hasFromQ) || (sit->fromQ==0)) && (sit != strat->S.begin()))
        {
          redSi=pHead(sit->p);
          sit->p = redMora(sit->p,sit,strat);
          if (sit->p==NULL)
          {
            sit = strat->S.erase_and_next(sit);
            pLmDelete(&redSi);
            kTest(strat);
            continue;
          }
          else if (pCmp(sit->p,redSi)!=0)
          {
            any_change=TRUE;
            h.p = sit->p;
            strat->initEcart(&h);
            sit->ecart = h.ecart;
            if (TEST_OPT_INTSTRATEGY)
            {
              if (TEST_OPT_CONTENTSB)
              {
                number n;
                p_Cleardenom_n(sit->p, currRing, n);// also does remove Content
                if (!nIsOne(n))
                {
                  denominator_list denom=(denominator_list)omAlloc(sizeof(denominator_list_s));
                  denom->n=nInvers(n);
                  denom->next=DENOMINATOR_LIST;
                  DENOMINATOR_LIST=denom;
                }
                nDelete(&n);
              }
              else
              {
                sit->p=p_Cleardenom(sit->p, currRing);// also does remove Content
              }
            }
            else
            {
              pNorm(sit->p); // == h.p
            }
            h.sev =  pGetShortExpVector(h.p);
            sit->sev = h.sev;
          }
          pLmDelete(&redSi);
          kTest(strat);
        }
        ++sit;
      }
#ifdef KDEBUG
      kTest(strat);
#endif
      if (any_change) strat->S.reorder(&suc,strat);
      else { suc=-1; break; }
      if (h.p!=NULL)
      {
        if (!strat->kAllAxis)
        {
          /*strat->kAllAxis =*/ HEckeTest(h.p,strat);
        }
        if (strat->kAllAxis)
          newHEdge(strat);
      }
    }
    for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
    {
      if ((!strat->hasFromQ) || (sit->fromQ==0))
      {
        sit->p = h.p = redtail(sit->p,strat->S.end(),strat);
        strat->initEcart(&h);
        sit->ecart = h.ecart;
        h.sev = pGetShortExpVector(h.p);
        sit->sev = h.sev;
      }
      else
      {
        h.p = sit->p;
        h.ecart=sit->ecart;
        h.sev = sit->sev;
        h.length = h.pLength = pLength(h.p);
      }
      if ((!strat->hasFromQ) || (sit->fromQ==0))
        cancelunit1(&h,&suc,strat->S.size()-1,strat);
      h.SetpFDeg();
      /*puts the elements of S also to T*/
      enterT(h,strat);
      sit->s_2_r = strat->T.size()-1;
#ifdef HAVE_SHIFTBBA
      if (currRing->isLPring)
        enterTShift(h, strat);
#endif
    }
    if (suc!= -1) updateS(toT,strat);
  }
#ifdef KDEBUG
  kTest(strat);
#endif
}


ideal skStrategy::getShdl()
{
  int n = S.size();
  if (n < 1) n = 1;  // always at least 1 slot (Singular convention for zero ideal)
  ideal I = idInit(n, Srank);
  int i = 0;
  for (auto sit = S.begin(); sit != S.end(); ++sit)
    I->m[i++] = sit->p;
  return I;
}

/*2
* -puts p to the standardbasis s at position at
* -saves the result in S
*/
sBasisSet::iterator enterSBba (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS)
{
  return strat->S.enter_bba(p, strat, atR, atS);
}

/*2
* sBasisSet::enter_bba — find sorted position and insert a new basis element.
* Replaces enterSBba for iterator-based usage.
* atS == end() means "compute position via find_pos internally"
* (was int atS == -1). Otherwise atS pins the insertion position.
*/
sBasisSet::iterator sBasisSet::enter_bba(LObject &p, kStrategy strat, int atR, iterator atS)
{
  strat->news = TRUE;
  ensure_capacity(size() + 1);

  int pos;
  if (atS != end())
    pos = atS.index();
  else if (order_ == SORDER_APPEND)
    pos = size();
  else
    pos = find_pos(p.p, p.ecart).index();

  SElement sobj;
  sobj.p = p.p;
  if (p.sev == 0)
    p.sev = pGetShortExpVector(p.p);
  else
    assume(p.sev == pGetShortExpVector(p.p));
  sobj.sev = p.sev;
  sobj.ecart = p.ecart;
  sobj.s_2_r = (atR >= 0) ? atR : p.i_r;
  sobj.length = 0;
  sobj.wlength = 0;
  sobj.fromQ = p.fromQ;
  // Capture arrival order.  Under parallel drain, this is taken from
  // strat->arrival_counter.fetch_add so that phase-1 iteration can
  // filter on "arrival_id < my_arrival".  Serial callers also bump the
  // counter to keep ordering monotonic.
  sobj.arrival_id = strat->arrival_counter.fetch_add(1, std::memory_order_relaxed);

  return insert_at(pos, sobj);
}

#ifdef HAVE_SHIFTBBA
sBasisSet::iterator enterSBbaShift (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS)
{
  auto it = strat->S.enter_bba(p, strat, atR, atS);

  int maxPossibleShift = p_mLPmaxPossibleShift(p.p, strat->tailRing);
  for (int i = maxPossibleShift; i > 0; i--)
  {
    LObject qq(p_Copy(p.p, strat->tailRing));
    p_mLPshift(qq.p, i, strat->tailRing);
    qq.shift = i;
    strat->initEcart(&qq);
    strat->S.enter_bba(qq, strat, -1, strat->S.end());
  }
  return it;
}
#endif

/*2
* -puts p to the standardbasis s at position at
* -saves the result in S
*/
sBasisSet::iterator enterSSba (LObject &p, kStrategy strat, int atR, sBasisSet::iterator atS)
{
  auto it = strat->S.enter_sba(p, strat, atR, atS);
#ifdef DEBUGF5
  Print("--- LIST S: %d ---\n",strat->S.size()-1);
  for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
  {
    pWrite(sit->sig);
  }
  PrintS("--- LIST S END ---\n");
#endif
  return it;
}

/*2
* sBasisSet::enter_sba — find sorted position and insert for signature-based algorithms.
* Replaces enterSSba for iterator-based usage.
*/
sBasisSet::iterator sBasisSet::enter_sba(LObject &p, kStrategy strat, int atR, iterator atS)
{
  strat->news = TRUE;
  ensure_capacity(size() + 1);

  int pos;
  if (atS != end())
    pos = atS.index();
  else if (order_ == SORDER_APPEND)
    pos = size();
  else
    pos = find_pos(p.p, p.ecart).index();

  SElement sobj;
  sobj.p = p.p;
  if (p.sev == 0)
    p.sev = pGetShortExpVector(p.p);
  else
    assume(p.sev == pGetShortExpVector(p.p));
  sobj.sev = p.sev;
  sobj.ecart = p.ecart;
  sobj.s_2_r = (atR >= 0) ? atR : p.i_r;
  sobj.length = 0;
  sobj.wlength = 0;
  sobj.fromQ = p.fromQ;
  sobj.arrival_id = strat->arrival_counter.fetch_add(1, std::memory_order_relaxed);

  // Original enterSSba inserts FIRST then sets sig/sevSig.
  auto it = insert_at(pos, sobj);
  elem(pos).sig = p.sig;
  if (p.sig != NULL)
  {
    if (p.sevSig == 0)
      p.sevSig = pGetShortExpVector(p.sig);
    else
      assume(p.sevSig == pGetShortExpVector(p.sig));
    elem(pos).sevSig = p.sevSig;
  }

#ifdef DEBUGF5
  int k;
  Print("--- LIST S: %d ---\n",size()-1);
  for(k=0;k < size();k++)
  {
    pWrite(elem(k).sig);
  }
  PrintS("--- LIST S END ---\n");
#endif
  return it;
}

void replaceInLAndSAndT(LObject &p, int tj, kStrategy strat)
{
  p.GetP(strat->lmBin);
  if (strat->homog) strat->initEcart(&p);
      strat->redTailChange=FALSE;
  if (TEST_OPT_INTSTRATEGY)
  {
    p.pCleardenom();
    if ((TEST_OPT_REDSB)||(TEST_OPT_REDTAIL))
    {
#ifdef HAVE_SHIFTBBA
      if (rIsLPRing(currRing))
        p.p = redtailBba(&p,strat->S.end(),strat, TRUE,!TEST_OPT_CONTENTSB);
      else
#endif
      {
        p.p = redtailBba(&p,strat->S.end(),strat, FALSE,!TEST_OPT_CONTENTSB);
      }
      p.pCleardenom();
      if (strat->redTailChange)
        p.t_p=NULL;
      if (strat->P.p!=NULL) strat->P.sev=p_GetShortExpVector(strat->P.p,currRing);
      else strat->P.sev=0;
    }
  }

  assume(strat->tailRing == p.tailRing);
  assume(p.pLength == 0 || pLength(p.p) == p.pLength || rIsSyzIndexRing(currRing)); // modulo syzring

  poly tp = strat->T[tj].p;

  /* enter p to T set */
  enterT(p, strat);

  auto sit_match = strat->S.end();
  for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
  {
    if (pLtCmp(tp, sit->p) == 0)
    {
      sit_match = sit;
      break;
    }
  }
  /* it may be that the exchanged element
   * is until now only in T and not in S */
  if (sit_match != strat->S.end())
  {
    strat->S.erase(sit_match);
  }

  auto pos = strat->S.find_pos(p.p, p.ecart);

  pp_Test(p.p, currRing, p.tailRing);
  assume(p.FDeg == p.pFDeg());

  /* remove useless pairs from L set */
  for (auto it = strat->L.begin(); it != strat->L.end(); )
  {
    if (it->p1 != NULL && pLtCmp(tp, it->p1) == 0)
    {
      it = strat->L.erase(it);
    }
    else if (it->p2 != NULL && pLtCmp(tp, it->p2) == 0)
    {
      it = strat->L.erase(it);
    }
    else
      it ++;
  }
#ifdef HAVE_SHIFTBBA
  if (rIsLPRing(currRing))
    enterpairsShift(p.p, strat->S.size()-1, p.ecart, pos, strat, strat->T.size()-1); // TODO LP
  else
#endif
  {
    /* generate new pairs with p, probably removing older, now useless pairs */
    superenterpairs(p.p, strat->S.size()-1, p.ecart, pos, strat, strat->T.size()-1);
  }
  /* enter p to S set */
  strat->enterS(p, strat, strat->T.size()-1, strat->S.end());

#ifdef HAVE_SHIFTBBA
  /* do this after enterS so that the index in R (which is strat->T.size()-1) is correct */
  if (rIsLPRing(currRing) && !strat->rightGB)
    enterTShift(p,strat);
#endif
}

/*2
* puts p to the set T at position atT
*/
void enterT(LObject &p, kStrategy strat, int atT)
{
  int i;

#ifdef PDEBUG
#ifdef HAVE_SHIFTBBA
  if (currRing->isLPring && p.shift > 0)
  {
    // in this case, the order is not correct. test LM and tail separately
    p_LmTest(p.p, currRing);
    p_Test(pNext(p.p), currRing);
  }
  else
#endif
  {
    pp_Test(p.p, currRing, p.tailRing);
  }
#endif
  assume(strat->tailRing == p.tailRing);
  // redMoraNF complains about this -- but, we don't really
  // need this so far
  assume(p.pLength == 0 || pLength(p.p) == p.pLength || rIsSyzIndexRing(currRing)); // modulo syzring
  assume(!strat->homog || (p.FDeg == p.pFDeg()));
  assume(!p.is_normalized || nIsOne(pGetCoeff(p.p)));

#ifdef KDEBUG
  // do not put an LObject twice into T:
  for(i=strat->T.size()-1;i>=0;i--)
  {
    if (p.p==strat->T[i].p)
    {
      printf("already in T at pos %d of %d, atT=%d\n",i,strat->T.size()-1,atT);
      return;
    }
  }
#endif

#ifdef HAVE_TAIL_RING
  if (currRing!=strat->tailRing)
  {
    p.t_p=p.GetLmTailRing();
  }
#endif
  strat->newt = TRUE;
  if (atT < 0)
    atT = strat->posInT(strat->T, strat->T.size()-1, p);

  // Probe (off-by-one pair-pop investigation): catch the case where
  // enterT targets a slot whose .p is already non-NULL (OVERWRITE)
  // or where the shift-up loop below will run (SHIFT).  Both
  // invalidate pair.i_r2 values that reference affected slots.
  {
    int pre_T_size = strat->T.size();
    poly pre_T_atT_p =
      (atT >= 0 && atT < pre_T_size) ? strat->T[atT].p : NULL;
    kt_debug_enterT_slot_probe(atT, pre_T_size,
                               (void*)pre_T_atT_p, (void*)p.p);
  }

  // Ensure capacity for tl+2 elements (current tl+1, plus the new one)
  strat->T.ensure_capacity(strat->T.size()-1 + 2);
  strat->sevT.ensure_capacity(strat->T.size()-1 + 2);
  strat->R.ensure_capacity(strat->T.size()-1 + 2);
  if (atT < strat->T.size())
  {
    for (i=strat->T.size(); i>=atT+1; i--)
    {
      strat->T[i] = strat->T[i-1];
      strat->sevT[i] = strat->sevT[i-1];
      // R-write probe: shift-loop rewrite.  Expected to be benign
      // (same ->p, new T-address); probe silent in that case.
      {
        int k = strat->T[i].i_r;
        TObject *old_target = (k >= 0) ? strat->R[k] : NULL;
        poly old_p = (old_target != NULL) ? old_target->p : NULL;
        TObject *new_target = strat->T.addr(i);
        poly new_p = (new_target != NULL) ? new_target->p : NULL;
        kt_debug_R_write_probe("enterT:shift", k,
                               (void*)old_target, (void*)new_target,
                               (void*)old_p, (void*)new_p);
      }
      strat->R[strat->T[i].i_r] = strat->T.addr(i);
    }
  }

  // Skip p_ShallowCopyDelete when tailBin == tailRing->PolyBin (always true
  // with --disable-omalloc, since omGetStickyBinOfBin is the identity).
  // This eliminates unnecessary malloc/memcpy/free of every monomial chain.
  if ((strat->tailBin != NULL)
      && (strat->tailRing == NULL || strat->tailBin != strat->tailRing->PolyBin)
      && (pNext(p.p) != NULL))
  {
    // PROBE: if we're here with --disable-omalloc, something's off —
    // tailBin should equal tailRing->PolyBin in that config.  Also
    // this path frees the source chain as it copies, so if p.p's
    // tail aliases with any live T entry's tail, that T entry gets
    // corrupted.
    extern void (*kbucket_debug_tag)(const char *op, void *lm, int slot, int arg);
    if (kbucket_debug_tag != NULL)
      kbucket_debug_tag("enterT:ShallowCopyDelete(active)",
                        (void*)pNext(p.p), -1, 0);
#ifdef HAVE_SHIFTBBA
    // letterplace: if p.shift > 0 then pNext(p.p) is already in the tailBin
    if (!(currRing->isLPring && p.shift > 0))
#endif
    {
      pNext(p.p)=p_ShallowCopyDelete(pNext(p.p),
          (strat->tailRing != NULL ?
           strat->tailRing : currRing),
          strat->tailBin);
      if (p.t_p != NULL) pNext(p.t_p) = pNext(p.p);
    }
  }
  strat->T[atT] = (TObject) p;
  //printf("\nenterT: add new: length = %i, ecart = %i\n",p.length,p.ecart);

  if ((pNext(p.p) != NULL) && (!rIsLPRing(currRing)))
    strat->T[atT].max_exp = p_GetMaxExpP(pNext(p.p), strat->tailRing);
  else
    strat->T[atT].max_exp = NULL;

  // Compute pLength inline (task 510 t-iterator-published).  Moves the
  // work that the parallel main-loop used to do in a post-drain refresh
  // loop (kthread.cc post-drain "pLength <= 0 -> set it") into enterT
  // itself.  The LObject-to-TObject assignment above copies p.pLength
  // (often 0 or -1 from sLObject::PrepareRed), so we fill it here if
  // still unset.  The next task removes the post-drain refresh loop
  // once all T readers have migrated to the iterator.
  if (strat->T[atT].pLength <= 0)
    strat->T[atT].pLength =
      ::pLength(strat->T[atT].p != NULL ? strat->T[atT].p
                                         : strat->T[atT].t_p);

  // Write sevT, R, and i_r BEFORE incrementing tl, so concurrent
  // readers (parallel bba workers) see fully initialized data when
  // they observe the new tl value.  The compiler barrier prevents
  // reordering of the tl++ past the data writes.
  assume((p.sev == 0) || (pGetShortExpVector(p.p) == p.sev));
  strat->sevT[atT] = (p.sev == 0 ? pGetShortExpVector(p.p) : p.sev);
  // R-write probe: fresh-slot write at end of enterT.  This is the
  // KEY SUSPECT for the remaining real bug: if strat->R[T.size()]
  // already held a non-NULL TObject whose .p differs from the new
  // T[atT].p, we're reassigning an R-slot that a prior pair still
  // references via its i_r2.
  {
    int k = (int)strat->T.size();
    TObject *old_target = (k >= 0) ? strat->R[k] : NULL;
    poly old_p = (old_target != NULL) ? old_target->p : NULL;
    TObject *new_target = strat->T.addr(atT);
    poly new_p = (new_target != NULL) ? new_target->p : NULL;
    kt_debug_R_write_probe("enterT:fresh_slot", k,
                           (void*)old_target, (void*)new_target,
                           (void*)old_p, (void*)new_p);
  }
  strat->R[strat->T.size()] = strat->T.addr(atT);
  strat->T[atT].i_r = strat->T.size();
  p.i_r = strat->T.size();  // propagate back so enterS can use it

  // Release-publish the slot (task 510 t-iterator-published).  All
  // field writes above must complete before any reader observes
  // published=true via the acquire-loading iterator.  In serial code
  // (THREADS=1 or main-thread drain) this is a no-op logically —
  // nobody else is reading — but keeps the invariant "T[i] published
  // iff T[i] fully initialised".  The existing __asm__ volatile
  // barrier below already orders this write with the later setsize();
  // the release store subsumes that ordering for the publish-observing
  // reader path.
  tobject_publish(strat->T[atT]);

  // Register every chain node of this new T entry into the T-node
  // registry so later pNext-writes that target these addresses fire
  // a "TNODE:mutation" tag, revealing the permanent-corruption
  // mutator.  Walk both representations (p in currRing, t_p in
  // tailRing) since they're independent allocations.
  for (poly q = strat->T[atT].p; q != NULL; q = pNext(q))
    kt_debug_register_tnode((void*)q, atT);
  for (poly q = strat->T[atT].t_p; q != NULL; q = pNext(q))
    kt_debug_register_tnode((void*)q, atT);

  // Snapshot T[atT].p so we can detect any post-enterT overwrite.
  kt_debug_snapshot_T_head(atT, (void*)strat->T[atT].p);

  __asm__ __volatile__("" ::: "memory");  // compiler barrier (x86 has strong HW ordering)
  strat->T.setsize(strat->T.size()+1);
  strat->R.setsize(strat->T.size());
  kTest_T(&(strat->T[atT]),strat);
}

/*2
* puts p to the set T at position atT
*/
void enterT_strong(LObject &p, kStrategy strat, int atT)
{
  assume(rField_is_Ring(currRing));
  int i;

  pp_Test(p.p, currRing, p.tailRing);
  assume(strat->tailRing == p.tailRing);
  // redMoraNF complains about this -- but, we don't really
  // need this so far
  assume(p.pLength == 0 || (int)pLength(p.p) == p.pLength || rIsSyzIndexRing(currRing)); // modulo syzring
  assume(p.FDeg == p.pFDeg());
  assume(!p.is_normalized || nIsOne(pGetCoeff(p.p)));

#ifdef KDEBUG
  // do not put an LObject twice into T:
  for(i=strat->T.size()-1;i>=0;i--)
  {
    if (p.p==strat->T[i].p)
    {
      printf("already in T at pos %d of %d, atT=%d\n",i,strat->T.size()-1,atT);
      return;
    }
  }
#endif

#ifdef HAVE_TAIL_RING
  if (currRing!=strat->tailRing)
  {
    p.t_p=p.GetLmTailRing();
  }
#endif
  strat->newt = TRUE;
  if (atT < 0)
    atT = strat->posInT(strat->T, strat->T.size()-1, p);
  // Ensure capacity for tl+2 elements
  strat->T.ensure_capacity(strat->T.size()-1 + 2);
  strat->sevT.ensure_capacity(strat->T.size()-1 + 2);
  strat->R.ensure_capacity(strat->T.size()-1 + 2);
  if (atT < strat->T.size())
  {
    for (i=strat->T.size(); i>=atT+1; i--)
    {
      strat->T[i] = strat->T[i-1];
      strat->sevT[i] = strat->sevT[i-1];
      // R-write probe: enterT_strong shift-loop rewrite.
      {
        int k = strat->T[i].i_r;
        TObject *old_target = (k >= 0) ? strat->R[k] : NULL;
        poly old_p = (old_target != NULL) ? old_target->p : NULL;
        TObject *new_target = strat->T.addr(i);
        poly new_p = (new_target != NULL) ? new_target->p : NULL;
        kt_debug_R_write_probe("enterT_strong:shift", k,
                               (void*)old_target, (void*)new_target,
                               (void*)old_p, (void*)new_p);
      }
      strat->R[strat->T[i].i_r] = strat->T.addr(i);
    }
  }

  if ((strat->tailBin != NULL) && (pNext(p.p) != NULL))
  {
    pNext(p.p)=p_ShallowCopyDelete(pNext(p.p),
                                   (strat->tailRing != NULL ?
                                    strat->tailRing : currRing),
                                   strat->tailBin);
    if (p.t_p != NULL) pNext(p.t_p) = pNext(p.p);
  }
  strat->T[atT] = (TObject) p;
  //printf("\nenterT_strong: add new: length = %i, ecart = %i\n",p.length,p.ecart);

  if (pNext(p.p) != NULL)
    strat->T[atT].max_exp = p_GetMaxExpP(pNext(p.p), strat->tailRing);
  else
    strat->T[atT].max_exp = NULL;

  // Compute pLength inline (task 510 t-iterator-published) — parallel
  // with enterT above.
  if (strat->T[atT].pLength <= 0)
    strat->T[atT].pLength =
      ::pLength(strat->T[atT].p != NULL ? strat->T[atT].p
                                         : strat->T[atT].t_p);

  strat->T.setsize(strat->T.size()+1);
  // R-write probe: enterT_strong fresh-slot write.
  {
    int k = (int)strat->T.size()-1;
    TObject *old_target = (k >= 0) ? strat->R[k] : NULL;
    poly old_p = (old_target != NULL) ? old_target->p : NULL;
    TObject *new_target = strat->T.addr(atT);
    poly new_p = (new_target != NULL) ? new_target->p : NULL;
    kt_debug_R_write_probe("enterT_strong:fresh_slot", k,
                           (void*)old_target, (void*)new_target,
                           (void*)old_p, (void*)new_p);
  }
  strat->R[strat->T.size()-1] = strat->T.addr(atT);
  strat->T[atT].i_r = strat->T.size()-1;
  assume(p.sev == 0 || pGetShortExpVector(p.p) == p.sev);
  strat->sevT[atT] = (p.sev == 0 ? pGetShortExpVector(p.p) : p.sev);
  // Release-publish the slot (task 510 t-iterator-published).
  tobject_publish(strat->T[atT]);
  #if 1
  if(rHasLocalOrMixedOrdering(currRing)
  && !n_IsUnit(p.p->coef, currRing->cf))
  {
    for(i=strat->T.size()-1;i>=0;i--)
    {
      if(strat->T[i].ecart <= p.ecart && pLmDivisibleBy(strat->T[i].p,p.p))
      {
        SElement tmp_si;
        tmp_si.p = strat->T[i].p;
        enterOneStrongPoly(tmp_si,p.p,p.ecart,0,strat,0 , TRUE);
      }
    }
  }
  /*
  printf("\nThis is T:\n");
  for(i=strat->T.size()-1;i>=0;i--)
  {
    pWrite(strat->T[i].p);
  }
  //getchar();*/
  #endif
  kTest_T(&(strat->T[atT]),strat);
}

/*2
* puts signature p.sig to the set syz
*/
void enterSyz(LObject &p, kStrategy strat, int atT)
{
  int i;
  strat->newt = TRUE;
  if (strat->syzl == strat->syzmax-1)
  {
    pEnlargeSet(&strat->syz,strat->syzmax,setmax);
    strat->sevSyz = (unsigned long*) omRealloc0Size(strat->sevSyz,
                                    (strat->syzmax)*sizeof(unsigned long),
                                    ((strat->syzmax)+setmax)
                                                  *sizeof(unsigned long));
    strat->syzmax += setmax;
  }
  if (atT < strat->syzl)
  {
#ifdef ENTER_USE_MEMMOVE
    memmove(&(strat->syz[atT+1]), &(strat->syz[atT]),
            (strat->syzl-atT+1)*sizeof(poly));
    memmove(&(strat->sevSyz[atT+1]), &(strat->sevSyz[atT]),
            (strat->syzl-atT+1)*sizeof(unsigned long));
#endif
    for (i=strat->syzl; i>=atT+1; i--)
    {
#ifndef ENTER_USE_MEMMOVE
      strat->syz[i] = strat->syz[i-1];
      strat->sevSyz[i] = strat->sevSyz[i-1];
#endif
    }
  }
  //i = strat->syzl;
  i = atT;
  //Makes sure the syz saves just the signature
  if(rField_is_Ring(currRing))
    pNext(p.sig) = NULL;
  strat->syz[atT] = p.sig;
  strat->sevSyz[atT] = p.sevSig;
  strat->syzl++;
#if F5DEBUG
  Print("element in strat->syz: %d--%d  ",atT+1,strat->syzmax);
  pWrite(strat->syz[atT]);
#endif
  // recheck pairs in strat->L with new syzygy rule and delete correspondingly.
  // Use the sevSig_flat_ array for cache-friendly pre-filtering: the filtered
  // iterator only visits elements whose sevSig passes the divisibility sev check,
  // then we do the full monomial divisibility test on those candidates.
  {
    unsigned long sev_syz = strat->sevSyz[atT];
    for (auto it = strat->L.ufbegin_sig(sev_syz); it != strat->L.ufend_sig(); )
    {
      if (p_LmDivisibleBy(strat->syz[atT], it->sig, currRing)
          &&((!rField_is_Ring(currRing))
          || (n_DivBy(pGetCoeff(it->sig),pGetCoeff(strat->syz[atT]),currRing->cf) && (pLtCmp(it->sig,strat->syz[atT])==1)))
          )
      {
        it = strat->L.erase(it);
      }
      else
        ++it;
    }
  }

//#if 1
#ifdef DEBUGF5
    PrintS("--- Syzygies ---\n");
    Print("syzl   %d\n",strat->syzl);
    Print("syzmax %d\n",strat->syzmax);
    PrintS("--------------------------------\n");
    for(i=0;i<=strat->syzl-1;i++)
    {
      Print("%d - ",i);
      pWrite(strat->syz[i]);
    }
    PrintS("--------------------------------\n");
#endif
}


void initHilbCrit(ideal/*F*/, ideal /*Q*/, bigintmat **hilb,kStrategy strat)
{

  //if the ordering is local, then hilb criterion
  //can be used also if the ideal is not homogeneous
  if((rHasLocalOrMixedOrdering(currRing)) && (rHasMixedOrdering(currRing)==FALSE))
  {
    if(rField_is_Ring(currRing))
      *hilb=NULL;
    else
      return;
  }
  if (strat->homog!=isHomog)
  {
    *hilb=NULL;
  }
}

void initBuchMoraCrit(kStrategy strat)
{
  strat->enterOnePair=enterOnePairNormal;
  strat->chainCrit=chainCritNormal;
  if (TEST_OPT_SB_1)
    strat->chainCrit=chainCritOpt_1;
  if (rField_is_Ring(currRing))
  {
    strat->enterOnePair=enterOnePairRing;
    strat->chainCrit=chainCritRing;
  }
#ifdef HAVE_RATGRING
  if (rIsRatGRing(currRing))
  {
     strat->chainCrit=chainCritPart;
     /* enterOnePairNormal get rational part in it */
  }
#endif
  if (TEST_OPT_IDLIFT
  && (strat->syzComp==1)
  && (!rIsPluralRing(currRing)))
    strat->enterOnePair=enterOnePairLift;

  strat->sugarCrit =        TEST_OPT_SUGARCRIT;
  strat->Gebauer =          strat->homog || strat->sugarCrit;
  strat->honey =            !strat->homog || strat->sugarCrit || TEST_OPT_WEIGHTM;
  if (TEST_OPT_NOT_SUGAR) strat->honey = FALSE;
  /* always use tailreduction, except:
  * - in local rings, - in lex order case, -in ring over extensions */
  strat->noTailReduction = !TEST_OPT_REDTAIL;
  //if(rHasMixedOrdering(currRing)==2)
  //{
  // strat->noTailReduction =TRUE;
  //}

#ifdef HAVE_PLURAL
  // and r is plural_ring
  //  hence this holds for r a rational_plural_ring
  if( rIsPluralRing(currRing) || (rIsSCA(currRing) && !strat->z2homog) )
  {    //or it has non-quasi-comm type... later
    strat->sugarCrit = FALSE;
    strat->Gebauer = FALSE;
    strat->honey = FALSE;
  }
#endif

  // Coefficient ring?
  if (rField_is_Ring(currRing))
  {
    strat->sugarCrit = FALSE;
    strat->Gebauer = FALSE;
    strat->honey = FALSE;
  }
  #ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    if (strat->homog) PrintS("ideal/module is homogeneous\n");
    else              PrintS("ideal/module is not homogeneous\n");
  }
  #endif
}

void initSbaCrit(kStrategy strat)
{
  //strat->enterOnePair=enterOnePairNormal;
  strat->enterOnePair = enterOnePairNormal;
  //strat->chainCrit=chainCritNormal;
  strat->chainCrit    = chainCritSig;
  /******************************************
   * rewCrit1 and rewCrit2 are already set in
   * kSba() in kstd1.cc
   *****************************************/
  //strat->rewCrit1     = faugereRewCriterion;
  if (strat->sbaOrder == 1)
  {
    strat->syzCrit  = syzCriterionInc;
  }
  else
  {
    strat->syzCrit  = syzCriterion;
  }
  if (rField_is_Ring(currRing))
  {
    strat->enterOnePair=enterOnePairRing;
    strat->chainCrit=chainCritRing;
  }
#ifdef HAVE_RATGRING
  if (rIsRatGRing(currRing))
  {
     strat->chainCrit=chainCritPart;
     /* enterOnePairNormal get rational part in it */
  }
#endif

  strat->sugarCrit =        TEST_OPT_SUGARCRIT;
  strat->Gebauer =          strat->homog || strat->sugarCrit;
  strat->honey =            !strat->homog || strat->sugarCrit || TEST_OPT_WEIGHTM;
  if (TEST_OPT_NOT_SUGAR) strat->honey = FALSE;
  /* always use tailreduction, except:
  * - in local rings, - in lex order case, -in ring over extensions */
  strat->noTailReduction = !TEST_OPT_REDTAIL;
  if(rHasMixedOrdering(currRing)) strat->noTailReduction =TRUE;

#ifdef HAVE_PLURAL
  // and r is plural_ring
  //  hence this holds for r a rational_plural_ring
  if( rIsPluralRing(currRing) || (rIsSCA(currRing) && !strat->z2homog) )
  {    //or it has non-quasi-comm type... later
    strat->sugarCrit = FALSE;
    strat->Gebauer = FALSE;
    strat->honey = FALSE;
  }
#endif

  // Coefficient ring?
  if (rField_is_Ring(currRing))
  {
    strat->sugarCrit = FALSE;
    strat->Gebauer = FALSE ;
    strat->honey = FALSE;
  }
  #ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    if (strat->homog) PrintS("ideal/module is homogeneous\n");
    else              PrintS("ideal/module is not homogeneous\n");
  }
  #endif
}

BOOLEAN kCompareLDependsOnLength(int (*compare_in_l)
                               (const LObject &lhs, const LObject &rhs, const kStrategy strat))
{
  if (compare_in_l == compareL110
  ||  compare_in_l == compareL10
  ||  compare_in_l == compareL110Ring
  )
    return TRUE;

  return FALSE;
}

void initBuchMoraPos (kStrategy strat)
{
  if (rHasGlobalOrdering(currRing))
  {
    if (strat->honey)
    {
      strat->compareL = compareL15;
      // ok -- here is the deal: from my experiments for Singular-2-0
      // I conclude that that posInT_EcartpLength is the best of
      // posInT15, posInT_EcartFDegpLength, posInT_FDegLength, posInT_pLength
      // see the table at the end of this file
      if (TEST_OPT_OLDSTD)
        strat->posInT = posInT15;
      else
        strat->posInT = posInT_EcartpLength;
    }
    else if (currRing->pLexOrder && !TEST_OPT_INTSTRATEGY)
    {
      strat->compareL = compareL11;
      strat->posInT = posInT11;
    }
    else if (TEST_OPT_INTSTRATEGY)
    {
      strat->compareL = compareL11;
      strat->posInT = posInT11;
    }
    else
    {
      strat->compareL = compareL0;
      strat->posInT = posInT0;
    }
    if (strat->homog)
    {
      strat->compareL = compareL110;
      strat->posInT = posInT110;
    }
  }
  else /* local/mixed ordering */
  {
    if (strat->homog)
    {
      strat->compareL = compareL11;
      strat->posInT = posInT11;
    }
    else
    {
      if ((currRing->order[0]==ringorder_c)
      ||(currRing->order[0]==ringorder_C))
      {
        strat->compareL = compareL17_c;
        strat->posInT = posInT17_c;
      }
      else
      {
        strat->compareL = compareL17;
        strat->posInT = posInT17;
      }
    }
  }
  if (strat->minim>0) strat->compareL = compareLSpecial;
  // for further tests only
  if ((BTEST1(11)) || (BTEST1(12)))
    strat->compareL = compareL11;
  else if ((BTEST1(13)) || (BTEST1(14)))
    strat->compareL = compareL13;
  else if ((BTEST1(15)) || (BTEST1(16)))
    strat->compareL = compareL15;
  else if ((BTEST1(17)) || (BTEST1(18)))
    strat->compareL = compareL17;
  if (BTEST1(11))
    strat->posInT = posInT11;
  else if (BTEST1(13))
    strat->posInT = posInT13;
  else if (BTEST1(15))
    strat->posInT = posInT15;
  else if ((BTEST1(17)))
    strat->posInT = posInT17;
  else if ((BTEST1(19)))
    strat->posInT = posInT19;
  else if (BTEST1(12) || BTEST1(14) || BTEST1(16) || BTEST1(18))
    strat->posInT = posInT1;
  strat->compareLDependsOnLength = kCompareLDependsOnLength(strat->compareL);
}

void initBuchMoraPosRing (kStrategy strat)
{
  if (rHasGlobalOrdering(currRing))
  {
    if (strat->honey)
    {
      strat->compareL = compareL15Ring;
      // ok -- here is the deal: from my experiments for Singular-2-0
      // I conclude that that posInT_EcartpLength is the best of
      // posInT15, posInT_EcartFDegpLength, posInT_FDegLength, posInT_pLength
      // see the table at the end of this file
      if (TEST_OPT_OLDSTD)
        strat->posInT = posInT15Ring;
      else
        strat->posInT = posInT_EcartpLength;
    }
    else if (currRing->pLexOrder && !TEST_OPT_INTSTRATEGY)
    {
      strat->compareL = compareL11Ring;
      strat->posInT = posInT11;
    }
    else if (TEST_OPT_INTSTRATEGY)
    {
      strat->compareL = compareL11Ring;
      strat->posInT = posInT11;
    }
    else
    {
      strat->compareL = compareL0Ring;
      strat->posInT = posInT0;
    }
    if (strat->homog)
    {
      strat->compareL = compareL110Ring;
      strat->posInT = posInT110Ring;
    }
  }
  else
  {
    if (strat->homog)
    {
      //printf("\nHere 3\n");
      strat->compareL = compareL11Ring;
      strat->posInT = posInT11Ring;
    }
    else
    {
      if ((currRing->order[0]==ringorder_c)
      ||(currRing->order[0]==ringorder_C))
      {
        strat->compareL = compareL17_cRing;
        strat->posInT = posInT17_cRing;
      }
      else
      {
        strat->compareL = compareL11Ringls;
        strat->posInT = posInT17Ring;
      }
    }
  }
  if (strat->minim>0) strat->compareL = compareLSpecial;
  // for further tests only
  if ((BTEST1(11)) || (BTEST1(12)))
    strat->compareL = compareL11Ring;
  else if ((BTEST1(13)) || (BTEST1(14)))
    strat->compareL = compareL13;
  else if ((BTEST1(15)) || (BTEST1(16)))
    strat->compareL = compareL15Ring;
  else if ((BTEST1(17)) || (BTEST1(18)))
    strat->compareL = compareL17Ring;
  if (BTEST1(11))
    strat->posInT = posInT11Ring;
  else if (BTEST1(13))
    strat->posInT = posInT13;
  else if (BTEST1(15))
    strat->posInT = posInT15Ring;
  else if ((BTEST1(17)))
    strat->posInT = posInT17Ring;
  else if ((BTEST1(19)))
    strat->posInT = posInT19;
  else if (BTEST1(12) || BTEST1(14) || BTEST1(16) || BTEST1(18))
    strat->posInT = posInT1;
  strat->compareLDependsOnLength = kCompareLDependsOnLength(strat->compareL);
}

void initBuchMora (ideal F,ideal Q,kStrategy strat)
{
  strat->interpt = BTEST1(OPT_INTERRUPT);
  /*- creating temp data structures------------------- -*/
  //strat->cp = 0; // already by skStragy()
  //strat->c3 = 0; // already by skStragy()
#ifdef HAVE_SHIFTBBA
  strat->cv = 0; // already by skStragy()
#endif
  strat->tail = pInit();
  /*- set s -*/
  strat->S.setsize(0);
  /*- set T -*/
  strat->T.setsize(0);
  initT(strat->T);
  initR(strat->R);
  initsevT(strat->sevT);
  /*- init local data struct.---------------------------------------- -*/
  //strat->P.ecart=0; // already by skStragy()
  //strat->P.length=0; // already by skStragy()
  //strat->P.pLength=0; // already by skStragy()
  if (rHasLocalOrMixedOrdering(currRing))
  {
    if (strat->kNoether!=NULL)
    {
      pSetComp(strat->kNoether, strat->ak);
      pSetComp(strat->kNoetherTail(), strat->ak);
    }
  }
  if(rField_is_Ring(currRing))
  {
    /*Shdl=*/initSL(F, Q,strat); /*sets also S, ecartS, fromQ */
  }
  else
  {
    if(TEST_OPT_SB_1)
    {
      int i;
      ideal P=idInit(IDELEMS(F)-strat->newIdeal,F->rank);
      for (i=strat->newIdeal;i<IDELEMS(F);i++)
      {
        P->m[i-strat->newIdeal] = F->m[i];
        F->m[i] = NULL;
      }
      initSSpecial(F,Q,P,strat);
      for (i=strat->newIdeal;i<IDELEMS(F);i++)
      {
        F->m[i] = P->m[i-strat->newIdeal];
        P->m[i-strat->newIdeal] = NULL;
      }
      idDelete(&P);
    }
    else
    {
      /*Shdl=*/initSL(F, Q,strat); /*sets also S, ecartS, fromQ */
      // /*Shdl=*/initS(F, Q,strat); /*sets also S, ecartS, fromQ */
    }
    if(errorreported) return;
  }
  strat->fromT = FALSE;
  strat->noTailReduction = !TEST_OPT_REDTAIL;
  if ((!TEST_OPT_SB_1)
  || (rField_is_Ring(currRing))
  )
  {
    updateS(TRUE,strat);
  }
#ifdef HAVE_SHIFTBBA
  if (!(rIsLPRing(currRing) && strat->rightGB)) // for right GB, we need to check later whether a poly is from Q
#endif
  {
    strat->hasFromQ=FALSE;
  }
  #ifdef KDEBUG
  assume(kTest_TS(strat));
  #endif
}

void exitBuchMora (kStrategy strat)
{
  /*- release temp data -*/
  // Dump L tombstone stats if requested (SINGULAR_LSET_STATS=1).  This is
  // the analogue of SINGULAR_SBASIS_STATS from task 503.
  if (getenv("SINGULAR_LSET_STATS")) {
    strat->L.debug_print_stats("exitBuchMora");
    strat_B(strat).debug_print_stats("exitBuchMora.B");
  }
  // Compact L before tearing down: forces poly cleanup on any tombstoned
  // entries left over from the last batch of chain-criterion erases.
  strat->L.compact();
  strat_B(strat).compact();
  cleanT(strat);
  strat->T.free_all();
  strat->R.free_all();
  strat->sevT.free_all();
  pLmFree(&strat->tail);
  strat->syzComp=0;

#ifdef HAVE_SHIFTBBA
  if (rIsLPRing(currRing) && strat->rightGB)
  {
    strat->hasFromQ=FALSE;
  }
#endif
}

void initSbaPos (kStrategy strat)
{
  if (rHasGlobalOrdering(currRing))
  {
    if (strat->honey)
    {
      if (TEST_OPT_OLDSTD)
        strat->posInT = posInT15;
      else
        strat->posInT = posInT_EcartpLength;
    }
    else if (currRing->pLexOrder && !TEST_OPT_INTSTRATEGY)
    {
      strat->posInT = posInT11;
    }
    else if (TEST_OPT_INTSTRATEGY)
    {
      strat->posInT = posInT11;
    }
    else
    {
      strat->posInT = posInT0;
    }
    if (strat->homog)
    {
      strat->posInT = posInT110;
    }
  }
  else
  {
    if (strat->homog)
    {
      strat->posInT = posInT11;
    }
    else
    {
      if ((currRing->order[0]==ringorder_c)
      ||(currRing->order[0]==ringorder_C))
      {
        strat->posInT = posInT17_c;
      }
      else
      {
        strat->posInT = posInT17;
      }
    }
  }
  // for further tests only
  if (BTEST1(11))
    strat->posInT = posInT11;
  else if (BTEST1(13))
    strat->posInT = posInT13;
  else if (BTEST1(15))
    strat->posInT = posInT15;
  else if ((BTEST1(17)))
    strat->posInT = posInT17;
  else if ((BTEST1(19)))
    strat->posInT = posInT19;
  else if (BTEST1(12) || BTEST1(14) || BTEST1(16) || BTEST1(18))
    strat->posInT = posInT1;
  if (rField_is_Ring(currRing))
  {
    strat->posInT = posInT11;
  }
  strat->compareLDependsOnLength = FALSE;
  strat->compareL = compareLSig;
}

void initSbaBuchMora (ideal F,ideal Q,kStrategy strat)
{
  strat->interpt = BTEST1(OPT_INTERRUPT);
  //strat->kNoether=NULL; // done by skStrategy
  /*- creating temp data structures------------------- -*/
  //strat->cp = 0; // done by skStrategy
  //strat->c3 = 0; // done by skStrategy
  strat->tail = pInit();
  /*- set s -*/
  strat->S.setsize(0);
  /*- set ps -*/
  strat->syzl = -1;
  /*- set T -*/
  strat->T.setsize(0);
  initT(strat->T);
  initR(strat->R);
  initsevT(strat->sevT);
  /*- init local data struct.---------------------------------------- -*/
  //strat->P.ecart=0;  // done by skStrategy
  //strat->P.length=0;  // done by skStrategy
  if (rHasLocalOrMixedOrdering(currRing))
  {
    if (strat->kNoether!=NULL)
    {
      pSetComp(strat->kNoether, strat->ak);
      pSetComp(strat->kNoetherTail(), strat->ak);
    }
  }
  if(rField_is_Ring(currRing))
  {
    /*Shdl=*/initSLSba(F, Q,strat); /*sets also S, ecartS, fromQ */
  }
  else
  {
    if(TEST_OPT_SB_1)
    {
      int i;
      ideal P=idInit(IDELEMS(F)-strat->newIdeal,F->rank);
      for (i=strat->newIdeal;i<IDELEMS(F);i++)
      {
        P->m[i-strat->newIdeal] = F->m[i];
        F->m[i] = NULL;
      }
      initSSpecialSba(F,Q,P,strat);
      for (i=strat->newIdeal;i<IDELEMS(F);i++)
      {
        F->m[i] = P->m[i-strat->newIdeal];
        P->m[i-strat->newIdeal] = NULL;
      }
      idDelete(&P);
    }
    else
    {
      initSLSba(F, Q,strat); /*sets also S, ecartS, fromQ */
    }
  }
  //strat->fromT = FALSE;  // done by skStrategy
  if (!TEST_OPT_SB_1)
  {
    if(!rField_is_Ring(currRing)) updateS(TRUE,strat);
  }
  //if (strat->hasFromQ) omFreeSize(strat->fromQ,strat->S.size()*sizeof(int));
  //strat->hasFromQ=FALSE;
  assume(kTest_TS(strat));
}

void exitSba (kStrategy strat)
{
  /*- release temp data -*/
  if(rField_is_Ring(currRing))
    cleanTSbaRing(strat);
  else
    cleanT(strat);
  strat->T.free_all();
  strat->R.free_all();
  strat->sevT.free_all();
  if(strat->syzmax>0)
  {
    omFreeSize((ADDRESS)strat->syz,(strat->syzmax)*sizeof(poly));
    omFreeSize((ADDRESS)strat->sevSyz,(strat->syzmax)*sizeof(unsigned long));
    if (strat->sbaOrder == 1)
    {
      omFreeSize(strat->syzIdx,(strat->syzidxmax)*sizeof(int));
    }
  }
  pLmDelete(&strat->tail);
  strat->syzComp=0;
}

/*2
* in the case of a standardbase of a module over a qring:
* replace polynomials in i by ak vectors,
* (the polynomial * unit vectors gen(1)..gen(ak)
* in every case (also for ideals:)
* deletes divisible vectors/polynomials
*/
void updateResult(ideal Q, kStrategy strat)
{
  if (strat->ak>0)
  {
    // Walk S in reverse raw-position order. These loops run before
    // final compaction; slots may have sit->p==NULL from earlier
    // pDelete in this same routine (semantic NULL, not tombstone).
    for (auto rit = strat->S.rbegin(); rit != strat->S.rend(); ++rit)
    {
      auto& se = *rit;
      if ((se.p!=NULL) && (pGetComp(se.p)==0))
      {
        pDelete(&se.p); // and set it to NULL
                      }
    }
    int q;
    poly p;
    if(!rField_is_Ring(currRing))
    {
      for (auto rit = strat->S.rbegin(); rit != strat->S.rend(); ++rit)
      {
        auto sit = rit.base(); --sit;
        if ((sit->p!=NULL)
        //&& (strat->syzComp>0)
        //&& (pGetComp(strat->S[l].p)<=strat->syzComp)
        )
        {
          for(q=IDELEMS(Q)-1; q>=0;q--)
          {
            if ((Q->m[q]!=NULL)
            &&(pLmDivisibleBy(Q->m[q],sit->p)))
            {
              if (TEST_OPT_REDSB)
              {
                p=sit->p;
                sit->p=kNF(Q,NULL,p);
                                pDelete(&p);
              }
              else
              {
                pDelete(&sit->p); // and set it to NULL
                              }
              break;
            }
          }
        }
      }
    }
    else
    {
      for (auto rit = strat->S.rbegin(); rit != strat->S.rend(); ++rit)
      {
        auto sit = rit.base(); --sit;
        if ((sit->p!=NULL)
        //&& (strat->syzComp>0)
        //&& (pGetComp(strat->S[l].p)<=strat->syzComp)
        )
        {
          for(q=IDELEMS(Q)-1; q>=0;q--)
          {
            if ((Q->m[q]!=NULL)
            &&(pLmDivisibleBy(Q->m[q],sit->p)))
            {
              if(n_DivBy(sit->p->coef, Q->m[q]->coef, currRing->cf))
              {
                if (TEST_OPT_REDSB)
                {
                  p=sit->p;
                  sit->p=kNF(Q,NULL,p);
                                    pDelete(&p);
                }
                else
                {
                  pDelete(&sit->p); // and set it to NULL
                                  }
                break;
              }
            }
          }
        }
      }
    }
  }
  else
  {
    int q;
    poly p;
    BOOLEAN reduction_found=FALSE;
    if (!rField_is_Ring(currRing))
    {
      for (auto rit = strat->S.rbegin(); rit != strat->S.rend(); ++rit)
      {
        auto sit = rit.base(); --sit;
        if (sit->p!=NULL)
        {
          for(q=IDELEMS(Q)-1; q>=0;q--)
          {
            if ((Q->m[q]!=NULL)&&(pLmEqual(Q->m[q],sit->p)))
            {
              if (TEST_OPT_REDSB)
              {
                p=sit->p;
                sit->p=kNF(Q,NULL,p);
                                pDelete(&p);
                reduction_found=TRUE;
              }
              else
              {
                pDelete(&sit->p); // and set it to NULL
                              }
              break;
            }
          }
        }
      }
    }
    //Also need divisibility of the leading coefficients
    else
    {
      for (auto rit = strat->S.rbegin(); rit != strat->S.rend(); ++rit)
      {
        auto sit = rit.base(); --sit;
        if (sit->p!=NULL)
        {
          for(q=IDELEMS(Q)-1; q>=0;q--)
          {
            if(n_DivBy(sit->p->coef, Q->m[q]->coef, currRing->cf))
            {
              if ((Q->m[q]!=NULL)&&(pLmEqual(Q->m[q],sit->p)) && pDivisibleBy(Q->m[q],sit->p))
              {
                if (TEST_OPT_REDSB)
                {
                  p=sit->p;
                  sit->p=kNF(Q,NULL,p);
                                    pDelete(&p);
                  reduction_found=TRUE;
                }
                else
                {
                  pDelete(&sit->p); // and set it to NULL
                                  }
                break;
              }
            }
          }
        }
      }
    }
    if (/*TEST_OPT_REDSB &&*/ reduction_found)
    {
      if(rField_is_Ring(currRing))
      {
        for (auto rit_l = strat->S.rbegin(); rit_l != strat->S.rend(); ++rit_l)
        {
          auto sit_l = rit_l.base(); --sit_l;
          if (sit_l->p!=NULL)
          {
            for (auto rit_q = strat->S.rbegin(); rit_q != strat->S.rend(); ++rit_q)
            {
              auto sit_q = rit_q.base(); --sit_q;
              if ((sit_l.index() != sit_q.index())
              && (sit_q->p!=NULL)
              &&(pLmDivisibleBy(sit_l->p,sit_q->p))
              &&(n_DivBy(sit_q->p->coef, sit_l->p->coef, currRing->cf))
              )
              {
                //If they are equal then take the one with the smallest length
                if(pLmDivisibleBy(sit_q->p,sit_l->p)
                && n_DivBy(sit_q->p->coef, sit_l->p->coef, currRing->cf)
                && (pLength(sit_q->p) < pLength(sit_l->p) ||
                (pLength(sit_q->p) == pLength(sit_l->p) && nGreaterZero(sit_q->p->coef))))
                {
                  pDelete(&sit_l->p);
                                    break;
                }
                else
                {
                  pDelete(&sit_q->p);
                                  }
              }
            }
          }
        }
      }
      else
      {
        for (auto rit_l = strat->S.rbegin(); rit_l != strat->S.rend(); ++rit_l)
        {
          auto sit_l = rit_l.base(); --sit_l;
          if (sit_l->p!=NULL)
          {
            for (auto rit_q = strat->S.rbegin(); rit_q != strat->S.rend(); ++rit_q)
            {
              auto sit_q = rit_q.base(); --sit_q;
              if ((sit_l.index() != sit_q.index())
              && (sit_q->p!=NULL)
              &&(pLmDivisibleBy(sit_l->p,sit_q->p))
              )
              {
                //If they are equal then take the one with the smallest length
                if(pLmDivisibleBy(sit_q->p,sit_l->p)
                &&(pLength(sit_q->p) < pLength(sit_l->p) ||
                (pLength(sit_q->p) == pLength(sit_l->p) && nGreaterZero(sit_q->p->coef))))
                {
                  pDelete(&sit_l->p);
                                    break;
                }
                else
                {
                  pDelete(&sit_q->p);
                                  }
              }
            }
          }
        }
      }
    }
  }
  // Compact S: remove NULL gaps left by pDelete, matching old idSkipZeroes
  // behavior. Bba continues using S after updateResult returns.
  strat->S.compact_null_p();
}

void completeReduce (kStrategy strat, BOOLEAN withT)
{
  int i;
  int low = (((rHasGlobalOrdering(currRing)) && (strat->ak==0)) ? 1 : 0);
  LObject L;

#ifdef KDEBUG
  // need to set this: during tailreductions of T[i], T[i].max is out of
  // sync
  sloppy_max = TRUE;
#endif

  strat->noTailReduction = FALSE;
  //if(rHasMixedOrdering(currRing)) strat->noTailReduction = TRUE;
  if (TEST_OPT_PROT)
  {
    PrintLn();
//    if (timerv) writeTime("standard base computed:");
  }
  if (TEST_OPT_PROT)
  {
    Print("(S:%d)",strat->S.size()-1);mflush();
  }
  // Reverse walk S in raw-position order; post-compaction, no tombstones.
  for (auto rit = strat->S.rbegin(); rit != strat->S.rend() && rit.index() >= low; ++rit)
  {
    i = rit.index();
    auto sit = rit.base(); --sit;  // forward iterator at the current element
    // Exclusive end iterator for redtailBba:
    //   strat->ak == 0  -> reduce against [0..i-1] inclusive = [begin, sit) exclusive
    //   otherwise       -> reduce against [0..size-1] inclusive = [begin, end()) exclusive
    sBasisSet::const_iterator end_it = strat->S.end();
    if ((strat->hasFromQ) && (sit->fromQ)) continue; // do not reduce Q_i
    if (strat->ak==0) end_it = sit;
    TObject* T_j = strat->S.s_2_t(sit, strat);
    if ((T_j != NULL)&&(T_j->p==sit->p))
    {
      L = *T_j;
      #ifdef KDEBUG
      if (TEST_OPT_DEBUG)
      {
        Print("test S[%d]:",i);
        p_wrp(L.p,currRing,strat->tailRing);
        PrintLn();
      }
      #endif
      if (rHasGlobalOrdering(currRing))
      {
        sit->p = redtailBba(&L, end_it, strat, withT,FALSE /*no normalize*/);
      }
      else
      {
        sit->p = redtail(&L, strat->S.end(), strat);
      }
      #ifdef KDEBUG
      if (TEST_OPT_DEBUG)
      {
        Print("to (tailR) S[%d]:",i);
        p_wrp(sit->p,currRing,strat->tailRing);
        PrintLn();
      }
      #endif

      if (strat->redTailChange)
      {
        if (T_j->max_exp != NULL) p_LmFree(T_j->max_exp, strat->tailRing);
        if (pNext(T_j->p) != NULL)
          T_j->max_exp = p_GetMaxExpP(pNext(T_j->p), strat->tailRing);
        else
          T_j->max_exp = NULL;
      }
      if (TEST_OPT_INTSTRATEGY)
        T_j->pCleardenom();
    }
    else
    {
      assume(currRing == strat->tailRing);
      #ifdef KDEBUG
      if (TEST_OPT_DEBUG)
      {
        Print("test S[%d]:",i);
        p_wrp(sit->p,currRing,strat->tailRing);
        PrintLn();
      }
      #endif
      if (rHasGlobalOrdering(currRing))
      {
        sit->p = redtailBba(sit->p, end_it, strat, withT);
      }
      else
      {
        sit->p = redtail(sit->p, strat->S.end(), strat);
      }
      if (TEST_OPT_INTSTRATEGY)
      {
        if (TEST_OPT_CONTENTSB)
        {
          number n;
          p_Cleardenom_n(sit->p, currRing, n);// also does remove Content
          if (!nIsOne(n))
          {
            denominator_list denom=(denominator_list)omAlloc(sizeof(denominator_list_s));
            denom->n=nInvers(n);
            denom->next=DENOMINATOR_LIST;
            DENOMINATOR_LIST=denom;
          }
          nDelete(&n);
        }
        else
        {
          sit->p=p_Cleardenom(sit->p, currRing);// also does remove Content
        }
      }
      #ifdef KDEBUG
      if (TEST_OPT_DEBUG)
      {
        Print("to (-tailR) S[%d]:",i);
        p_wrp(sit->p,currRing,strat->tailRing);
        PrintLn();
      }
      #endif
    }
    if (TEST_OPT_PROT)
      PrintS("-");
  }
  if (TEST_OPT_PROT) PrintLn();
#ifdef KDEBUG
  sloppy_max = FALSE;
#endif
}


/*2
* computes the new strat->kNoether and the new pNoether,
* returns TRUE, if pNoether has changed
*/
BOOLEAN newHEdge(kStrategy strat)
{
  if (currRing->pLexOrder || rHasMixedOrdering(currRing))
    return FALSE;
  int i,j;
  poly newNoether;

#if 0
  if (currRing->weight_all_1)
    scComputeHC(strat->getShdl(),NULL,strat->ak,strat->kNoether);
  else
    scComputeHCw(strat->getShdl(),NULL,strat->ak,strat->kNoether);
#else
  scComputeHC(strat->getShdl(),NULL,strat->ak,strat->kNoether);
#endif
  if (strat->kNoether==NULL) return FALSE;
  if (strat->t_kNoether != NULL)
  {
    p_LmFree(strat->t_kNoether, strat->tailRing);
    strat->t_kNoether=NULL;
  }
  if (strat->tailRing != currRing)
    strat->t_kNoether = k_LmInit_currRing_2_tailRing(strat->kNoether, strat->tailRing);
  /* compare old and new noether*/
  newNoether = pLmInit(strat->kNoether);
  pSetCoeff0(newNoether,nInit(1));
  j = p_FDeg(newNoether,currRing);
  // newNoether is now the new highest edge (on the boundary)
  // now find the next monomial after highest monomial in R/I
  // (the smallest monomial not in R/I)
  for (i=currRing->N-1;i>0; i--)
  {
    int e;
    if ((e=pGetExp(newNoether, i)) > 0)
    {
      e--;
      pSetExp(newNoether,i,e);
    }
  }
  pSetm(newNoether);
  if (j < HCord) /*- statistics -*/
  {
    if (TEST_OPT_PROT)
    {
      Print("H(%d)",j);
      mflush();
    }
    HCord=j;
    #ifdef KDEBUG
    if (TEST_OPT_DEBUG)
    {
      Print("H(%d):",j);
      wrp(strat->kNoether);
      PrintLn();
    }
    #endif
  }
  if (pCmp(strat->kNoether,newNoether)!=1)
  {
    if (strat->kNoether!=NULL) p_LmDelete0(strat->kNoether,currRing);
    strat->kNoether=newNoether;
    if (strat->t_kNoether != NULL)
    {
      p_LmFree(strat->t_kNoether, strat->tailRing);
      strat->t_kNoether=NULL;
    }
    if (strat->tailRing != currRing)
      strat->t_kNoether = k_LmInit_currRing_2_tailRing(strat->kNoether, strat->tailRing);

    return TRUE;
  }
  pLmDelete(newNoether);
  return FALSE;
}

/***************************************************************
 *
 * Routines related for ring changes during std computations
 *
 ***************************************************************/
BOOLEAN kCheckSpolyCreation(LObject *L, kStrategy strat, poly &m1, poly &m2)
{
  if (strat->overflow) return FALSE;
  assume(L->p1 != NULL && L->p2 != NULL);
  // shift changes: from 0 to -1
  assume(L->i_r1 >= -1 && L->i_r1 < strat->T.size());
  assume(L->i_r2 >= -1 && L->i_r2 < strat->T.size());

  if (! k_GetLeadTerms(L->p1, L->p2, currRing, m1, m2, strat->tailRing))
    return FALSE;
  // shift changes: extra case inserted
  if ((L->i_r1 == -1) || (L->i_r2 == -1) )
  {
    return TRUE;
  }
  poly p1_max=NULL;
  if ((L->i_r1>=0)&&(strat->R[L->i_r1]!=NULL)) p1_max = (strat->R[L->i_r1])->max_exp;
  poly p2_max=NULL;
  if ((L->i_r2>=0)&&(strat->R[L->i_r2]!=NULL)) p2_max = (strat->R[L->i_r2])->max_exp;

  if (((p1_max != NULL) && !p_LmExpVectorAddIsOk(m1, p1_max, strat->tailRing)) ||
      ((p2_max != NULL) && !p_LmExpVectorAddIsOk(m2, p2_max, strat->tailRing)))
  {
    p_LmFree(m1, strat->tailRing);
    p_LmFree(m2, strat->tailRing);
    m1 = NULL;
    m2 = NULL;
    return FALSE;
  }
  return TRUE;
}

/***************************************************************
 *
 * Checks, if we can compute the gcd poly / strong pair
 * gcd-poly = m1 * R[atR] + m2 * S[atS]
 *
 ***************************************************************/
BOOLEAN kCheckStrongCreation(int atR, poly m1, sBasisSet::const_iterator atS, poly m2, kStrategy strat)
{
  assume(atS->s_2_r >= -1 && atS->s_2_r < strat->T.size());
  //assume(strat->tailRing != currRing);

  poly p1_max = (strat->R[atR])->max_exp;
  poly p2_max = (strat->R[atS->s_2_r])->max_exp;

  if (((p1_max != NULL) && !p_LmExpVectorAddIsOk(m1, p1_max, strat->tailRing)) ||
      ((p2_max != NULL) && !p_LmExpVectorAddIsOk(m2, p2_max, strat->tailRing)))
  {
    return FALSE;
  }
  return TRUE;
}

/*!
  used for GB over ZZ: look for constant and monomial elements in the ideal
  background: any known constant element of ideal suppresses
              intermediate coefficient swell
*/
poly preIntegerCheck(const ideal Forig, const ideal Q)
{
  assume(nCoeff_is_Z(currRing->cf));
  ideal F = idCopy(Forig);
  idSkipZeroes(F);
  poly pmon;
  ring origR = currRing;
  ideal monred = idInit(1,1);
  for(int i=0; i<idElem(F); i++)
  {
    if(pNext(F->m[i]) == NULL)
        idInsertPoly(monred, pCopy(F->m[i]));
  }
  int posconst = idPosConstant(F);
  if((posconst != -1) && (!nIsZero(F->m[posconst]->coef)))
  {
    idDelete(&F);
    idDelete(&monred);
    return NULL;
  }
  int idelemQ = 0;
  if(Q!=NULL)
  {
    idelemQ = IDELEMS(Q);
    for(int i=0; i<idelemQ; i++)
    {
      if(pNext(Q->m[i]) == NULL)
        idInsertPoly(monred, pCopy(Q->m[i]));
    }
    idSkipZeroes(monred);
    posconst = idPosConstant(monred);
    //the constant, if found, will be from Q
    if((posconst != -1) && (!nIsZero(monred->m[posconst]->coef)))
    {
      pmon = pCopy(monred->m[posconst]);
      idDelete(&F);
      idDelete(&monred);
      return pmon;
    }
  }
  ring QQ_ring = rCopy0(currRing,FALSE);
  nKillChar(QQ_ring->cf);
  QQ_ring->cf = nInitChar(n_Q, NULL);
  rComplete(QQ_ring,1);
  QQ_ring = rAssure_c_dp(QQ_ring);
  rChangeCurrRing(QQ_ring);
  nMapFunc nMap = n_SetMap(origR->cf, QQ_ring->cf);
  ideal II = idInit(IDELEMS(F)+idelemQ+2,id_RankFreeModule(F, origR));
  for(int i = 0, j = 0; i<IDELEMS(F); i++)
    II->m[j++] = prMapR(F->m[i], nMap, origR, QQ_ring);
  for(int i = 0, j = IDELEMS(F); i<idelemQ; i++)
    II->m[j++] = prMapR(Q->m[i], nMap, origR, QQ_ring);
  ideal one = kStd2(II, NULL, isNotHomog, NULL,(bigintmat*)NULL);
  idSkipZeroes(one);
  if(idIsConstant(one))
  {
    //one should be <1>
    for(int i = IDELEMS(II)-1; i>=0; i--)
      if(II->m[i] != NULL)
        II->m[i+1] = II->m[i];
    II->m[0] = pOne();
    ideal syz = idSyzygies(II, isNotHomog, NULL);
    poly integer = NULL;
    for(int i = IDELEMS(syz)-1;i>=0; i--)
    {
      if(pGetComp(syz->m[i]) == 1)
      {
        pSetComp(syz->m[i],0);
        if(pIsConstant(pHead(syz->m[i])))
        {
          integer = pHead(syz->m[i]);
          break;
        }
      }
    }
    rChangeCurrRing(origR);
    nMapFunc nMap2 = n_SetMap(QQ_ring->cf, origR->cf);
    pmon = prMapR(integer, nMap2, QQ_ring, origR);
    idDelete(&monred);
    idDelete(&F);
    id_Delete(&II,QQ_ring);
    id_Delete(&one,QQ_ring);
    id_Delete(&syz,QQ_ring);
    p_Delete(&integer,QQ_ring);
    rDelete(QQ_ring);
    return pmon;
  }
  else
  {
    if(idIs0(monred))
    {
      poly mindegmon = NULL;
      for(int i = 0; i<IDELEMS(one); i++)
      {
        if(pNext(one->m[i]) == NULL)
        {
          if(mindegmon == NULL)
            mindegmon = pCopy(one->m[i]);
          else
          {
            if(p_Deg(one->m[i], QQ_ring) < p_Deg(mindegmon, QQ_ring))
              mindegmon = pCopy(one->m[i]);
          }
        }
      }
      if(mindegmon != NULL)
      {
        for(int i = IDELEMS(II)-1; i>=0; i--)
          if(II->m[i] != NULL)
            II->m[i+1] = II->m[i];
        II->m[0] = pCopy(mindegmon);
        ideal syz = idSyzygies(II, isNotHomog, NULL);
        bool found = FALSE;
        for(int i = IDELEMS(syz)-1;i>=0; i--)
        {
          if(pGetComp(syz->m[i]) == 1)
          {
            pSetComp(syz->m[i],0);
            if(pIsConstant(pHead(syz->m[i])))
            {
              pSetCoeff(mindegmon, nCopy(syz->m[i]->coef));
              found = TRUE;
              break;
            }
          }
        }
        id_Delete(&syz,QQ_ring);
        if (found == FALSE)
        {
          rChangeCurrRing(origR);
          idDelete(&monred);
          idDelete(&F);
          id_Delete(&II,QQ_ring);
          id_Delete(&one,QQ_ring);
          rDelete(QQ_ring);
          return NULL;
        }
        rChangeCurrRing(origR);
        nMapFunc nMap2 = n_SetMap(QQ_ring->cf, origR->cf);
        pmon = prMapR(mindegmon, nMap2, QQ_ring, origR);
        idDelete(&monred);
        idDelete(&F);
        id_Delete(&II,QQ_ring);
        id_Delete(&one,QQ_ring);
        id_Delete(&syz,QQ_ring);
        rDelete(QQ_ring);
        return pmon;
      }
    }
  }
  rChangeCurrRing(origR);
  idDelete(&monred);
  idDelete(&F);
  id_Delete(&II,QQ_ring);
  id_Delete(&one,QQ_ring);
  rDelete(QQ_ring);
  return NULL;
}

/*!
  used for GB over ZZ: intermediate reduction by monomial elements
  background: any known constant element of ideal suppresses
              intermediate coefficient swell
*/
void postReduceByMon(LObject* h, kStrategy strat)
{
  if(!nCoeff_is_Z(currRing->cf))
      return;
  poly pH = h->GetP();
  poly p,pp;
  p = pH;
  bool deleted = FALSE, ok = FALSE;
  for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
  {
    p = pH;
    if(pNext(sit->p) == NULL)
    {
      //pWrite(p);
      //pWrite(sit->p);
      while(ok == FALSE && p != NULL)
      {
        if(pLmDivisibleBy(sit->p, p)
#ifdef HAVE_SHIFTBBA
            || (rIsLPRing(currRing) && pLPLmDivisibleBy(sit->p, p))
#endif
          )
        {
          number dummy = n_IntMod(p->coef, sit->p->coef, currRing->cf);
          p_SetCoeff(p,dummy,currRing);
        }
        if(nIsZero(p->coef))
        {
          pLmDelete(&p);
          h->p = p;
          deleted = TRUE;
        }
        else
        {
          ok = TRUE;
        }
      }
      if (p!=NULL)
      {
        pp = pNext(p);
        while(pp != NULL)
        {
          if(pLmDivisibleBy(sit->p, pp)
#ifdef HAVE_SHIFTBBA
            || (rIsLPRing(currRing) && pLPLmDivisibleBy(sit->p, pp))
#endif
            )
          {
            number dummy = n_IntMod(pp->coef, sit->p->coef, currRing->cf);
            p_SetCoeff(pp,dummy,currRing);
            if(nIsZero(pp->coef))
            {
              pLmDelete(&pNext(p));
              pp = pNext(p);
              deleted = TRUE;
            }
            else
            {
              p = pp;
              pp = pNext(p);
            }
          }
          else
          {
            p = pp;
            pp = pNext(p);
          }
        }
      }
    }
  }
  h->SetLmCurrRing();
  if((deleted)&&(h->p!=NULL))
    strat->initEcart(h);
}

void postReduceByMonSig(LObject* h, kStrategy strat)
{
  if(!nCoeff_is_Z(currRing->cf))
      return;
  poly hSig = h->sig;
  poly pH = h->GetP();
  poly p,pp;
  p = pH;
  bool deleted = FALSE, ok = FALSE;
  for (auto sit = strat->S.begin(); sit != strat->S.end(); ++sit)
  {
    p = pH;
    if(pNext(sit->p) == NULL)
    {
      while(ok == FALSE && p!=NULL)
      {
        if(pLmDivisibleBy(sit->p, p))
        {
          poly sigMult = pDivideM(pHead(p),pHead(sit->p));
          sigMult = ppMult_mm(sigMult,pCopy(sit->sig));
          if(sigMult!= NULL && pLtCmp(hSig,sigMult) == 1)
          {
            number dummy = n_IntMod(p->coef, sit->p->coef, currRing->cf);
            p_SetCoeff(p,dummy,currRing);
          }
          pDelete(&sigMult);
        }
        if(nIsZero(p->coef))
        {
          pLmDelete(&p);
          h->p = p;
          deleted = TRUE;
        }
        else
        {
          ok = TRUE;
        }
      }
      if(p == NULL)
        return;
      pp = pNext(p);
      while(pp != NULL)
      {
        if(pLmDivisibleBy(sit->p, pp))
        {
          poly sigMult = pDivideM(pHead(p),pHead(sit->p));
          sigMult = ppMult_mm(sigMult,pCopy(sit->sig));
          if(sigMult!= NULL && pLtCmp(hSig,sigMult) == 1)
          {
            number dummy = n_IntMod(pp->coef, sit->p->coef, currRing->cf);
            p_SetCoeff(pp,dummy,currRing);
            if(nIsZero(pp->coef))
            {
              pLmDelete(&pNext(p));
              pp = pNext(p);
              deleted = TRUE;
            }
            else
            {
              p = pp;
              pp = pNext(p);
            }
          }
          else
          {
            p = pp;
            pp = pNext(p);
          }
          pDelete(&sigMult);
        }
        else
        {
          p = pp;
          pp = pNext(p);
        }
      }
    }
  }
  h->SetLmCurrRing();
  if(deleted)
    strat->initEcart(h);

}

/*!
  used for GB over ZZ: final reduction by constant elements
  background: any known constant element of ideal suppresses
              intermediate coefficient swell and beautifies output
*/
void finalReduceByMon(kStrategy strat)
{
  assume(strat->T.empty()); /* can only be called with no elements in T:
                          i.e. after exitBuchMora */
  /* do not use strat->S, strat->S.size()-1 as they may be out of sync*/
  if(!nCoeff_is_Z(currRing->cf))
      return;
  poly p,pp;
  for (auto sit_j = strat->S.begin(); sit_j != strat->S.end(); ++sit_j)
  {
    if((sit_j->p!=NULL)&&(pNext(sit_j->p) == NULL))
    {
      for (auto sit_i = strat->S.begin(); sit_i != strat->S.end(); ++sit_i)
      {
        if((sit_i != sit_j) && (sit_i->p != NULL))
        {
          p = sit_i->p;
          while((p!=NULL) && (pLmDivisibleBy(sit_j->p, p)
#if HAVE_SHIFTBBA
                || (rIsLPRing(currRing) && pLPLmDivisibleBy(sit_j->p, p))
#endif
                ))
          {
            number dummy = n_IntMod(p->coef, sit_j->p->coef, currRing->cf);
            if (!nEqual(dummy,p->coef))
            {
              if (nIsZero(dummy))
              {
                nDelete(&dummy);
                pLmDelete(&sit_i->p);
                p=sit_i->p;
              }
              else
              {
                p_SetCoeff(p,dummy,currRing);
                break;
              }
            }
            else
            {
              nDelete(&dummy);
              break;
            }
          }
          if (p!=NULL)
          {
            pp = pNext(p);
            while(pp != NULL)
            {
              if(pLmDivisibleBy(sit_j->p, pp)
#if HAVE_SHIFTBBA
                  || (rIsLPRing(currRing) && pLPLmDivisibleBy(sit_j->p, pp))
#endif
                )
              {
                number dummy = n_IntMod(pp->coef, sit_j->p->coef, currRing->cf);
                if (!nEqual(dummy,pp->coef))
                {
                  p_SetCoeff(pp,dummy,currRing);
                  if(nIsZero(pp->coef))
                  {
                    pLmDelete(&pNext(p));
                    pp = pNext(p);
                  }
                  else
                  {
                    p = pp;
                    pp = pNext(p);
                  }
                }
                else
                {
                  nDelete(&dummy);
                  p = pp;
                  pp = pNext(p);
                }
              }
              else
              {
                p = pp;
                pp = pNext(p);
              }
            }
          }
        }
      }
      //idPrint removed (no Shdl)
    }
  }
}


BOOLEAN kStratChangeTailRing(kStrategy strat, LObject *L, TObject* T, unsigned long expbound)
{
  // Audit task 308 suspect #2: bulk T rewrite on exp-bound overflow.
  // Mutates every T entry's t_p chain without a lock.  If this fires
  // during parallel phase, concurrent readers see torn T entries.
  kt_debug_tag("kStratChangeTailRing:entry", (void*)strat,
               strat->T.size(), 0);
  assume((strat->tailRing == currRing) || (strat->tailRing->bitmask <= currRing->bitmask));
  /* initial setup or extending */

  if (rIsLPRing(currRing)) return TRUE;
  if (expbound == 0) expbound = strat->tailRing->bitmask << 1;
  if (expbound >= currRing->bitmask) return FALSE;
  strat->overflow=FALSE;
  ring new_tailRing = rModifyRing(currRing,
  // Hmmm .. the condition pFDeg == p_Deg
  // might be too strong
  (strat->homog && currRing->pFDeg == p_Deg && !(rField_is_Ring(currRing))), // omit degree
  (strat->ak==0), // omit_comp if the input is an ideal
  expbound); // exp_limit

  if (new_tailRing == currRing) return TRUE;

  strat->pOrigFDeg_TailRing = new_tailRing->pFDeg;
  strat->pOrigLDeg_TailRing = new_tailRing->pLDeg;

  if (currRing->pFDeg != currRing->pFDegOrig)
  {
    new_tailRing->pFDeg = currRing->pFDeg;
    new_tailRing->pLDeg = currRing->pLDeg;
  }

  if (TEST_OPT_PROT)
    Print("[%lu:%d", (unsigned long) new_tailRing->bitmask, new_tailRing->ExpL_Size);
  #ifdef KDEBUG
  kTest_TS(strat);
  assume(new_tailRing != strat->tailRing);
  #endif
  pShallowCopyDeleteProc p_shallow_copy_delete
    = pGetShallowCopyDeleteProc(strat->tailRing, new_tailRing);

  omBin new_tailBin = omGetStickyBinOfBin(new_tailRing->PolyBin);

  int i;
  for (i=0; i < strat->T.size(); i++)
  {
    strat->T[i].ShallowCopyDelete(new_tailRing, new_tailBin,
                                  p_shallow_copy_delete);
  }
  for (auto& Lp: strat->L) {
    assume(Lp.p != NULL);
    if (pNext(Lp.p) != strat->tail)
      Lp.ShallowCopyDelete(new_tailRing, p_shallow_copy_delete);
  }
  if ((strat->P.t_p != NULL) ||
      ((strat->P.p != NULL) && pNext(strat->P.p) != strat->tail))
    strat->P.ShallowCopyDelete(new_tailRing, p_shallow_copy_delete);

  if ((L != NULL) && (L->tailRing != new_tailRing))
  {
    if (L->i_r < 0)
      L->ShallowCopyDelete(new_tailRing, p_shallow_copy_delete);
    else
    {
      assume(L->i_r < strat->T.size());
      TObject* t_l = strat->R[L->i_r];
      assume(t_l != NULL);
      L->tailRing = new_tailRing;
      L->p = t_l->p;
      L->t_p = t_l->t_p;
      L->max_exp = t_l->max_exp;
    }
  }

  if ((T != NULL) && (T->tailRing != new_tailRing && T->i_r < 0))
    T->ShallowCopyDelete(new_tailRing, new_tailBin, p_shallow_copy_delete);

  omMergeStickyBinIntoBin(strat->tailBin, strat->tailRing->PolyBin);
  if (strat->tailRing != currRing)
    rKillModifiedRing(strat->tailRing);

  strat->tailRing = new_tailRing;
  strat->tailBin = new_tailBin;
  strat->p_shallow_copy_delete
    = pGetShallowCopyDeleteProc(currRing, new_tailRing);

  if (strat->kNoether != NULL)
  {
    if (strat->t_kNoether != NULL)
      p_LmFree(strat->t_kNoether, strat->tailRing);
    strat->t_kNoether=k_LmInit_currRing_2_tailRing(strat->kNoether, new_tailRing);
  }

  #ifdef KDEBUG
  kTest_TS(strat);
  #endif
  if (TEST_OPT_PROT)
    PrintS("]");
  return TRUE;
}

void kStratInitChangeTailRing(kStrategy strat)
{
  unsigned long l = 0;
  int i;
  long e;

  assume(strat->tailRing == currRing);

  for (auto& Lp: strat->L)
  {
    l = p_GetMaxExpL(Lp.p, currRing, l);
  }
  for (i=0; i < strat->T.size(); i++)
  {
    // Hmm ... this we could do in one Step
    l = p_GetMaxExpL(strat->T[i].p, currRing, l);
  }
  if (rField_is_Ring(currRing))
  {
    l *= 2;
  }
  e = p_GetMaxExp(l, currRing);
  if (e <= 1) e = 2;
  if (rIsLPRing(currRing)) e = 1;

  kStratChangeTailRing(strat, NULL, NULL, e);
}

ring sbaRing (kStrategy strat, const ring r, BOOLEAN /*complete*/, int /*sgn*/)
{
  int n = rBlocks(r); // Including trailing zero!
  // if sbaOrder == 1 => use (C,monomial order from r)
  if (strat->sbaOrder == 1)
  {
    if (r->order[0] == ringorder_C || r->order[0] == ringorder_c)
    {
      return r;
    }
    ring res = rCopy0(r, TRUE, FALSE);
    res->order  = (rRingOrder_t *)omAlloc0((n+1)*sizeof(rRingOrder_t));
    res->block0 = (int *)omAlloc0((n+1)*sizeof(int));
    res->block1 = (int *)omAlloc0((n+1)*sizeof(int));
    int **wvhdl = (int **)omAlloc0((n+1)*sizeof(int*));
    res->wvhdl  = wvhdl;
    for (int i=1; i<n; i++)
    {
      res->order[i]   = r->order[i-1];
      res->block0[i]  = r->block0[i-1];
      res->block1[i]  = r->block1[i-1];
      res->wvhdl[i]   = r->wvhdl[i-1];
    }

    // new 1st block
    res->order[0]   = ringorder_C; // Prefix
    // removes useless secondary component order if defined in old ring
    for (int i=rBlocks(res); i>0; --i)
    {
      if (res->order[i] == ringorder_C || res->order[i] == ringorder_c)
      {
        res->order[i] = (rRingOrder_t)0;
      }
    }
    rComplete(res, 1);
#ifdef HAVE_PLURAL
    if (rIsPluralRing(r))
    {
      if ( nc_rComplete(r, res, false) ) // no qideal!
      {
#ifndef SING_NDEBUG
        WarnS("error in nc_rComplete");
#endif
        // cleanup?

        //      rDelete(res);
        //      return r;

        // just go on..
      }
    }
#endif
    strat->tailRing = res;
    return (res);
  }
  // if sbaOrder == 3 => degree - position - ring order
  if (strat->sbaOrder == 3)
  {
    ring res = rCopy0(r, TRUE, FALSE);
    res->order  = (rRingOrder_t*)omAlloc0((n+2)*sizeof(rRingOrder_t));
    res->block0 = (int *)omAlloc0((n+2)*sizeof(int));
    res->block1 = (int *)omAlloc0((n+2)*sizeof(int));
    int **wvhdl = (int **)omAlloc0((n+2)*sizeof(int*));
    res->wvhdl  = wvhdl;
    for (int i=2; i<n+2; i++)
    {
      res->order[i]   = r->order[i-2];
      res->block0[i]  = r->block0[i-2];
      res->block1[i]  = r->block1[i-2];
      res->wvhdl[i]   = r->wvhdl[i-2];
    }

    // new 1st block
    res->order[0]   = ringorder_a; // Prefix
    res->block0[0]  = 1;
    res->wvhdl[0]   = (int *)omAlloc(res->N*sizeof(int));
    for (int i=0; i<res->N; ++i)
      res->wvhdl[0][i]  = 1;
    res->block1[0]  = si_min(res->N, rVar(res));
    // new 2nd block
    res->order[1]   = ringorder_C; // Prefix
    res->wvhdl[1]   = NULL;
    // removes useless secondary component order if defined in old ring
    for (int i=rBlocks(res); i>1; --i)
    {
      if (res->order[i] == ringorder_C || res->order[i] == ringorder_c)
      {
        res->order[i] = (rRingOrder_t)0;
      }
    }
    rComplete(res, 1);
#ifdef HAVE_PLURAL
    if (rIsPluralRing(r))
    {
      if ( nc_rComplete(r, res, false) ) // no qideal!
      {
#ifndef SING_NDEBUG
        WarnS("error in nc_rComplete");
#endif
        // cleanup?

        //      rDelete(res);
        //      return r;

        // just go on..
      }
    }
#endif
    strat->tailRing = res;
    return (res);
  }

  // not sbaOrder == 1 => use Schreyer order
  // this is done by a trick when initializing the signatures
  // in initSLSba():
  // Instead of using the signature 1e_i for F->m[i], we start
  // with the signature LM(F->m[i])e_i for F->m[i]. Doing this we get a
  // Schreyer order w.r.t. the underlying monomial order.
  // => we do not need to change the underlying polynomial ring at all!

  // UPDATE/NOTE/TODO: use induced Schreyer ordering 'IS'!!!!????

  /*
  else
  {
    ring res = rCopy0(r, FALSE, FALSE);
    // Create 2 more blocks for prefix/suffix:
    res->order=(int *)omAlloc0((n+2)*sizeof(int)); // 0  ..  n+1
    res->block0=(int *)omAlloc0((n+2)*sizeof(int));
    res->block1=(int *)omAlloc0((n+2)*sizeof(int));
    int ** wvhdl =(int **)omAlloc0((n+2)*sizeof(int**));

    // Encapsulate all existing blocks between induced Schreyer ordering markers: prefix and suffix!
    // Note that prefix and suffix have the same ringorder marker and only differ in block[] parameters!

    // new 1st block
    int j = 0;
    res->order[j] = ringorder_IS; // Prefix
    res->block0[j] = res->block1[j] = 0;
    // wvhdl[j] = NULL;
    j++;

    for(int i = 0; (i < n) && (r->order[i] != 0); i++, j++) // i = [0 .. n-1] <- non-zero old blocks
    {
      res->order [j] = r->order [i];
      res->block0[j] = r->block0[i];
      res->block1[j] = r->block1[i];

      if (r->wvhdl[i] != NULL)
      {
        wvhdl[j] = (int*) omMemDup(r->wvhdl[i]);
      } // else wvhdl[j] = NULL;
    }

    // new last block
    res->order [j] = ringorder_IS; // Suffix
    res->block0[j] = res->block1[j] = sgn; // Sign of v[o]: 1 for C, -1 for c
    // wvhdl[j] = NULL;
    j++;

    // res->order [j] = 0; // The End!
    res->wvhdl = wvhdl;

    // j == the last zero block now!
    assume(j == (n+1));
    assume(res->order[0]==ringorder_IS);
    assume(res->order[j-1]==ringorder_IS);
    assume(res->order[j]==0);

    if (complete)
    {
      rComplete(res, 1);

#ifdef HAVE_PLURAL
      if (rIsPluralRing(r))
      {
        if ( nc_rComplete(r, res, false) ) // no qideal!
        {
        }
      }
      assume(rIsPluralRing(r) == rIsPluralRing(res));
#endif


#ifdef HAVE_PLURAL
      ring old_ring = r;

#endif

      if (r->qideal!=NULL)
      {
        res->qideal= idrCopyR_NoSort(r->qideal, r, res);

        assume(idRankFreeModule(res->qideal, res) == 0);

#ifdef HAVE_PLURAL
        if( rIsPluralRing(res) )
          if( nc_SetupQuotient(res, r, true) )
          {
            //          WarnS("error in nc_SetupQuotient"); // cleanup?      rDelete(res);       return r;  // just go on...?
          }

#endif
        assume(idRankFreeModule(res->qideal, res) == 0);
      }

#ifdef HAVE_PLURAL
      assume((res->qideal==NULL) == (old_ring->qideal==NULL));
      assume(rIsPluralRing(res) == rIsPluralRing(old_ring));
      assume(rIsSCA(res) == rIsSCA(old_ring));
      assume(ncRingType(res) == ncRingType(old_ring));
#endif
    }
    strat->tailRing = res;
    return res;
  }
  */

  assume(FALSE);
  return(NULL);
}

skStrategy::skStrategy()
{
  L.key_comp().strat = this;
  B.key_comp().strat = this;
  strat_nr++;
  nr=strat_nr;
  tailRing = currRing;
  P.tailRing = currRing;
  T.setsize(0);
  S.setsize(0);
#ifdef HAVE_LM_BIN
  lmBin = omGetStickyBinOfBin(currRing->PolyBin);
#endif
#ifdef HAVE_TAIL_BIN
  tailBin = omGetStickyBinOfBin(tailRing->PolyBin);
#endif
  pOrigFDeg = currRing->pFDeg;
  pOrigLDeg = currRing->pLDeg;
}


skStrategy::~skStrategy()
{
  if(kNoether!=NULL) pLmFree(&kNoether);
  if (lmBin != NULL)
    omMergeStickyBinIntoBin(lmBin, currRing->PolyBin);
  if (tailBin != NULL)// && !rField_is_Ring(currRing))
    omMergeStickyBinIntoBin(tailBin,
                            ((tailRing != NULL) ? tailRing->PolyBin:
                             currRing->PolyBin));
  if (t_kNoether != NULL)
    p_LmFree(t_kNoether, tailRing);

  if (currRing != tailRing)
    rKillModifiedRing(tailRing);
  pRestoreDegProcs(currRing,pOrigFDeg, pOrigLDeg);
}

#if 0
Timings for the different possibilities of posInT:
            T15           EDL         DL          EL            L         1-2-3
Gonnet      43.26       42.30       38.34       41.98       38.40      100.04
Hairer_2_1   1.11        1.15        1.04        1.22        1.08        4.7
Twomat3      1.62        1.69        1.70        1.65        1.54       11.32
ahml         4.48        4.03        4.03        4.38        4.96       26.50
c7          15.02       13.98       15.16       13.24       17.31       47.89
c8         505.09      407.46      852.76      413.21      499.19        n/a
f855        12.65        9.27       14.97        8.78       14.23       33.12
gametwo6    11.47       11.35       14.57       11.20       12.02       35.07
gerhard_3    2.73        2.83        2.93        2.64        3.12        6.24
ilias13     22.89       22.46       24.62       20.60       23.34       53.86
noon8       40.68       37.02       37.99       36.82       35.59      877.16
rcyclic_19  48.22       42.29       43.99       45.35       51.51      204.29
rkat9       82.37       79.46       77.20       77.63       82.54      267.92
schwarz_11  16.46       16.81       16.76       16.81       16.72       35.56
test016     16.39       14.17       14.40       13.50       14.26       34.07
test017     34.70       36.01       33.16       35.48       32.75       71.45
test042     10.76       10.99       10.27       11.57       10.45       23.04
test058      6.78        6.75        6.51        6.95        6.22        9.47
test066     10.71       10.94       10.76       10.61       10.56       19.06
test073     10.75       11.11       10.17       10.79        8.63       58.10
test086     12.23       11.81       12.88       12.24       13.37       66.68
test103      5.05        4.80        5.47        4.64        4.89       11.90
test154     12.96       11.64       13.51       12.46       14.61       36.35
test162     65.27       64.01       67.35       59.79       67.54      196.46
test164      7.50        6.50        7.68        6.70        7.96       17.13
virasoro     3.39        3.50        3.35        3.47        3.70        7.66
#endif


//#ifdef HAVE_MORE_POS_IN_T
#if 1
// determines the position based on: 1.) Ecart 2.) FDeg 3.) pLength
int posInT_EcartFDegpLength(const BlockArray<TObject> &set,const int length,LObject &p)
{

  if (length==-1) return 0;

  int o = p.ecart;
  int op=p.GetpFDeg();
  int ol = p.GetpLength();

  if (set[length].ecart < o)
    return length+1;
  if (set[length].ecart == o)
  {
     int oo=set[length].GetpFDeg();
     if ((oo < op) || ((oo==op) && (set[length].length < ol)))
       return length+1;
  }

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      if (set[an].ecart > o)
        return an;
      if (set[an].ecart == o)
      {
         int oo=set[an].GetpFDeg();
         if((oo > op)
         || ((oo==op) && (set[an].pLength > ol)))
           return an;
      }
      return en;
    }
    i=(an+en) / 2;
    if (set[i].ecart > o)
      en=i;
    else if (set[i].ecart == o)
    {
       int oo=set[i].GetpFDeg();
       if ((oo > op)
       || ((oo == op) && (set[i].pLength > ol)))
         en=i;
       else
        an=i;
    }
    else
      an=i;
  }
}

// determines the position based on: 1.) FDeg 2.) pLength
int posInT_FDegpLength(const BlockArray<TObject> &set,const int length,LObject &p)
{

  if (length==-1) return 0;

  int op=p.GetpFDeg();
  int ol = p.GetpLength();

  int oo=set[length].GetpFDeg();
  if ((oo < op) || ((oo==op) && (set[length].length < ol)))
    return length+1;

  int i;
  int an = 0;
  int en= length;
  loop
  {
    if (an >= en-1)
    {
      int oo=set[an].GetpFDeg();
      if((oo > op)
         || ((oo==op) && (set[an].pLength > ol)))
        return an;
      return en;
    }
    i=(an+en) / 2;
    int oo=set[i].GetpFDeg();
    if ((oo > op)
        || ((oo == op) && (set[i].pLength > ol)))
      en=i;
    else
      an=i;
  }
}


// determines the position based on: 1.) pLength
int posInT_pLength(const BlockArray<TObject> &set,const int length,LObject &p)
{
  int ol = p.GetpLength();
  if (length==-1)
    return 0;
  if (set[length].length<p.length)
    return length+1;

  int i;
  int an = 0;
  int en= length;

  loop
  {
    if (an >= en-1)
    {
      if (set[an].pLength>ol) return an;
      return en;
    }
    i=(an+en) / 2;
    if (set[i].pLength>ol) en=i;
    else                        an=i;
  }
}
#endif

// ../Singular/misc.cc:
extern char *  showOption();

void kDebugPrint(kStrategy strat)
{
  printf("red: ");
    if (strat->red==redFirst) printf("redFirst\n");
    else if (strat->red==redHoney) printf("redHoney\n");
    else if (strat->red==redEcart) printf("redEcart\n");
    else if (strat->red==redHomog) printf("redHomog\n");
    else if (strat->red==redLazy) printf("redLazy\n");
    else if (strat->red==redLiftstd) printf("redLiftstd\n");
    else  printf("%p\n",(void*)strat->red);
  printf("posInT: ");
    if (strat->posInT==posInT0) printf("posInT0\n");
    else if (strat->posInT==posInT1) printf("posInT1\n");
    else if (strat->posInT==posInT11) printf("posInT11\n");
    else if (strat->posInT==posInT110) printf("posInT110\n");
    else if (strat->posInT==posInT13) printf("posInT13\n");
    else if (strat->posInT==posInT15) printf("posInT15\n");
    else if (strat->posInT==posInT17) printf("posInT17\n");
    else if (strat->posInT==posInT17_c) printf("posInT17_c\n");
    else if (strat->posInT==posInT19) printf("posInT19\n");
    else if (strat->posInT==posInT2) printf("posInT2\n");
    else if (strat->posInT==posInT11Ring) printf("posInT11Ring\n");
    else if (strat->posInT==posInT110Ring) printf("posInT110Ring\n");
    else if (strat->posInT==posInT15Ring) printf("posInT15Ring\n");
    else if (strat->posInT==posInT17Ring) printf("posInT17Ring\n");
    else if (strat->posInT==posInT17_cRing) printf("posInT17_cRing\n");
#ifdef HAVE_MORE_POS_IN_T
    else if (strat->posInT==posInT_EcartFDegpLength) printf("posInT_EcartFDegpLength\n");
    else if (strat->posInT==posInT_FDegpLength) printf("posInT_FDegpLength\n");
    else if (strat->posInT==posInT_pLength) printf("posInT_pLength\n");
#endif
    else if (strat->posInT==posInT_EcartpLength) printf("posInT_EcartpLength\n");
    else  printf("%p\n",(void*)strat->posInT);
  printf("compareL: ");
    if (strat->compareL==compareL0) printf("compareL0\n");
    else if (strat->compareL==compareL10) printf("compareL10\n");
    else if (strat->compareL==compareL11) printf("compareL11\n");
    else if (strat->compareL==compareL110) printf("compareL110\n");
    else if (strat->compareL==compareL13) printf("compareL13\n");
    else if (strat->compareL==compareL15) printf("compareL15\n");
    else if (strat->compareL==compareL17) printf("compareL17\n");
    else if (strat->compareL==compareL17_c) printf("compareL17_c\n");
    else if (strat->compareL==compareL0) printf("compareL0Ring\n");
    else if (strat->compareL==compareL11Ring) printf("compareL11Ring\n");
    else if (strat->compareL==compareL11Ringls) printf("compareL11Ringls\n");
    else if (strat->compareL==compareL110Ring) printf("compareL110Ring\n");
    else if (strat->compareL==compareL15Ring) printf("compareL15Ring\n");
    else if (strat->compareL==compareL17Ring) printf("compareL17Ring\n");
    else if (strat->compareL==compareL17_cRing) printf("compareL17_cRing\n");
    else if (strat->compareL==compareLSpecial) printf("compareLSpecial\n");
    else  printf("%p\n",(void*)strat->compareL);
  printf("enterS: ");
    if (strat->enterS==enterSBba) printf("enterSBba\n");
    else if (strat->enterS==enterSMora) printf("enterSMora\n");
    else if (strat->enterS==enterSMoraNF) printf("enterSMoraNF\n");
    else  printf("%p\n",(void*)strat->enterS);
  printf("initEcart: ");
    if (strat->initEcart==initEcartBBA) printf("initEcartBBA\n");
    else if (strat->initEcart==initEcartNormal) printf("initEcartNormal\n");
    else  printf("%p\n",(void*)strat->initEcart);
  printf("initEcartPair: ");
    if (strat->initEcartPair==initEcartPairBba) printf("initEcartPairBba\n");
    else if (strat->initEcartPair==initEcartPairMora) printf("initEcartPairMora\n");
    else  printf("%p\n",(void*)strat->initEcartPair);
  printf("homog=%d, LazyDegree=%d, LazyPass=%d, ak=%d,\n",
    strat->homog, strat->LazyDegree,strat->LazyPass, strat->ak);
  printf("honey=%d, sugarCrit=%d, Gebauer=%d, noTailReduction=%d, use_buckets=%d\n",
    strat->honey,strat->sugarCrit,strat->Gebauer,strat->noTailReduction,strat->use_buckets);
  printf("chainCrit: ");
    if (strat->chainCrit==chainCritNormal) printf("chainCritNormal\n");
    else if (strat->chainCrit==chainCritOpt_1) printf("chainCritOpt_1\n");
    else  printf("%p\n",(void*)strat->chainCrit);
  printf("compareLDependsOnLength=%d\n",
         strat->compareLDependsOnLength);
  printf("%s\n",showOption());
  printf("LDeg: ");
    if (currRing->pLDeg==pLDeg0) printf("pLDeg0");
    else if (currRing->pLDeg==pLDeg0c) printf("pLDeg0c");
    else if (currRing->pLDeg==pLDegb) printf("pLDegb");
    else if (currRing->pLDeg==pLDeg1) printf("pLDeg1");
    else if (currRing->pLDeg==pLDeg1c) printf("pLDeg1c");
    else if (currRing->pLDeg==pLDeg1_Deg) printf("pLDeg1_Deg");
    else if (currRing->pLDeg==pLDeg1c_Deg) printf("pLDeg1c_Deg");
    else if (currRing->pLDeg==pLDeg1_Totaldegree) printf("pLDeg1_Totaldegree");
    else if (currRing->pLDeg==pLDeg1c_Totaldegree) printf("pLDeg1c_Totaldegree");
    else if (currRing->pLDeg==pLDeg1_WFirstTotalDegree) printf("pLDeg1_WFirstTotalDegree");
    else if (currRing->pLDeg==pLDeg1c_WFirstTotalDegree) printf("pLDeg1c_WFirstTotalDegree");
    else if (currRing->pLDeg==maxdegreeWecart) printf("maxdegreeWecart");
    else printf("? (%lx)", (long)currRing->pLDeg);
    printf(" / ");
    if (strat->tailRing->pLDeg==pLDeg0) printf("pLDeg0");
    else if (strat->tailRing->pLDeg==pLDeg0c) printf("pLDeg0c");
    else if (strat->tailRing->pLDeg==pLDegb) printf("pLDegb");
    else if (strat->tailRing->pLDeg==pLDeg1) printf("pLDeg1");
    else if (strat->tailRing->pLDeg==pLDeg1c) printf("pLDeg1c");
    else if (strat->tailRing->pLDeg==pLDeg1_Deg) printf("pLDeg1_Deg");
    else if (strat->tailRing->pLDeg==pLDeg1c_Deg) printf("pLDeg1c_Deg");
    else if (strat->tailRing->pLDeg==pLDeg1_Totaldegree) printf("pLDeg1_Totaldegree");
    else if (strat->tailRing->pLDeg==pLDeg1c_Totaldegree) printf("pLDeg1c_Totaldegree");
    else if (strat->tailRing->pLDeg==pLDeg1_WFirstTotalDegree) printf("pLDeg1_WFirstTotalDegree");
    else if (strat->tailRing->pLDeg==pLDeg1c_WFirstTotalDegree) printf("pLDeg1c_WFirstTotalDegree");
    else if (strat->tailRing->pLDeg==maxdegreeWecart) printf("maxdegreeWecart");
    else printf("? (%lx)", (long)strat->tailRing->pLDeg);
    printf("\n");
  printf("currRing->pFDeg: ");
    if (currRing->pFDeg==p_Totaldegree) printf("p_Totaldegree");
    else if (currRing->pFDeg==p_WFirstTotalDegree) printf("pWFirstTotalDegree");
    else if (currRing->pFDeg==p_Deg) printf("p_Deg");
    else if (currRing->pFDeg==kHomModDeg) printf("kHomModDeg");
    else if (currRing->pFDeg==kModDeg) printf("kModDeg");
    else if (currRing->pFDeg==totaldegreeWecart) printf("totaldegreeWecart");
    else if (currRing->pFDeg==p_WTotaldegree) printf("p_WTotaldegree");
    else printf("? (%lx)", (long)currRing->pFDeg);
    printf("\n");
    printf(" syzring:%d, syzComp(strat):%d limit:%d\n",rIsSyzIndexRing(currRing),strat->syzComp,rGetCurrSyzLimit(currRing));
    if(TEST_OPT_DEGBOUND)
      printf(" degBound: %d\n", Kstd1_deg);

    if( ecartWeights != NULL )
    {
       printf("ecartWeights: ");
       for (int i = rVar(currRing); i > 0; i--)
         printf("%hd ", ecartWeights[i]);
       printf("\n");
       assume( TEST_OPT_WEIGHTM );
    }

#ifndef SING_NDEBUG
    rDebugPrint(currRing);
#endif
}

//LObject pCopyp2L(poly p, kStrategy strat)
//{
    /* creates LObject from the poly in currRing */
  /* actually put p into L.p and make L.t_p=NULL : does not work */

//}

/*2
* put the  lcm(q,p)  into the set B, q is the shift of some s[i]
*/
#ifdef HAVE_SHIFTBBA
static BOOLEAN enterOneStrongPolyShift (poly q, poly p, int /*ecart*/, int /*isFromQ*/, kStrategy strat, int atR, int /*ecartq*/, int /*qisFromQ*/, int shiftcount, sBasisSet::const_iterator ifromS)
{
  number d, s, t;
  /* assume(atR >= 0); */
  assume(ifromS == strat->S.end() || ifromS.index() < strat->S.size());
  assume(rField_is_Ring(currRing));
  poly m1, m2, gcd;
  //printf("\n--------------------------------\n");
  //pWrite(p);pWrite(si);
  d = n_ExtGcd(pGetCoeff(p), pGetCoeff(q), &s, &t, currRing->cf);

  if (nIsZero(s) || nIsZero(t))  // evtl. durch divBy tests ersetzen
  {
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
    return FALSE;
  }

  assume(pIsInV(p));

  k_GetStrongLeadTerms(p, q, currRing, m1, m2, gcd, strat->tailRing);

  /* the V criterion */
  if (!pmIsInV(gcd))
  {
    strat->cv++;
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
    pLmFree(gcd);
    return FALSE;
  }

  // disabled for Letterplace because it is not so easy to check
  /* if (!rHasLocalOrMixedOrdering(currRing)) { */
  /*   unsigned long sev = pGetShortExpVector(gcd); */

  /*   for (int j = 0; j < strat->S.size()-1; j++) { */
  /*     if (j == i) */
  /*       continue; */

  /*     if (n_DivBy(d, pGetCoeff(strat->S.iterator_at(j)->p), currRing->cf) && */
  /*         !(strat->S.iterator_at(j)->sev & ~sev) && */
  /*         p_LmDivisibleBy(strat->S.iterator_at(j)->p, gcd, currRing)) { */
  /*       nDelete(&d); */
  /*       nDelete(&s); */
  /*       nDelete(&t); */
  /*       return FALSE; */
  /*     } */
  /*   } */
  /* } */

  poly m12, m22;
  assume(p_mFirstVblock(p, currRing) <= 1 || p_mFirstVblock(q, currRing) <= 1);
  k_SplitFrame(m1, m12, si_max(p_mFirstVblock(p, currRing), 1), currRing);
  k_SplitFrame(m2, m22, si_max(p_mFirstVblock(q, currRing), 1), currRing);
  // manually free the coeffs, because pSetCoeff0 is used in the next step
  n_Delete(&(m1->coef), currRing->cf);
  n_Delete(&(m2->coef), currRing->cf);

  //p_Test(m1,strat->tailRing);
  //p_Test(m2,strat->tailRing);
  /*if(!enterTstrong)
  {
    while (! kCheckStrongCreation(atR, m1, i, m2, strat) )
    {
      memset(&(strat->P), 0, sizeof(strat->P));
      kStratChangeTailRing(strat);
      strat->P = *(strat->R[atR]);
      p_LmFree(m1, strat->tailRing);
      p_LmFree(m2, strat->tailRing);
      p_LmFree(gcd, currRing);
      k_GetStrongLeadTerms(p, si, currRing, m1, m2, gcd, strat->tailRing);
    }
  }*/
  pSetCoeff0(m1, s);
  pSetCoeff0(m2, t);
  pSetCoeff0(gcd, d);
  p_Test(m1,strat->tailRing);
  p_Test(m2,strat->tailRing);
  p_Test(m12,strat->tailRing);
  p_Test(m22,strat->tailRing);
  assume(pmIsInV(m1));
  assume(pmIsInV(m2));
  assume(pmIsInV(m12));
  assume(pmIsInV(m22));
  //printf("\n===================================\n");
  //pWrite(m1);pWrite(m2);pWrite(gcd);
#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    // Print("t = %d; s = %d; d = %d\n", nInt(t), nInt(s), nInt(d));
    PrintS("m1 = ");
    p_wrp(m1, strat->tailRing);
    PrintS("m12 = ");
    p_wrp(m12, strat->tailRing);
    PrintS(" ; m2 = ");
    p_wrp(m2, strat->tailRing);
    PrintS(" ; m22 = ");
    p_wrp(m22, strat->tailRing);
    PrintS(" ; gcd = ");
    wrp(gcd);
    PrintS("\n--- create strong gcd poly: ");
    PrintS("\n p: ");
    wrp(p);
    Print("\n q (strat->S[%d].p): ", ifromS == strat->S.end() ? -1 : ifromS.index());
    wrp(q);
    PrintS(" ---> ");
  }
#endif

  pNext(gcd) = p_Add_q(pp_Mult_mm(pp_mm_Mult(pNext(p), m1, strat->tailRing), m12, strat->tailRing), pp_Mult_mm(pp_mm_Mult(pNext(q), m2, strat->tailRing), m22, strat->tailRing), strat->tailRing);
  p_LmDelete(m1, strat->tailRing);
  p_LmDelete(m2, strat->tailRing);
  p_LmDelete(m12, strat->tailRing);
  p_LmDelete(m22, strat->tailRing);

  assume(pIsInV(gcd));

#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    wrp(gcd);
    PrintLn();
  }
#endif

  LObject h;
  h.p = gcd;
  h.tailRing = strat->tailRing;
  strat->initEcart(&h);
  h.sev = pGetShortExpVector(h.p);
  h.i_r1 = -1;h.i_r2 = -1;
  if (currRing!=strat->tailRing)
    h.t_p = k_LmInit_currRing_2_tailRing(h.p, strat->tailRing);
#if 1
  h.p1 = p;
  h.p2 = q;
#endif
  if (atR >= 0 && shiftcount == 0 && ifromS != strat->S.end())
  {
    h.i_r2 = kFindInT(h.p1, strat);
    h.i_r1 = atR;
  }
  else
  {
    h.i_r1 = -1;
    h.i_r2 = -1;
  }

  assume(pIsInV(h.p));
  assume(pIsInV(h.p1));

  strat->L.push(h);
  return TRUE;
}
#endif


/*2
* put the pair (q,p)  into the set B, ecart=ecart(p), q is the shift of some s[i] (ring case)
*/
#ifdef HAVE_SHIFTBBA
static void enterOnePairRingShift (poly q, poly p, int /*ecart*/, int isFromQ, kStrategy strat, int atR, int /*ecartq*/, int qisFromQ, int shiftcount, sBasisSet::const_iterator ifromS)
{
  /* assume(atR >= 0); */
  /* assume(i < strat->S.size()); */
  assume(p!=NULL);
  assume(rField_is_Ring(currRing));
  assume(pIsInV(p));
  #if ALL_VS_JUST
  //Over rings, if we construct the strong pair, do not add the spair
  if(rField_is_Ring(currRing))
  {
    number s,t,d;
    d = n_ExtGcd(pGetCoeff(p), pGetCoeff(q, &s, &t, currRing->cf);

    if (!nIsZero(s) && !nIsZero(t))  // evtl. durch divBy tests ersetzen
    {
      nDelete(&d);
      nDelete(&s);
      nDelete(&t);
      return;
    }
    nDelete(&d);
    nDelete(&s);
    nDelete(&t);
  }
  #endif
  int      compare,compareCoeff;
  LObject  h;

#ifdef KDEBUG
  h.ecart=0; h.length=0;
#endif
  /*- computes the lcm(s[i],p) -*/
  if(pHasNotCFRing(p,q))
  {
      strat->cp++;
      return;
  }
  h.lcm = p_Lcm(p,q,currRing);
  pSetCoeff0(h.lcm, n_Lcm(pGetCoeff(p), pGetCoeff(q), currRing->cf));
  if (nIsZero(pGetCoeff(h.lcm)))
  {
      strat->cp++;
      pLmDelete(h.lcm);
      return;
  }

  /* the V criterion */
  if (!pmIsInV(h.lcm))
  {
    strat->cv++;
    pLmDelete(h.lcm);
    return;
  }
  // basic chain criterion
  /*
  *the set B collects the pairs of type (S[j],p)
  *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p) != lcm(r,p)
  *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
  *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
  */

  for(auto jt = strat_B(strat).ubegin(); jt != strat_B(strat).uend(); )
  {
    bool j_deleted = false;
    compare=pDivCompRing(jt->lcm,h.lcm);
    compareCoeff = n_DivComp(pGetCoeff(jt->lcm), pGetCoeff(h.lcm), currRing->cf);
    if(compare == pDivComp_EQUAL)
    {
      //They have the same LM
      if(compareCoeff == pDivComp_LESS)
      {
        if ((!strat->hasFromQ) || (isFromQ==0) || (qisFromQ==0))
        {
          strat->c3++;
          pLmDelete(h.lcm);
          return;
        }
        break;
      }
      if(compareCoeff == pDivComp_GREATER)
      {
        jt = strat_B(strat).erase(jt);
        j_deleted = true;
        strat->c3++;
      }
      if(compareCoeff == pDivComp_EQUAL)
      {
        if ((!strat->hasFromQ) || (isFromQ==0) || (qisFromQ==0))
        {
          strat->c3++;
          pLmDelete(h.lcm);
          return;
        }
        break;
      }
    }
    if(compareCoeff == compare || compareCoeff == pDivComp_EQUAL)
    {
      if(compare == pDivComp_LESS)
      {
        if ((!strat->hasFromQ) || (isFromQ==0) || (qisFromQ==0))
        {
          strat->c3++;
          pLmDelete(h.lcm);
          return;
        }
        break;
      }
      if(compare == pDivComp_GREATER)
      {
        jt = strat_B(strat).erase(jt);
        j_deleted = true;
        strat->c3++;
      }
    }
    if (! j_deleted) ++jt;
  }
  number s, t;
  poly m1, m2, gcd = NULL;
  s = pGetCoeff(q);
  t = pGetCoeff(p);
  k_GetLeadTerms(p,q,currRing,m1,m2,currRing);

  poly m12, m22;
  assume(p_mFirstVblock(p, currRing) <= 1 || p_mFirstVblock(q, currRing) <= 1);
  k_SplitFrame(m1, m12, si_max(p_mFirstVblock(p, currRing), 1), currRing);
  k_SplitFrame(m2, m22, si_max(p_mFirstVblock(q, currRing), 1), currRing);
  // manually free the coeffs, because pSetCoeff0 is used in the next step
  n_Delete(&(m1->coef), currRing->cf);
  n_Delete(&(m2->coef), currRing->cf);

  ksCheckCoeff(&s, &t, currRing->cf);
  pSetCoeff0(m1, s);
  pSetCoeff0(m2, t);
  m2 = pNeg(m2);
  p_Test(m1,strat->tailRing);
  p_Test(m2,strat->tailRing);
  p_Test(m12,strat->tailRing);
  p_Test(m22,strat->tailRing);
  assume(pmIsInV(m1));
  assume(pmIsInV(m2));
  assume(pmIsInV(m12));
  assume(pmIsInV(m22));
  poly pm1 = pp_Mult_mm(pp_mm_Mult(pNext(p), m1, strat->tailRing), m12, strat->tailRing);
  poly sim2 = pp_Mult_mm(pp_mm_Mult(pNext(q), m2, strat->tailRing), m22, strat->tailRing);
  assume(pIsInV(pm1));
  assume(pIsInV(sim2));
  p_LmDelete(m1, currRing);
  p_LmDelete(m2, currRing);
  p_LmDelete(m12, currRing);
  p_LmDelete(m22, currRing);
  if(sim2 == NULL)
  {
    if(pm1 == NULL)
    {
      if(h.lcm != NULL)
      {
        pLmDelete(h.lcm);
        h.lcm=NULL;
      }
      h.Clear();
      /* TEMPORARILY DISABLED FOR SHIFTS because there is no i*/
      /* strat->pairtest[i] = TRUE; */
      return;
    }
    else
    {
      gcd = pm1;
      pm1 = NULL;
    }
  }
  else
  {
    if((pGetComp(q) == 0) && (0 != pGetComp(p)))
    {
      p_SetCompP(sim2, pGetComp(p), strat->tailRing);
      pSetmComp(sim2);
    }
    //p_Write(pm1,strat->tailRing);p_Write(sim2,strat->tailRing);
    gcd = p_Add_q(pm1, sim2, strat->tailRing);
  }
  p_Test(gcd, strat->tailRing);
  assume(pIsInV(gcd));
#ifdef KDEBUG
  if (TEST_OPT_DEBUG)
  {
    wrp(gcd);
    PrintLn();
  }
#endif
  h.p = gcd;
  h.i_r = -1;
  if(h.p == NULL)
  {
    /* TEMPORARILY DISABLED FOR SHIFTS because there is no i*/
    /* strat->pairtest[i] = TRUE; */
    return;
  }
  h.tailRing = strat->tailRing;
  //h.pCleardenom();
  //pSetm(h.p);
  h.i_r1 = -1;h.i_r2 = -1;
  strat->initEcart(&h);
  #if 1
  h.p1 = p;
  h.p2 = q;
  #endif
  #if 1
  /* TEMPORARILY DISABLED FOR SHIFTS because there's no i*/
  /* at the beginning we DO NOT set atR = -1 ANYMORE*/
  if (atR >= 0 && shiftcount == 0 && ifromS != strat->S.end())
  {
    h.i_r2 = kFindInT(h.p1, strat); //strat->S[i].s_2_r;
    h.i_r1 = atR;
  }
  else
  {
    /* END _ TEMPORARILY DISABLED FOR SHIFTS */
    h.i_r1 = -1;
    h.i_r2 = -1;
  }
  #endif
  h.sev = pGetShortExpVector(h.p);
  if (currRing!=strat->tailRing)
    h.t_p = k_LmInit_currRing_2_tailRing(h.p, strat->tailRing);

  assume(pIsInV(h.p));
  assume(pIsInV(h.p1));
  assume(h.lcm != NULL);
  assume(pIsInV(h.lcm));

  strat_B(strat).push(h);
  kTest_TS(strat);
}
#endif

#ifdef HAVE_SHIFTBBA
// adds the strong pair and the normal pair for rings (aka gpoly and spoly)
static BOOLEAN enterOneStrongPolyAndEnterOnePairRingShift(poly q, poly p, int ecart, int isFromQ, kStrategy strat, int atR, int ecartq, int qisFromQ, int shiftcount, sBasisSet::const_iterator ifromS)
{
  enterOneStrongPolyShift(q, p, ecart, isFromQ, strat, atR, ecartq, qisFromQ, shiftcount, ifromS); // "gpoly"
  enterOnePairRingShift(q, p, ecart, isFromQ, strat, atR, ecartq, qisFromQ, shiftcount, ifromS); // "spoly"
  return FALSE; // TODO: delete q?
}
#endif

#ifdef HAVE_SHIFTBBA
// creates if possible (q,p), (shifts(q),p)
static BOOLEAN enterOnePairWithShifts (sBasisSet::const_iterator q_inS /*also i*/, poly q, poly p, int ecartp, int p_isFromQ, kStrategy strat, int /*atR*/, int p_lastVblock, int q_lastVblock)
{
  // note: ecart and isFromQ is for p
  assume(q_inS == strat->S.end() || q_inS->p == q); // if q is from S, q_inS should be the index of q in S
  assume(pmFirstVblock(p) == 1);
  assume(pmFirstVblock(q) == 1);
  assume(p_lastVblock == pmLastVblock(p));
  assume(q_lastVblock == pmLastVblock(q));

  // TODO: is ecartq = 0 still ok?
  int ecartq = 0; //Hans says it's ok; we're in the homog case, no ecart

  int q_isFromQ = 0;
  if (strat->hasFromQ && q_inS != strat->S.end())
    q_isFromQ = q_inS->fromQ;

  BOOLEAN (*enterPair)(poly, poly, int, int, kStrategy, int, int, int, int, sBasisSet::const_iterator);
  if (rField_is_Ring(currRing))
    enterPair = enterOneStrongPolyAndEnterOnePairRingShift;
  else
    enterPair = enterOnePairShift;

  int degbound = currRing->N/currRing->isLPring;
  int neededShift = p_lastVblock - ((pGetComp(p) > 0 || pGetComp(q) > 0) ? 0 : 1); // in the module case, product criterion does not hold
  int maxPossibleShift = degbound - q_lastVblock;
  int maxShift = si_min(neededShift, maxPossibleShift);
  int firstShift = (q == p ? 1 : 0); // do not add (q,p) if q=p
  BOOLEAN delete_pair=TRUE;
  for (int j = firstShift; j <= maxShift; j++)
  {
    poly qq = pLPCopyAndShiftLM(q, j);
    if (enterPair(qq, p, ecartp, p_isFromQ, strat, -1, ecartq, q_isFromQ, j, q_inS))
    {
      if (j>0) pLmDelete(qq);
      // delete qq, if not it does not enter the pair set
    }
    else
      delete_pair=FALSE;
  }

  if (rField_is_Ring(currRing) && p_lastVblock >= firstShift && p_lastVblock <= maxPossibleShift)
  {
    // add pairs (m*shifts(q), p) where m is a monomial and the pair has no overlap
    for (int j = p_lastVblock; j <= maxPossibleShift; j++)
    {
      ideal fillers = id_MaxIdeal(j - p_lastVblock, currRing);
      for (int k = 0; k < IDELEMS(fillers); k++)
      {
        poly qq = pLPCopyAndShiftLM(pp_mm_Mult(q, fillers->m[k], currRing), p_lastVblock);
        enterPair(qq, p, ecartp, p_isFromQ, strat, -1, ecartq, q_isFromQ, p_lastVblock, q_inS);
      }
      idDelete(&fillers);
    }
  }
  return delete_pair;
}
#endif

#ifdef HAVE_SHIFTBBA
// creates (q,p), use it when q is already shifted
// return TRUE, if (q,p) is discarded
static BOOLEAN enterOnePairWithoutShifts (sBasisSet::const_iterator p_inS /*also i*/, poly q, poly p, int ecartq, int q_isFromQ, kStrategy strat, int /*atR*/, int p_lastVblock, int q_shift)
{
  // note: ecart and isFromQ is for p
  assume(p_inS == strat->S.end() || p_inS->p == p); // if p is from S, p_inS should be the index of p in S
  assume(pmFirstVblock(p) == 1);
  assume(p_lastVblock == pmLastVblock(p));
  assume(q_shift == pmFirstVblock(q) - 1);

  // TODO: is ecartp = 0 still ok?
  int ecartp = 0; //Hans says it's ok; we're in the homog e:, no ecart

  int p_isFromQ = 0;
  if (strat->hasFromQ && p_inS != strat->S.end())
    p_isFromQ = p_inS->fromQ;

  if (rField_is_Ring(currRing))
  {
    assume(q_shift <= p_lastVblock); // we allow the special case where there is no overlap
    return enterOneStrongPolyAndEnterOnePairRingShift(q, p, ecartp, p_isFromQ, strat, -1, ecartq, q_isFromQ, q_shift, strat->S.end());
  }
  else
  {
    assume(q_shift <= p_lastVblock - ((pGetComp(q) > 0 || pGetComp(p) > 0) ? 0 : 1)); // there should be an overlap (in the module case epsilon overlap is also allowed)
    return enterOnePairShift(q, p, ecartp, p_isFromQ, strat, -1, ecartq, q_isFromQ, q_shift, strat->S.end());
  }
}
#endif


#ifdef KDEBUG
// enable to print which pairs are considered or discarded and why
/* #define CRITERION_DEBUG */
#endif
/*2
* put the pair (q,p)  into the set B, ecart=ecart(p), q is the shift of some s[i]
* return TRUE, if (q,p) does not enter B
*/
#ifdef HAVE_SHIFTBBA
BOOLEAN enterOnePairShift (poly q, poly p, int ecart, int isFromQ, kStrategy strat, int atR, int ecartq, int qisFromQ, int shiftcount, sBasisSet::const_iterator ifromS)
{
#ifdef CRITERION_DEBUG
  if (TEST_OPT_DEBUG)
  {
    PrintS("Consider pair ("); wrp(q); PrintS(", "); wrp(p); PrintS(")"); PrintLn();
    // also write the LMs in separate lines:
    poly lmq = pHead(q);
    poly lmp = pHead(p);
    pSetCoeff(lmq, n_Init(1, currRing->cf));
    pSetCoeff(lmp, n_Init(1, currRing->cf));
    Print("    %s\n", pString(lmq));
    Print("    %s\n", pString(lmp));
    pLmDelete(lmq);
    pLmDelete(lmp);
  }
#endif

  /* Format: q and p are like strat->P.p, so lm in CR, tail in TR */

  /* check this Formats: */
  assume(p_LmCheckIsFromRing(q,currRing));
  assume(p_CheckIsFromRing(pNext(q),strat->tailRing));
  assume(p_LmCheckIsFromRing(p,currRing));
  assume(p_CheckIsFromRing(pNext(p),strat->tailRing));

  /* poly q stays for s[i], ecartq = ecart(q), qisFromQ = applies to q */

  int qfromQ = qisFromQ;

  /* need additionally: int up_to_degree, poly V0 with the variables in (0)  or just the number lV = the length of the first block */

  int      compare;
  LObject  Lp;
  Lp.i_r = -1;

#ifdef KDEBUG
  Lp.ecart=0; Lp.length=0;
#endif
  /*- computes the lcm(s[i],p) -*/
  Lp.lcm = p_Lcm(p,q, currRing); // q is what was strat->S[i].p, so a poly in LM/TR presentation
  Lp.sev_lcm = p_GetShortExpVector(Lp.lcm, currRing);

  /* the V criterion */
  if (!pmIsInV(Lp.lcm))
  {
    strat->cv++; // counter for applying the V criterion
    pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
    if (TEST_OPT_DEBUG) PrintS("--- V crit\n");
#endif
    return TRUE;
  }

  if (strat->sugarCrit && ALLOW_PROD_CRIT(strat))
  {
    if((!((ecartq>0)&&(ecart>0)))
    && pHasNotCF(p,q))
    {
    /*
    *the product criterion has applied for (s,p),
    *i.e. lcm(s,p)=product of the leading terms of s and p.
    *Suppose (s,r) is in L and the leading term
    *of p divides lcm(s,r)
    *(==> the leading term of p divides the leading term of r)
    *but the leading term of s does not divide the leading term of r
    *(notice that this condition is automatically satisfied if r is still
    *in S), then (s,r) can be cancelled.
    *This should be done here because the
    *case lcm(s,r)=lcm(s,p) is not covered by chainCrit.
    *
    *Moreover, skipping (s,r) holds also for the noncommutative case.
    */
      strat->cp++;
      pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
      if (TEST_OPT_DEBUG) PrintS("--- prod crit\n");
#endif
      return TRUE;
    }
    else
      Lp.ecart = si_max(ecart,ecartq);
    if (strat->fromT && (ecartq>ecart))
    {
      pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
      if (TEST_OPT_DEBUG) PrintS("--- ecartq > ecart\n");
#endif
      return TRUE;
      /*the pair is (s[i],t[.]), discard it if the ecart is too big*/
    }
    /*
    *the set B collects the pairs of type (S[j],p)
    *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p)#lcm(r,p)
    *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
    *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
    */
    {
      for (auto jt = strat_B(strat).ubegin(); jt != strat_B(strat).uend(); )
      {
        compare=pLPDivComp(jt->lcm,Lp.lcm);
        if ((compare==1)
        &&(sugarDivisibleBy(jt->ecart,Lp.ecart)))
        {
          strat->c3++;
          if ((!strat->hasFromQ) || (isFromQ==0) || (qfromQ==0))
          {
            pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
            if (TEST_OPT_DEBUG)
            {
              Print("--- chain crit using B[j].lcm=%s\n", pString(jt->lcm));
            }
#endif
            return TRUE;
          }
          break;
        }
        else
        if ((compare ==-1)
        && sugarDivisibleBy(Lp.ecart,jt->ecart))
        {
#ifdef CRITERION_DEBUG
          if (TEST_OPT_DEBUG)
          {
            Print("--- chain crit using pair to remove B[j].lcm=%s\n", pString(jt->lcm));
          }
#endif
          jt = strat_B(strat).erase(jt);
          strat->c3++;
        }
        else
          ++jt;
      }
    }
  }
  else /*sugarcrit*/
  {
    if (ALLOW_PROD_CRIT(strat))
    {
      // if currRing->nc_type!=quasi (or skew)
      // TODO: enable productCrit for super commutative algebras...
      if(/*(strat->ak==0) && productCrit(p,strat->S.iterator_at(i)->p)*/
      pHasNotCF(p,q))
      {
      /*
      *the product criterion has applied for (s,p),
      *i.e. lcm(s,p)=product of the leading terms of s and p.
      *Suppose (s,r) is in L and the leading term
      *of p divides lcm(s,r)
      *(==> the leading term of p divides the leading term of r)
      *but the leading term of s does not divide the leading term of r
      *(notice that tis condition is automatically satisfied if r is still
      *in S), then (s,r) can be canceled.
      *This should be done here because the
      *case lcm(s,r)=lcm(s,p) is not covered by chainCrit.
      */
          strat->cp++;
          pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
          if (TEST_OPT_DEBUG) PrintS("--- prod crit\n");
#endif
          return TRUE;
      }
      if (strat->fromT && (ecartq>ecart))
      {
        pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
        if (TEST_OPT_DEBUG) PrintS("--- ecartq > ecart\n");
#endif
        return TRUE;
        /*the pair is (s[i],t[.]), discard it if the ecart is too big*/
      }
      /*
      *the set B collects the pairs of type (S[j],p)
      *suppose (r,p) is in B and (s,p) is the new pair and lcm(s,p)#lcm(r,p)
      *if the leading term of s divides lcm(r,p) then (r,p) will be canceled
      *if the leading term of r divides lcm(s,p) then (s,p) will not enter B
      */
      for(auto jt = strat_B(strat).ubegin(); jt != strat_B(strat).uend(); )
      {
        compare=pLPDivComp(jt->lcm,Lp.lcm);
        if (compare==1)
        {
          strat->c3++;
          if ((!strat->hasFromQ) || (isFromQ==0) || (qfromQ==0))
          {
            pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
            if (TEST_OPT_DEBUG)
            {
              Print("--- chain crit using B[j].lcm=%s\n", pString(jt->lcm));
            }
#endif
            return TRUE;
          }
          break;
        }
        else
        if (compare ==-1)
        {
#ifdef CRITERION_DEBUG
          if (TEST_OPT_DEBUG)
          {
            Print("--- chain crit using pair to remove B[j].lcm=%s\n", pString(jt->lcm));
          }
#endif
          jt = strat_B(strat).erase(jt);
          strat->c3++;
        }
        else
          ++jt;
      }
    }
  }
  /*
  *the pair (S[i],p) enters B if the spoly != 0
  */
  /*-  compute the short s-polynomial -*/
  if (strat->fromT && !TEST_OPT_INTSTRATEGY)
    pNorm(p);
  if ((q==NULL) || (p==NULL))
  {
#ifdef CRITERION_DEBUG
    if (TEST_OPT_DEBUG) PrintS("--- q == NULL || p == NULL\n");
#endif
    return FALSE;
  }
  if ((strat->hasFromQ) && (isFromQ!=0) && (qfromQ!=0))
  {
    Lp.p=NULL;
#ifdef CRITERION_DEBUG
    if (TEST_OPT_DEBUG) PrintS("--- pair is from Q\n");
#endif
  }
  else
  {
//     if ( rIsPluralRing(currRing) )
//     {
//       if(pHasNotCF(p, q))
//       {
//         if(ncRingType(currRing) == nc_lie)
//         {
//             // generalized prod-crit for lie-type
//             strat->cp++;
//             Lp.p = nc_p_Bracket_qq(pCopy(p),q, currRing);
//         }
//         else
//         if( ALLOW_PROD_CRIT(strat) )
//         {
//             // product criterion for homogeneous case in SCA
//             strat->cp++;
//             Lp.p = NULL;
//         }
//         else
//           Lp.p = nc_CreateSpoly(q,p,currRing); // ?
//       }
//       else  Lp.p = nc_CreateSpoly(q,p,currRing);
//     }
//     else
//     {

    /* ksCreateShortSpoly needs two Lobject-kind presentations */
    /* p is already in this form, so convert q */
    Lp.p = ksCreateShortSpoly(q, p, strat->tailRing);
      //  }
  }
  if (Lp.p == NULL)
  {
    /*- the case that the s-poly is 0 -*/
    // TODO: currently ifromS is only > 0 if called from enterOnePairWithShifts
    if (ifromS != strat->S.end() && ifromS.index() > 0)
    {
      record_pairtest_hit(*ifromS, strat);/*- hint for spoly(S^[i],p)=0 -*/
    }
      //if (TEST_OPT_DEBUG){Print("!");} // option teach
    /* END _ TEMPORARILY DISABLED FOR SHIFTS */
    /*hint for spoly(S[i],p) == 0 for some i,0 <= i <= sl*/
    /*
    *suppose we have (s,r),(r,p),(s,p) and spoly(s,p) == 0 and (r,p) is
    *still in B (i.e. lcm(r,p) == lcm(s,p) or the leading term of s does not
    *divide lcm(r,p)). In the last case (s,r) can be canceled if the leading
    *term of p divides the lcm(s,r)
    *(this canceling should be done here because
    *the case lcm(s,p) == lcm(s,r) is not covered in chainCrit)
    *the first case is handled in chainCrit
    */
    if (Lp.lcm!=NULL) pLmFree(Lp.lcm);
#ifdef CRITERION_DEBUG
    if (TEST_OPT_DEBUG) PrintS("--- S-poly = 0\n");
#endif
    return TRUE;
  }
  else
  {
    /*- the pair (S[i],p) enters B -*/
    /* both of them should have their LM in currRing and TAIL in tailring */
    Lp.p1 = q;  // already in the needed form
    Lp.p2 = p; // already in the needed form

    if ( !rIsPluralRing(currRing) )
      pNext(Lp.p) = strat->tail;

    /* TEMPORARILY DISABLED FOR SHIFTS because there's no i*/
    /* at the beginning we DO NOT set atR = -1 ANYMORE*/
    if ( (atR >= 0) && (shiftcount==0) && (ifromS != strat->S.end()) )
    {
      Lp.i_r1 = kFindInT(Lp.p1,strat); //strat->S[ifromS].s_2_r;
      Lp.i_r2 = atR;
    }
    else
    {
      /* END _ TEMPORARILY DISABLED FOR SHIFTS */
      Lp.i_r1 = -1;
      Lp.i_r2 = -1;
     }
    strat->initEcartPair(&Lp,q,p,ecartq,ecart);

    if (TEST_OPT_INTSTRATEGY)
    {
      if (!rIsPluralRing(currRing)
      && !rField_is_Ring(currRing)
      && (Lp.p->coef!=NULL))
        nDelete(&(Lp.p->coef));
    }

    strat_B(strat).push(Lp);
#ifdef CRITERION_DEBUG
    if (TEST_OPT_DEBUG) PrintS("+++ Entered pair\n");
#endif
  }
  return FALSE;
}
#endif

/*3
*(s[0], s \dot h),...,(s[k],s \dot h) will be put to the pairset L
* also the pairs (h, s\dot s[0]), ..., (h, s\dot s[k]) enter L
* additionally we put the pairs (h, s \sdot h) for s>=1 to L
*/
#ifdef HAVE_SHIFTBBA
void initenterpairsShift (poly h,int k,int ecart,int isFromQ, kStrategy strat, int atR)
{
  int h_lastVblock = pmLastVblock(h);
  assume(h_lastVblock != 0 || pLmIsConstantComp(h));
  // TODO: is it allowed to skip pairs with constants? also with constants from other components?
  if (h_lastVblock == 0) return;
  assume(pmFirstVblock(h) == 1);
  /* h comes from strat->P.p, that is LObject with LM in currRing and Tail in tailRing */
  //  atR = -1;
  if ((strat->syzComp==0)
  || (pGetComp(h)<=strat->syzComp))
  {
    int i,j;
    BOOLEAN new_pair=FALSE;

    int degbound = currRing->N/currRing->isLPring;
    int maxShift = degbound - h_lastVblock;

    if (pGetComp(h)==0)
    {
      if (strat->rightGB)
      {
        if (isFromQ)
        {
          // pairs (shifts(h),s[1..k]), (h, s[1..k])
          for (i=0; i<=maxShift; i++)
          {
            poly hh = pLPCopyAndShiftLM(h, i);
            BOOLEAN delete_hh=TRUE;
            for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
            {
              if (!strat->hasFromQ || !sjt->fromQ)
              {
                new_pair=TRUE;
                poly s = sjt->p;
                if (!enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, pmLastVblock(s), i))
                  delete_hh=FALSE;
              }
            }
            if (delete_hh) pLmDelete(hh);
          }
        }
        else
        {
          new_pair=TRUE;
          for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
          {
            poly s = sjt->p;
            if (strat->hasFromQ && sjt->fromQ)
            {
              // pairs (shifts(s[j]),h), (s[j],h)
              enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
            }
            else
            {
              // pair (h, s[j])
              enterOnePairWithoutShifts(sjt, h, s, ecart, isFromQ, strat, atR, pmLastVblock(s), 0);
            }
          }
        }
      }
      /* for Q!=NULL: build pairs (f,q),(f1,f2), but not (q1,q2)*/
      else if ((isFromQ)&&(strat->hasFromQ))
      {
        // pairs (shifts(s[1..k]),h), (s[1..k],h)
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if (!sjt->fromQ)
          {
            new_pair=TRUE;
            poly s = sjt->p;
            enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
          }
        }
        // pairs (shifts(h),s[1..k])
        if (new_pair)
        {
          for (i=1; i<=maxShift; i++)
          {
            BOOLEAN delete_hh=TRUE;
            poly hh = pLPCopyAndShiftLM(h, i);
            for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
            {
              if (!sjt->fromQ)
              {
                poly s = sjt->p;
                int s_lastVblock = pmLastVblock(s);
                if (i < s_lastVblock || (pGetComp(s) > 0 && i == s_lastVblock)) // in the module case, product criterion does not hold (note: comp h is always zero here)
                {
                  if(!enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, s_lastVblock, i))
                    delete_hh=FALSE;
                }
                else if (rField_is_Ring(currRing))
                {
                  assume(i >= s_lastVblock); // this is always the case, but just to be very sure
                  ideal fillers = id_MaxIdeal(i - s_lastVblock, currRing);
                  for (int k = 0; k < IDELEMS(fillers); k++)
                  {
                    poly hhh = pLPCopyAndShiftLM(pp_mm_Mult(h, fillers->m[k], currRing), s_lastVblock);
                    enterOnePairWithoutShifts(sjt, hhh, s, ecart, isFromQ, strat, atR, s_lastVblock, s_lastVblock);
                  }
                  idDelete(&fillers);
                }
              }
            }
            if (delete_hh) p_LmDelete(hh,currRing);
          }
        }
      }
      else
      {
        new_pair=TRUE;
        // pairs (shifts(s[1..k]),h), (s[1..k],h)
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          poly s = sjt->p;
          enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
        }
        // pairs (shifts(h),s[1..k]), (shifts(h), h)
        for (i=1; i<=maxShift; i++)
        {
          poly hh = pLPCopyAndShiftLM(h, i);
          BOOLEAN delete_hh=TRUE;
          for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
          {
            poly s = sjt->p;
            int s_lastVblock = pmLastVblock(s);
            if (i < s_lastVblock || (pGetComp(s) > 0 && i == s_lastVblock)) // in the module case, product criterion does not hold (note: comp h is always zero here)
              delete_hh=enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, s_lastVblock, i)
                && delete_hh;
            else if (rField_is_Ring(currRing))
            {
              assume(i >= s_lastVblock); // this is always the case, but just to be very sure
              ideal fillers = id_MaxIdeal(i - s_lastVblock, currRing);
              for (int k = 0; k < IDELEMS(fillers); k++)
              {
                poly hhh = pLPCopyAndShiftLM(pp_mm_Mult(h, fillers->m[k], currRing), s_lastVblock);
                enterOnePairWithoutShifts(sjt, hhh, s, ecart, isFromQ, strat, atR, s_lastVblock, s_lastVblock);
              }
              idDelete(&fillers);
            }
          }
          if (i < h_lastVblock) // in the module case, product criterion does not hold (note: comp h is always zero here)
            delete_hh=enterOnePairWithoutShifts(strat->S.end(), hh, h, ecart, isFromQ, strat, atR, h_lastVblock, i)
              && delete_hh;
          else if (rField_is_Ring(currRing))
          {
            assume(i >= h_lastVblock); // this is always the case, but just to be very sure
            ideal fillers = id_MaxIdeal(i - h_lastVblock, currRing);
            for (int k = 0; k < IDELEMS(fillers); k++)
            {
              poly hhh = pLPCopyAndShiftLM(pp_mm_Mult(h, fillers->m[k], currRing), h_lastVblock);
              enterOnePairWithoutShifts(strat->S.end(), hhh, h, ecart, isFromQ, strat, atR, h_lastVblock, h_lastVblock);
            }
            idDelete(&fillers);
          }
          if (delete_hh) pLmDelete(hh);
        }
      }
    }
    else
    {
      assume(isFromQ == 0); // an element from Q should always has 0 component
      new_pair=TRUE;
      if (strat->rightGB)
      {
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if ((pGetComp(h)==pGetComp(sjt->p))
              || (pGetComp(sjt->p)==0))
          {
            poly s = sjt->p;
            if (strat->hasFromQ && sjt->fromQ)
            {
              // pairs (shifts(s[j]),h), (s[j],h)
              enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
            }
            else
            {
              // pair (h, s[j])
              enterOnePairWithoutShifts(sjt, h, s, ecart, isFromQ, strat, atR, pmLastVblock(s), 0);
            }
          }
        }
      }
      else
      {
        // pairs (shifts(s[1..k]),h), (s[1..k],h)
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if ((pGetComp(h)==pGetComp(sjt->p))
              || (pGetComp(sjt->p)==0))
          {
            poly s = sjt->p;
            enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
          }
        }
        // pairs (shifts(h),s[1..k]), (shifts(h), h)
        for (i=1; i<=maxShift; i++)
        {
          poly hh = pLPCopyAndShiftLM(h, i);
          for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
          {
            if ((pGetComp(h)==pGetComp(sjt->p))
                || (pGetComp(sjt->p)==0))
            {
              poly s = sjt->p;
              int s_lastVblock = pmLastVblock(s);
              if (i <= s_lastVblock) // in the module case, product criterion does not hold
                enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, s_lastVblock, i);
              else if (rField_is_Ring(currRing))
              {
                assume(i >= s_lastVblock); // this is always the case, but just to be very sure
                ideal fillers = id_MaxIdeal(i - s_lastVblock, currRing);
                for (int k = 0; k < IDELEMS(fillers); k++)
                {
                  poly hhh = pLPCopyAndShiftLM(pp_mm_Mult(h, fillers->m[k], currRing), s_lastVblock);
                  enterOnePairWithoutShifts(sjt, hhh, s, ecart, isFromQ, strat, atR, s_lastVblock, s_lastVblock);
                }
                idDelete(&fillers);
              }
            }
          }
          if (i <= h_lastVblock) // in the module case, product criterion does not hold
            enterOnePairWithoutShifts(strat->S.end(), hh, h, ecart, isFromQ, strat, atR, h_lastVblock, i);
          else if (rField_is_Ring(currRing))
          {
            assume(i >= h_lastVblock); // this is always the case, but just to be very sure
            ideal fillers = id_MaxIdeal(i - h_lastVblock, currRing);
            for (int k = 0; k < IDELEMS(fillers); k++)
            {
              BOOLEAN delete_hhh=TRUE;
              poly hhh = pLPCopyAndShiftLM(pp_mm_Mult(h, fillers->m[k], currRing), h_lastVblock);
              if(!enterOnePairWithoutShifts(strat->S.end(), hhh, h, ecart, isFromQ, strat, atR, h_lastVblock, h_lastVblock))
                delete_hhh=FALSE;
              if (delete_hhh) p_LmDelete(hhh,currRing);
            }
            idDelete(&fillers);
          }
        }
      }
    }

    if (new_pair)
    {
      strat->chainCrit(h,ecart,strat);
    }
    kMergeBintoL(strat);
  }
}
#endif

/*3
*(s[0], s \dot h),...,(s[k],s \dot h) will be put to the pairset L
* also the pairs (h, s\dot s[0]), ..., (h, s\dot s[k]) enter L
* additionally we put the pairs (h, s \sdot h) for s>=1 to L
*/
#ifdef HAVE_SHIFTBBA
void initenterstrongPairsShift (poly h,int k,int ecart,int isFromQ, kStrategy strat, int atR)
{
  int h_lastVblock = pmLastVblock(h);
  assume(h_lastVblock != 0 || pLmIsConstantComp(h));
  // TODO: is it allowed to skip pairs with constants? also with constants from other components?
  if (h_lastVblock == 0) return;
  assume(pmFirstVblock(h) == 1);
  /* h comes from strat->P.p, that is LObject with LM in currRing and Tail in tailRing */
  //  atR = -1;
  if ((strat->syzComp==0)
  || (pGetComp(h)<=strat->syzComp))
  {
    int i,j;
    BOOLEAN new_pair=FALSE;

    int degbound = currRing->N/currRing->isLPring;
    int maxShift = degbound - h_lastVblock;

    if (pGetComp(h)==0)
    {
      if (strat->rightGB)
      {
        if (isFromQ)
        {
          // pairs (shifts(h),s[1..k]), (h, s[1..k])
          for (i=0; i<=maxShift; i++)
          {
            poly hh = pLPCopyAndShiftLM(h, i);
            for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
            {
              if (!strat->hasFromQ || !sjt->fromQ)
              {
                new_pair=TRUE;
                poly s = sjt->p;
                enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, pmLastVblock(s), i);
              }
            }
          }
        }
        else
        {
          new_pair=TRUE;
          for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
          {
            poly s = sjt->p;
            if (strat->hasFromQ && sjt->fromQ)
            {
              // pairs (shifts(s[j]),h), (s[j],h)
              enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
            }
            else
            {
              // pair (h, s[j])
              enterOnePairWithoutShifts(sjt, h, s, ecart, isFromQ, strat, atR, pmLastVblock(s), 0);
            }
          }
        }
      }
      /* for Q!=NULL: build pairs (f,q),(f1,f2), but not (q1,q2)*/
      else if ((isFromQ)&&(strat->hasFromQ))
      {
        // pairs (shifts(s[1..k]),h), (s[1..k],h)
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if (!sjt->fromQ)
          {
            new_pair=TRUE;
            poly s = sjt->p;
            enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
          }
        }
        // pairs (shifts(h),s[1..k])
        if (new_pair)
        {
          for (i=1; i<=maxShift; i++)
          {
            poly hh = pLPCopyAndShiftLM(h, i);
            for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
            {
              if (!sjt->fromQ)
              {
                poly s = sjt->p;
                enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, pmLastVblock(s), i);
              }
            }
          }
        }
      }
      else
      {
        new_pair=TRUE;
        // pairs (shifts(s[1..k]),h), (s[1..k],h)
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          poly s = sjt->p;
          // TODO: cache lastVblock of s[1..k] for later use
          enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
        }
        // pairs (shifts(h),s[1..k]), (shifts(h), h)
        for (i=1; i<=maxShift; i++)
        {
          poly hh = pLPCopyAndShiftLM(h, i);
          BOOLEAN delete_hh=TRUE;
          for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
          {
            poly s = sjt->p;
            if(!enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, pmLastVblock(s), i))
              delete_hh=FALSE;
          }
          if(!enterOnePairWithoutShifts(strat->S.end(), hh, h, ecart, isFromQ, strat, atR, h_lastVblock, i))
            delete_hh=FALSE;
          if (delete_hh) p_LmDelete(hh,currRing);
        }
      }
    }
    else
    {
      new_pair=TRUE;
      if (strat->rightGB)
      {
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if ((pGetComp(h)==pGetComp(sjt->p))
              || (pGetComp(sjt->p)==0))
          {
            assume(isFromQ == 0); // this case is not handled here and should also never happen
            poly s = sjt->p;
            if (strat->hasFromQ && sjt->fromQ)
            {
              // pairs (shifts(s[j]),h), (s[j],h)
              enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
            }
            else
            {
              // pair (h, s[j])
              enterOnePairWithoutShifts(sjt, h, s, ecart, isFromQ, strat, atR, pmLastVblock(s), 0);
            }
          }
        }
      }
      else
      {
        // pairs (shifts(s[1..k]),h), (s[1..k],h)
        for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
        {
          if ((pGetComp(h)==pGetComp(sjt->p))
              || (pGetComp(sjt->p)==0))
          {
            poly s = sjt->p;
            enterOnePairWithShifts(sjt, s, h, ecart, isFromQ, strat, atR, h_lastVblock, pmLastVblock(s));
          }
        }
        // pairs (shifts(h),s[1..k]), (shifts(h), h)
        for (i=1; i<=maxShift; i++)
        {
          poly hh = pLPCopyAndShiftLM(h, i);
          for (auto sjt = strat->S.begin(); sjt != strat->S.end(); ++sjt)
          {
            if ((pGetComp(h)==pGetComp(sjt->p))
                || (pGetComp(sjt->p)==0))
            {
              poly s = sjt->p;
              enterOnePairWithoutShifts(sjt, hh, s, ecart, isFromQ, strat, atR, pmLastVblock(s), i);
            }
          }
          enterOnePairWithoutShifts(strat->S.end(), hh, h, ecart, isFromQ, strat, atR, h_lastVblock, i);
        }
      }
    }

    if (new_pair)
    {
      strat->chainCrit(h,ecart,strat);
    }
    kMergeBintoL(strat);
  }
}
#endif

/*2
*(s[0],h),...,(s[k],h) will be put to the pairset L(via initenterpairs)
*superfluous elements in S will be deleted
*/
#ifdef HAVE_SHIFTBBA
void enterpairsShift (poly h,int k,int ecart,sBasisSet::iterator pos,kStrategy strat, int atR)
{
  /* h is strat->P.p, that is LObject with LM in currRing and Tail in tailRing */
  /* Q: what is exactly the strat->fromT ? A: a local case trick; don't need it yet*/

  /* if (!(rField_is_Domain(currRing))) enterExtendedSpoly(h, strat); */ // TODO: enterExtendedSpoly not for LP yet
  initenterpairsShift(h,k,ecart,0,strat, atR);
  if ( (!strat->fromT)
  && ((strat->syzComp==0)
    ||(pGetComp(h)<=strat->syzComp)))
  {
    // Caller contract: k == strat->S.size()-1 on entry.
    assume(k == strat->S.size() - 1);
    unsigned long h_sev = pGetShortExpVector(h);
    for (auto it = pos; it != strat->S.end(); ++it)
    {
      // TODO this currently doesn't clear all possible elements because of commutative division
      if (strat->rightGB && strat->hasFromQ && it->fromQ) continue;
      strat->S.clear_if_divisible(h, h_sev, it, strat);
    }
  }
}
#endif

/*2
* enteres all admissible shifts of p into T
* assumes that p is already in T!
*/
#ifdef HAVE_SHIFTBBA
void enterTShift(LObject p, kStrategy strat, int atT)
{
  /* determine how many elements we have to insert */
  /* x(0)y(1)z(2) : lastVblock-1=2, to add until lastVblock=uptodeg-1 */
  /* hence, a total number of elt's to add is: */
  /*  int toInsert = 1 + (uptodeg-1) - (pLastVblock(p.p, lV) -1);  */
  pAssume(p.p != NULL);

  int maxPossibleShift = p_mLPmaxPossibleShift(p.p, strat->tailRing);

  for (int i = 1; i <= maxPossibleShift; i++)
  {
    LObject qq;
    qq.p = pLPCopyAndShiftLM(p.p, i); // don't use Set() because it'll test the poly order
    qq.shift = i;
    strat->initEcart(&qq); // initEcartBBA sets length, pLength, FDeg and ecart

    enterT(qq, strat, atT); // enterT is modified, so it doesn't copy and delete the tail of shifted polys
  }
}
#endif

#ifdef HAVE_SHIFTBBA
poly redtailBbaShift (LObject* L, sBasisSet::const_iterator end, kStrategy strat, BOOLEAN withT, BOOLEAN normalize)
{
  /* for the shift case need to run it with withT = TRUE */
  strat->redTailChange=FALSE;
  if (strat->noTailReduction) return L->GetLmCurrRing();
  poly h, p;
  p = h = L->GetLmTailRing();
  if ((h==NULL) || (pNext(h)==NULL))
    return L->GetLmCurrRing();

  TObject* With;
  // placeholder in case strat->T.empty()
  TObject  With_s(strat->tailRing);

  LObject Ln(pNext(h), strat->tailRing);
  Ln.pLength = L->GetpLength() - 1;

  pNext(h) = NULL;
  if (L->p != NULL) pNext(L->p) = NULL;
  L->pLength = 1;

  Ln.PrepareRed(strat->use_buckets);

  while(!Ln.IsNull())
  {
    loop
    {
      Ln.SetShortExpVector();
      if (withT)
      {
        int j;
        j = kFindDivisibleByInT(strat, &Ln);
        if (j < 0) break;
        With = &(strat->T[j]);
      }
      else
      {
        With = kFindDivisibleByInS_T(strat, end, &Ln, &With_s);
        if (With == NULL) break;
      }
      if (normalize && (!TEST_OPT_INTSTRATEGY) && (!nIsOne(pGetCoeff(With->p))))
      {
        With->pNorm();
        //if (TEST_OPT_PROT) { PrintS("n"); mflush(); }
      }
      strat->redTailChange=TRUE;
      if (ksReducePolyTail(L, With, &Ln))
      {
        // reducing the tail would violate the exp bound
        //  set a flag and hope for a retry (in bba)
        strat->completeReduce_retry=TRUE;
        if ((Ln.p != NULL) && (Ln.t_p != NULL)) Ln.p=NULL;
        do
        {
          pNext(h) = Ln.LmExtractAndIter();
          pIter(h);
          L->pLength++;
        } while (!Ln.IsNull());
        goto all_done;
      }
      if (Ln.IsNull()) goto all_done;
      if (! withT) With_s.Init(currRing);
    }
    pNext(h) = Ln.LmExtractAndIter();
    pIter(h);
    L->pLength++;
  }

  all_done:
  Ln.Delete();
  if (L->p != NULL) pNext(L->p) = pNext(p);

  if (strat->redTailChange)
  {
    L->length = 0;
  }
  L->Normalize(); // HANNES: should have a test
  kTest_L(L,strat);
  return L->GetLmCurrRing();
}
#endif
