/**
 * @file kevlog_wrap.h
 * @brief pLmFree / pDelete wrappers for leak-on-trace mode (task 325).
 *
 * Include this AFTER kernel/polys.h in the three files that are
 * event-log hot:
 *   kernel/GBEngine/kthread.cc
 *   kernel/GBEngine/kutil.cc
 *   kernel/GBEngine/kstd2.cc
 *
 * When g_defer_frees is true (SINGULAR_EVENT_LOG=1), each pLmFree /
 * pDelete / p_LmFree / p_Delete call is turned into a NO-OP.  We
 * simply do not free.  Poly pointers captured in events therefore
 * stay valid until the dump reads them for the polys.txt file.
 * Process exit reclaims the leaked memory — fine for a 0.2s
 * reproducer.
 *
 * The wrappers preserve the exact signatures of the original
 * overloaded inline functions in kernel/polys.h:
 *     pLmFree(poly)    and  pLmFree(poly *)
 *     p_LmFree(poly, ring)      and  p_LmFree(poly *, ring)
 *     p_Delete(poly *, ring)    and  p_Delete(poly *, ring, ring)
 *     pDelete(poly *)  (macro)
 *
 * In the poly* overloads we still set *p = NULL so caller null-checks
 * continue to work — otherwise double-skipped frees leave a stale
 * pointer that the caller may re-free later.
 *
 * When g_defer_frees is false (untraced runs), these delegate to the
 * originals: a single-load branch, zero-cost when the log is disabled.
 */

#ifndef KEVLOG_WRAP_H
#define KEVLOG_WRAP_H

#include "kernel/GBEngine/kevlog.h"

/* ------------------------------------------------------------------ */
/*  pLmFree wrappers (overloaded)                                      */
/* ------------------------------------------------------------------ */
static inline void kt_pLmFree(poly p)
{
  if (g_defer_frees) return;   /* no-op: leak and keep the pointer live */
  p_LmFree(p, currRing);
}

static inline void kt_pLmFree(poly *p)
{
  if (g_defer_frees) {
    /* p_LmFree(poly*) both frees head AND advances *p to pNext(*p).
     * Callers rely on the advance.  We leak the head but still
     * perform the advance so the traversal semantics are preserved. */
    if (p != NULL && *p != NULL) *p = (*p)->next;
    return;
  }
  p_LmFree(p, currRing);
}

/* ------------------------------------------------------------------ */
/*  p_LmFree wrappers (overloaded, ring-aware)                         */
/* ------------------------------------------------------------------ */
static inline void kt_p_LmFree(poly p, ring r)
{
  if (g_defer_frees) return;
  p_LmFree(p, r);
}

static inline void kt_p_LmFree(poly *p, ring r)
{
  if (g_defer_frees) {
    if (p != NULL && *p != NULL) *p = (*p)->next;
    return;
  }
  p_LmFree(p, r);
}

/* ------------------------------------------------------------------ */
/*  pDelete / p_Delete wrappers (full-chain free)                      */
/* ------------------------------------------------------------------ */
static inline void kt_pDelete(poly *pp)
{
  if (g_defer_frees) {
    /* pDelete normally sets *pp = NULL after freeing the chain.
     * We mirror that to preserve caller expectations (they check
     * *pp for NULL to know the ideal has been cleared). */
    if (pp != NULL) *pp = NULL;
    return;
  }
  p_Delete(pp, currRing);
}

static inline void kt_p_Delete(poly *pp, const ring r)
{
  if (g_defer_frees) {
    if (pp != NULL) *pp = NULL;
    return;
  }
  p_Delete(pp, r);
}

static inline void kt_p_Delete(poly *pp, const ring lmRing, const ring tailRing)
{
  if (g_defer_frees) {
    if (pp != NULL) *pp = NULL;
    return;
  }
  p_Delete(pp, lmRing, tailRing);
}

/* ------------------------------------------------------------------ */
/*  Redirect the names used at call sites in kthread/kutil/kstd2.      */
/* ------------------------------------------------------------------ */
#define pLmFree  kt_pLmFree
#define pDelete  kt_pDelete
#define p_LmFree kt_p_LmFree
#define p_Delete kt_p_Delete

#endif /* KEVLOG_WRAP_H */
