/*
 * singrust — bridge between Singular and the rustgb Gröbner-basis
 * engine. Exposes one Singular procedure `rustgb_std(ideal) -> ideal`
 * that computes the reduced Gröbner basis via rustgb. The rustgb
 * engine supports a narrow slice of rings (Z/p, degrevlex, <=31 vars,
 * commutative, global ordering); unsupported inputs are rejected
 * with a clear error rather than a crash.
 *
 * Reference: ~/project/docs/rust-bba-port-plan.md §§11–12 and
 * ~/project/docs/rustgb-singular-ffi-report.md.
 *
 * This file is modelled on `singmathic.cc` — same shape of
 * input/output streaming, same iiAddCproc registration — but the
 * actual GB engine is the rustgb cdylib (see include/rustgb.h in
 * the rustgb crate).
 */

#include "kernel/mod2.h"

#ifdef HAVE_RUSTGB

#include "misc/options.h"

#include "kernel/ideals.h"
#include "kernel/polys.h"

#include "Singular/ipid.h"
#include "Singular/feOpt.h"
#include "Singular/mod_lib.h"

#include <rustgb.h>

#include <vector>
#include <cstring>

namespace {

/*
 * Validate the current ring against rustgb's supported subset.
 * Returns TRUE on success, FALSE on error (with WerrorS already
 * called).
 *
 * Requirements enforced here:
 *   * commutative (no Plural / Letterplace),
 *   * no module / quotient-ring stuff — the ideal is in a bare
 *     polynomial ring,
 *   * coefficient ring is Z/p with p prime and p < 2^31,
 *   * at least 1 variable, at most 31 (rustgb's packing limit),
 *   * ordering is plain dp (degrevlex), i.e. `rOrd_is_dp` holds
 *     OR nvars == 1 and the order block is dp (`rOrd_is_dp`
 *     requires nvars > 1 so we handle that corner separately).
 */
BOOLEAN rustgb_validate_ring(ring r)
{
  if (r == NULL)
  {
    WerrorS("rustgb_std: no current ring");
    return FALSE;
  }
  if (rIsPluralRing(r))
  {
    WerrorS("rustgb_std: non-commutative (Plural) rings are not supported");
    return FALSE;
  }
#ifdef HAVE_RINGS
  if (rField_is_Ring(r))
  {
    WerrorS("rustgb_std: coefficient rings (Z, Z/n for n composite, etc.) are not supported");
    return FALSE;
  }
#endif
  if (!rField_is_Zp(r))
  {
    WerrorS("rustgb_std: coefficient field must be Z/p (small prime). Q, Z, algebraic extensions, GF(q) are not supported");
    return FALSE;
  }
  const int ch = n_GetChar(r->cf);
  if (ch < 2)
  {
    Werror("rustgb_std: unexpected characteristic %d", ch);
    return FALSE;
  }
  // p < 2^31 — npInt stores the prime in `ch` as a positive int, so
  // this matches directly.
  if ((unsigned int)ch >= (1u << 31))
  {
    WerrorS("rustgb_std: prime must be < 2^31");
    return FALSE;
  }
  if (!rHasGlobalOrdering(r))
  {
    WerrorS("rustgb_std: ring ordering must be global (degrevlex / dp)");
    return FALSE;
  }
  const int nvars = r->N;
  if (nvars < 1 || nvars > 31)
  {
    Werror("rustgb_std: nvars = %d is out of range (need 1..=31)", nvars);
    return FALSE;
  }

  // `rOrd_is_dp` requires nvars > 1. For nvars == 1, any global
  // single-variable order is equivalent to dp, so we accept any
  // order block we can identify as dp / Dp / lp. To keep the
  // contract narrow ("dp only") we check the raw block tags.
  //
  // Skip an optional leading C/c component block, then the primary
  // block must be ringorder_dp spanning all variables.
  int ord = 0;
  if (r->order[0] == ringorder_C || r->order[0] == ringorder_c) ord = 1;
  if (r->order[ord] != ringorder_dp)
  {
    WerrorS("rustgb_std: ordering must be dp (degrevlex). Other orderings (lp, rp, Dp, weights, mixed, block) are not supported");
    return FALSE;
  }
  if (r->block0[ord] != 1 || r->block1[ord] != nvars)
  {
    WerrorS("rustgb_std: dp block must span all variables (no mixed block orders)");
    return FALSE;
  }
  // Trailing block should be the end sentinel or a C/c component.
  rRingOrder_t trailing = static_cast<rRingOrder_t>(r->order[ord + 1]);
  if (trailing != ringorder_no && trailing != ringorder_C && trailing != ringorder_c)
  {
    WerrorS("rustgb_std: ordering has a second block; only single-block dp is supported");
    return FALSE;
  }

  return TRUE;
}

/*
 * Extract a Zp coefficient in the canonical range [0, p).
 *
 * `n_Int` returns the representing integer in (-p/2, p/2] for Zp;
 * we add p if negative. This avoids reaching into Zp internals
 * (which the spec forbids).
 */
uint32_t extract_coeff_zp(number c, const coeffs cf, int p)
{
  number cc = c;
  long v = n_Int(cc, cf);
  if (v < 0) v += p;
  // Defensive: clip in case n_Int returned something unexpected.
  if (v < 0 || v >= p)
  {
    v = ((v % p) + p) % p;
  }
  return (uint32_t) v;
}

} // namespace

/*
 * rustgb_std(I): compute the reduced Gröbner basis of `I` using
 * the rustgb engine. Returns a new ideal on success, TRUE (error)
 * on failure.
 */
static BOOLEAN rustgb_std(leftv result, leftv arg)
{
  result->rtyp = NONE;
  result->data = NULL;

  if (arg == NULL || arg->next != NULL || arg->Typ() != IDEAL_CMD)
  {
    WerrorS("Syntax: rustgb_std(<ideal>)");
    return TRUE;
  }
  ring r = currRing;
  if (!rustgb_validate_ring(r))
  {
    return TRUE;
  }

  const int p = n_GetChar(r->cf);
  const int nvars = r->N;
  const ideal I = (ideal) arg->Data();

  // Build rustgb ring.
  rustgb_ring* rring = rustgb_ring_create((uint32_t) nvars, (uint32_t) p);
  if (rring == NULL)
  {
    Werror("rustgb_std: rustgb_ring_create failed: %s", rustgb_last_error());
    return TRUE;
  }

  rustgb_input* rin = rustgb_input_begin(rring);
  if (rin == NULL)
  {
    rustgb_ring_destroy(rring);
    Werror("rustgb_std: rustgb_input_begin failed: %s", rustgb_last_error());
    return TRUE;
  }

  // Scratch buffer for exponent vectors. Reused across terms.
  std::vector<int32_t> exps(nvars);

  const int size = IDELEMS(I);
  for (int i = 0; i < size; ++i)
  {
    if (rustgb_input_poly_begin(rin) != 0)
    {
      Werror("rustgb_std: poly_begin failed: %s", rustgb_last_error());
      rustgb_input_destroy(rin);
      rustgb_ring_destroy(rring);
      return TRUE;
    }

    for (poly q = I->m[i]; q != NULL; q = pNext(q))
    {
      // Exponents: Singular's pGetExp is 1-based; rustgb uses
      // index 0..nvars-1.
      for (int v = 1; v <= nvars; ++v)
      {
        exps[v - 1] = (int32_t) pGetExp(q, v);
      }
      number c = pGetCoeff(q);
      const uint32_t coeff = extract_coeff_zp(c, r->cf, p);
      if (coeff == 0) continue; // defensive; pGetCoeff should be nonzero
      if (rustgb_input_term(rin, exps.data(), coeff) != 0)
      {
        Werror("rustgb_std: input_term failed: %s", rustgb_last_error());
        rustgb_input_destroy(rin);
        rustgb_ring_destroy(rring);
        return TRUE;
      }
    }

    if (rustgb_input_poly_end(rin) != 0)
    {
      Werror("rustgb_std: poly_end failed: %s", rustgb_last_error());
      rustgb_input_destroy(rin);
      rustgb_ring_destroy(rring);
      return TRUE;
    }
  }

  // Run the GB. `rustgb_compute` consumes `rin` regardless of
  // success or failure — do NOT call rustgb_input_destroy after.
  rustgb_basis* rb = rustgb_compute(rin);
  if (rb == NULL)
  {
    Werror("rustgb_std: compute failed: %s", rustgb_last_error());
    rustgb_ring_destroy(rring);
    return TRUE;
  }

  // Convert rustgb_basis -> Singular ideal.
  const size_t n_poly = rustgb_basis_poly_count(rb);
  ideal out = idInit((int) n_poly, 1);

  for (size_t pi = 0; pi < n_poly; ++pi)
  {
    poly head = NULL;
    poly tail = NULL;

    // Open an opaque iterator over this poly's terms. The random-
    // access rustgb_basis_term was removed in favour of a cursor
    // so the rustgb side can later move to a linked-list poly
    // without reshaping the C surface.
    rustgb_term_iter* it = rustgb_term_iter_open(rb, pi);
    if (it == NULL)
    {
      Werror("rustgb_std: term_iter_open failed: %s", rustgb_last_error());
      idDelete(&out);
      rustgb_basis_destroy(rb);
      rustgb_ring_destroy(rring);
      return TRUE;
    }

    for (;;)
    {
      uint32_t coeff_out = 0;
      const int rc = rustgb_term_iter_next(it, exps.data(), &coeff_out);
      if (rc == 1) break;          // exhausted
      if (rc != 0)                  // error (rc == 2)
      {
        Werror("rustgb_std: term_iter_next failed: %s", rustgb_last_error());
        // Best-effort cleanup: drop partial polys, ideal.
        if (head != NULL) p_Delete(&head, r);
        rustgb_term_iter_close(it);
        idDelete(&out);
        rustgb_basis_destroy(rb);
        rustgb_ring_destroy(rring);
        return TRUE;
      }

      poly t = p_Init(r);
      // Exponents: rustgb gives 0-based, Singular stores 1-based.
      for (int v = 1; v <= nvars; ++v)
      {
        pSetExp(t, v, (long) exps[v - 1]);
      }
      // Coefficient: rustgb gives [0, p); n_Init takes long, Zp
      // handles canonicalisation.
      pSetCoeff0(t, n_Init((long) coeff_out, r->cf));
      p_Setm(t, r);

      if (head == NULL)
      {
        head = t;
        tail = t;
      }
      else
      {
        pNext(tail) = t;
        tail = t;
      }
    }

    rustgb_term_iter_close(it);
    out->m[pi] = head;
  }

  rustgb_basis_destroy(rb);
  rustgb_ring_destroy(rring);

  result->rtyp = IDEAL_CMD;
  result->data = (void*) out;
  return FALSE;
}

/*
 * rustgb_version(): returns the rustgb library's version string.
 * Useful smoke test from within Singular.
 */
static BOOLEAN rustgb_version_proc(leftv result, leftv arg)
{
  (void) arg;
  const char* v = rustgb_version();
  result->rtyp = STRING_CMD;
  result->data = (void*) omStrDup(v);
  return FALSE;
}

extern "C" int SI_MOD_INIT(singrust)(SModulFunctions* psModulFunctions)
{
  psModulFunctions->iiAddCproc(
    (currPack->libname ? currPack->libname : ""),
    "rustgb_std",
    FALSE,
    rustgb_std);
  psModulFunctions->iiAddCproc(
    (currPack->libname ? currPack->libname : ""),
    "rustgb_version",
    FALSE,
    rustgb_version_proc);
  return MAX_TOK;
}

#endif /* HAVE_RUSTGB */
