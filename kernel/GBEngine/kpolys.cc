


#include "kernel/mod2.h"

#include "kernel/polys.h"

#ifdef HAVE_AVX2
#include <immintrin.h>
#endif

/* Returns TRUE if
     * LM(p) | LM(lcm)
     * LC(p) | LC(lcm) only if ring
     * Exists i, j:
         * LE(p, i)  != LE(lcm, i)
         * LE(p1, i) != LE(lcm, i)   ==> LCM(p1, p) != lcm
         * LE(p, j)  != LE(lcm, j)
         * LE(p2, j) != LE(lcm, j)   ==> LCM(p2, p) != lcm
*/

#ifdef HAVE_AVX2
BOOLEAN pCompareChain_16bit_AVX2 (poly p,poly p1,poly p2,poly lcm, const ring R)
{
  __m256i * p_exp_ptr = (__m256i *) p->exp;
  __m256i * p1_exp_ptr = (__m256i *) p1->exp;
  __m256i * p2_exp_ptr = (__m256i *) p2->exp;
  __m256i * lcm_exp_ptr = (__m256i *) lcm->exp;
  __m256i * VarL_Bitmask_ptr = (__m256i *) R->VarL_Bitmask;
  int p_diff_count = 0;
  int p1_diff_count = 0;
  int p2_diff_count = 0;
  int p1_and_p2_common_diff_count = 0;

  for (int i=0; i < R->Exp_SIMD_Size; i++) {
    // unaligned loads because I haven't been able to get exponent fields aligned on a SIMD_VECTOR_SIZE boundary
    __m256i p_exp = _mm256_loadu_si256 (p_exp_ptr + i);
    __m256i p1_exp = _mm256_loadu_si256 (p1_exp_ptr + i);
    __m256i p2_exp = _mm256_loadu_si256 (p2_exp_ptr + i);
    __m256i lcm_exp = _mm256_loadu_si256 (lcm_exp_ptr + i);
    __m256i VarL_Bitmask = _mm256_loadu_si256 (VarL_Bitmask_ptr + i);

    // Basic divisibility check
    // if (p_exp[i] > lcm_exp[i]) return FALSE;
    __m256i result = _mm256_cmpgt_epi16(p_exp, lcm_exp);
    result = _mm256_and_si256(result, VarL_Bitmask);
    if (! _mm256_testz_si256(result, result)) return FALSE;

    // Compute difference indicators
    __m256i p_diff = _mm256_cmpeq_epi16(p_exp, lcm_exp);
    __m256i p1_diff = _mm256_cmpeq_epi16(p1_exp, lcm_exp);
    __m256i p2_diff = _mm256_cmpeq_epi16(p2_exp, lcm_exp);

    // invert the eq to neq, and mask off the variables
    p_diff = _mm256_andnot_si256(p_diff, VarL_Bitmask);
    p1_diff = _mm256_andnot_si256(p1_diff, VarL_Bitmask);
    p2_diff = _mm256_andnot_si256(p2_diff, VarL_Bitmask);

    p1_diff = _mm256_and_si256(p1_diff, p_diff);
    p2_diff = _mm256_and_si256(p2_diff, p_diff);

    __m256i p1_and_p2_common_diff = _mm256_and_si256(p1_diff, p2_diff);

    // we count high bits in epi8, but our exponents are epu16, so each TRUE gets counted twice, so we divide by 2

    p_diff_count += _mm_popcnt_u32(_mm256_movemask_epi8(p_diff)) / 2;
    p1_diff_count += _mm_popcnt_u32(_mm256_movemask_epi8(p1_diff)) / 2;
    p2_diff_count += _mm_popcnt_u32(_mm256_movemask_epi8(p2_diff)) / 2;

    p1_and_p2_common_diff_count += _mm_popcnt_u32(_mm256_movemask_epi8(p1_and_p2_common_diff)) / 2;
  }

  // Chain criterion needs at least 2 variables where p differs from lcm
  if (p_diff_count <= 1) return FALSE;

  // If p1 or p2 equals lcm everywhere, chain criterion cannot apply
  if (p1_diff_count == 0 || p2_diff_count == 0) return FALSE;

  // If at least one of them differs in two places (and the other differs in at least one place), chain criterion applies
  if (p1_diff_count > 1 && p2_diff_count > 1) return TRUE;

  // p1 and p2 differ from lcm in only one variable.  If it's the same variable, chain criteron cannot apply
  if (p1_and_p2_common_diff_count == 1) return FALSE;

  // p1 and p2 differ from lcm in only one variable and it's two different variables, chain criteron applies
  return TRUE;
}
#endif

#if 0
BOOLEAN pCompareChain_16bit (poly p,poly p1,poly p2,poly lcm, const ring R)
{
  int k, j;
  const int N = R->N;

  if (lcm==NULL) return FALSE;

  if (pGetComp(p) != pGetComp(lcm)) return FALSE;

  unsigned short * p_exp = p->exp;
  unsigned short * p1_exp = p1->exp;
  unsigned short * p2_exp = p2->exp;
  unsigned short * lcm_exp = lcm->exp;

  // Early termination optimization: Track which variables differ
  bool p_diff[4 * R->ExpL_Size];   // Where p != lcm
  bool p1_diff[4 * R->ExpL_Size];  // Where p1 != lcm
  bool p2_diff[4 * R->ExpL_Size];  // Where p2 != lcm

  int p_diff_count = 0;
  int p1_diff_count = 0;
  int p2_diff_count = 0;

  for (int i = 0; i < 4 * R->ExpL_Size; i++) {
    // Original divisibility check
    if (p_exp[i] > lcm_exp[i]) return FALSE;

    // Compute difference indicators
    p_diff[i] = (p_exp[i] != lcm_exp[i]);
    p1_diff[i] = (p1_exp[i] != lcm_exp[i]);
    p2_diff[i] = (p2_exp[i] != lcm_exp[i]);

    if (p_diff[i]) p_diff_count++;
    if (p1_diff[i]) p1_diff_count++;
    if (p2_diff[i]) p2_diff_count++;
  }

  // Early termination checks based on difference patterns
  // Chain criterion needs at least 2 variables where p differs from lcm
  if (p_diff_count <= 1) return FALSE;

  // If p1 or p2 equals lcm everywhere, chain criterion cannot apply
  if (p1_diff_count == 0 || p2_diff_count == 0) return FALSE;

  // Main chain criterion loop - use pre-computed difference arrays
  for (j=N; j; j--)
  {
    if (p1_diff[j])  // p1 differs from lcm at position j
    {
      if (p_diff[j])  // p also differs at j
      {
        // Search for k where both p and p2 differ from lcm
        for (k=N; k>j; k--)
        {
          if (p_diff[k] && p2_diff[k])
            return TRUE;
        }
        for (k=j-1; k; k--)
        {
          if (p_diff[k] && p2_diff[k])
            return TRUE;
        }
        return FALSE;
      }
    }
    else if (p2_diff[j])  // p2 differs from lcm at position j (but not p1)
    {
      if (p_diff[j])  // p also differs at j
      {
        // Search for k where both p and p1 differ from lcm
        for (k=N; k>j; k--)
        {
          if (p_diff[k] && p1_diff[k])
            return TRUE;
        }
        for (k=j-1; k!=0 ; k--)
        {
          if (p_diff[k] && p1_diff[k])
            return TRUE;
        }
        return FALSE;
      }
    }
  }
  return FALSE;
}
#endif

BOOLEAN pCompareChain (poly p,poly p1,poly p2,poly lcm, const ring R)
{
  int k, j;
  const int N = R->N;

  if (lcm==NULL) return FALSE;
  if (pGetComp(p) != pGetComp(lcm)) return FALSE;

#ifdef HAVE_AVX2
  if (R->BitsPerExp == 16) return pCompareChain_16bit_AVX2(p, p1, p2, lcm, R);
#endif

  // Optimization: Cache all exponents to eliminate redundant p_GetExp calls
  // Pre-allocate arrays on stack (typically N <= 100, so ~3.2KB max)
  long p_exp[N+1];
  long p1_exp[N+1];
  long p2_exp[N+1];
  long lcm_exp[N+1];

  // Early termination optimization: Track which variables differ
  bool p_diff[N+1];   // Where p != lcm
  bool p1_diff[N+1];  // Where p1 != lcm
  bool p2_diff[N+1];  // Where p2 != lcm

  int p_diff_count = 0;
  int p1_diff_count = 0;
  int p2_diff_count = 0;

  // Fetch each exponent exactly once and compute difference indicators
  for (int i = 1; i <= N; i++) {
    p_exp[i] = p_GetExp(p, i, R);
    p1_exp[i] = p_GetExp(p1, i, R);
    p2_exp[i] = p_GetExp(p2, i, R);
    lcm_exp[i] = p_GetExp(lcm, i, R);

    p_diff[i] = (p_exp[i] != lcm_exp[i]);
    p1_diff[i] = (p1_exp[i] != lcm_exp[i]);
    p2_diff[i] = (p2_exp[i] != lcm_exp[i]);

    if (p_diff[i]) p_diff_count++;
    if (p1_diff[i]) p1_diff_count++;
    if (p2_diff[i]) p2_diff_count++;
  }

  // Early termination checks based on difference patterns
  // Chain criterion needs at least 2 variables where p differs from lcm
  if (p_diff_count <= 1) return FALSE;

  // If p1 or p2 equals lcm everywhere, chain criterion cannot apply
  if (p1_diff_count == 0 || p2_diff_count == 0) return FALSE;

  // Original divisibility check
  for (j=N; j; j--)
    if (p_exp[j] > lcm_exp[j]) return FALSE;

  // Main chain criterion loop - use pre-computed difference arrays
  for (j=N; j; j--)
  {
    if (p1_diff[j])  // p1 differs from lcm at position j
    {
      if (p_diff[j])  // p also differs at j
      {
        // Search for k where both p and p2 differ from lcm
        for (k=N; k>j; k--)
        {
          if (p_diff[k] && p2_diff[k])
            return TRUE;
        }
        for (k=j-1; k; k--)
        {
          if (p_diff[k] && p2_diff[k])
            return TRUE;
        }
        return FALSE;
      }
    }
    else if (p2_diff[j])  // p2 differs from lcm at position j (but not p1)
    {
      if (p_diff[j])  // p also differs at j
      {
        // Search for k where both p and p1 differ from lcm
        for (k=N; k>j; k--)
        {
          if (p_diff[k] && p1_diff[k])
            return TRUE;
        }
        for (k=j-1; k!=0 ; k--)
        {
          if (p_diff[k] && p1_diff[k])
            return TRUE;
        }
        return FALSE;
      }
    }
  }
  return FALSE;
}

#ifdef HAVE_RATGRING
BOOLEAN pCompareChainPart (poly p,poly p1,poly p2,poly lcm, const ring R)
{
  int k, j;

  if (lcm==NULL) return FALSE;

  for (j=R->real_var_end; j>=R->real_var_start; j--)
    if ( p_GetExp(p,j, R) >  p_GetExp(lcm,j, R)) return FALSE;
  if ( pGetComp(p) !=  pGetComp(lcm)) return FALSE;
  for (j=R->real_var_end; j>=R->real_var_start; j--)
  {
    if (p_GetExp(p1,j, R)!=p_GetExp(lcm,j, R))
    {
      if (p_GetExp(p,j, R)!=p_GetExp(lcm,j, R))
      {
        for (k=(R->N); k>j; k--)
        for (k=R->real_var_end; k>j; k--)
        {
          if ((p_GetExp(p,k, R)!=p_GetExp(lcm,k, R))
          && (p_GetExp(p2,k, R)!=p_GetExp(lcm,k, R)))
            return TRUE;
        }
        for (k=j-1; k>=R->real_var_start; k--)
        {
          if ((p_GetExp(p,k, R)!=p_GetExp(lcm,k, R))
          && (p_GetExp(p2,k, R)!=p_GetExp(lcm,k, R)))
            return TRUE;
        }
        return FALSE;
      }
    }
    else if (p_GetExp(p2,j, R)!=p_GetExp(lcm,j, R))
    {
      if (p_GetExp(p,j, R)!=p_GetExp(lcm,j, R))
      {
        for (k=R->real_var_end; k>j; k--)
        {
          if ((p_GetExp(p,k, R)!=p_GetExp(lcm,k, R))
          && (p_GetExp(p1,k, R)!=p_GetExp(lcm,k, R)))
            return TRUE;
        }
        for (k=j-1; k>=R->real_var_start; k--)
        {
          if ((p_GetExp(p,k, R)!=p_GetExp(lcm,k, R))
          && (p_GetExp(p1,k, R)!=p_GetExp(lcm,k, R)))
            return TRUE;
        }
        return FALSE;
      }
    }
  }
  return FALSE;
}
#endif

