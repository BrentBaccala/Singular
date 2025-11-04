


#include "kernel/mod2.h"

#include "kernel/polys.h"

/* Returns TRUE if
     * LM(p) | LM(lcm)
     * LC(p) | LC(lcm) only if ring
     * Exists i, j:
         * LE(p, i)  != LE(lcm, i)
         * LE(p1, i) != LE(lcm, i)   ==> LCM(p1, p) != lcm
         * LE(p, j)  != LE(lcm, j)
         * LE(p2, j) != LE(lcm, j)   ==> LCM(p2, p) != lcm
*/
BOOLEAN pCompareChain (poly p,poly p1,poly p2,poly lcm, const ring R)
{
  int k, j;
  const int N = R->N;

  if (lcm==NULL) return FALSE;

  // Optimization: Cache all exponents to eliminate redundant p_GetExp calls
  // Pre-allocate arrays on stack (typically N <= 100, so ~3.2KB max)
  long p_exp[N+1];
  long p1_exp[N+1];
  long p2_exp[N+1];
  long lcm_exp[N+1];

  // Fetch each exponent exactly once - O(N) instead of O(N²) or O(N³)
  for (int i = 1; i <= N; i++) {
    p_exp[i] = p_GetExp(p, i, R);
    p1_exp[i] = p_GetExp(p1, i, R);
    p2_exp[i] = p_GetExp(p2, i, R);
    lcm_exp[i] = p_GetExp(lcm, i, R);
  }

  // Original logic, now using cached values (simple array lookups)
  for (j=N; j; j--)
    if (p_exp[j] > lcm_exp[j]) return FALSE;
  if (pGetComp(p) != pGetComp(lcm)) return FALSE;

  for (j=N; j; j--)
  {
    if (p1_exp[j] != lcm_exp[j])
    {
      if (p_exp[j] != lcm_exp[j])
      {
        for (k=N; k>j; k--)
        {
          if ((p_exp[k] != lcm_exp[k])
          && (p2_exp[k] != lcm_exp[k]))
            return TRUE;
        }
        for (k=j-1; k; k--)
        {
          if ((p_exp[k] != lcm_exp[k])
          && (p2_exp[k] != lcm_exp[k]))
            return TRUE;
        }
        return FALSE;
      }
    }
    else if (p2_exp[j] != lcm_exp[j])
    {
      if (p_exp[j] != lcm_exp[j])
      {
        for (k=N; k>j; k--)
        {
          if ((p_exp[k] != lcm_exp[k])
          && (p1_exp[k] != lcm_exp[k]))
            return TRUE;
        }
        for (k=j-1; k!=0 ; k--)
        {
          if ((p_exp[k] != lcm_exp[k])
          && (p1_exp[k] != lcm_exp[k]))
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

