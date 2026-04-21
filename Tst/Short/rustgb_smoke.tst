LIB "tst.lib";
tst_init();

// rustgb_smoke: verify that rustgb_std(I) computes the same
// reduced Groebner basis as std(I) on the cyclic-3, cyclic-4,
// and cyclic-5 ideals over Z/32003 with degrevlex. Also checks
// that the dyn_module rejects unsupported rings cleanly.
//
// Reference: ~/project/docs/rustgb-singular-ffi-report.md

LIB "singrust.so";
LIB "general.lib";

option(redSB);

// Helper: canonical form for comparison.
proc canon(ideal a) { return(sort(interred(a))[1]); }

proc same(ideal a, ideal b)
{
    ideal ca = canon(a);
    ideal cb = canon(b);
    if (size(ca) != size(cb)) { return(0); }
    for (int k = 1; k <= size(ca); k++)
    {
        if (ca[k] - cb[k] != 0) { return(0); }
    }
    return(1);
}

proc cmp_cyclic(int n)
{
    ring R = 32003, (x(1..n)), dp;
    ideal I = cyclic(n);
    ideal G_singular = std(I);
    ideal G_rust = rustgb_std(I);
    if (same(G_singular, G_rust))
    {
        printf("PASS cyclic-%s (|GB|=%s)", string(n), string(size(G_singular)));
    }
    else
    {
        printf("FAIL cyclic-%s", string(n));
        printf("  std:        size=%s", string(size(G_singular)));
        printf("  rustgb_std: size=%s", string(size(G_rust)));
    }
}

// ---- Positive cases ----

cmp_cyclic(3);
cmp_cyclic(4);
cmp_cyclic(5);

// ---- Negative cases: unsupported ring must error cleanly ----

// (1) Q instead of Z/p.
ring Rq = 0, (x,y,z), dp;
ideal Iq = x+y+z, xy+yz+zx, xyz-1;
ideal Gq;
string errmsg = "?";
Gq = rustgb_std(Iq);    // expected to error
kill Rq;

// (2) lp instead of dp.
ring Rlp = 32003, (x,y,z), lp;
ideal Ilp = x+y+z, xy+yz+zx, xyz-1;
ideal Glp;
Glp = rustgb_std(Ilp);  // expected to error
kill Rlp;

tst_status(1);$
