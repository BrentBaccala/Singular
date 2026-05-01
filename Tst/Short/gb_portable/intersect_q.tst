// Modified from Short/intersect_q.tst
// Removed option(prot) to make output algorithm-independent.
// Oscar #4249: intersect in quotient ring.

LIB "tst.lib";
tst_init();

ring R=QQ,(x,y),dp;
qring A=std(ideal(x^2-y^3, x-y));

intersect(ideal(y2),ideal(x));

tst_status(1);$
