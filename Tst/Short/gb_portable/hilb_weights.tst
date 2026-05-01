// Modified from Short/hilb_weights.tst
// Removed option(prot) to make output algorithm-independent.
// The test checks that hilbert-driven std produces the same result as plain std.

LIB "tst.lib";
tst_init();

intvec w=23,11,11*23,11*23;
ring s2=32003,(t,x,y,z),(wp(w));
ideal i=
9x8+y7t+5x4y2t2+2xy2z3t2,
9y8+7xy6t+2x5yt2+2x2yz3t2,
9z8+3x2y2z2t2;
ideal j = subst(i,t, t^11);
j=subst(j,x,x23);
ideal std_j = std(j);
bigintvec v_j = hilb(std_j, 1, w);
ideal k=std(j,v_j,w);

matrix(k) == matrix(std_j);

tst_status(1);$
