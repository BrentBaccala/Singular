// Modified from Short/bug_orderM.tst
// Changed to check invariants instead of raw GB output.
// Tests std with matrix ordering.

LIB "tst.lib";
tst_init();

ring r=0,(x,y,z),(M(1,1,1,0,1,0,1,0,0));
ideal i=x + y + z, x^2 + y^2 + z^3;
ideal g = std(i);
// Check generators reduce to 0 in GB
size(reduce(i, g));
// Check GB size and dimension
size(g);
dim(g);

tst_status(1);$
