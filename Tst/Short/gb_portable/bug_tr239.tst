// bug_tr239 - invariant-checking version
// Original prints finduni() results; this checks size and membership

LIB "tst.lib";
tst_init();

// finduni required id, should now also take data:
ring  r=0,(x,y,z), dp;
ideal i=y3+x2,x2y+x2,z4-x2-y;
option(redSB);
i=std(i);
ideal k=finduni(i);
// should find univariate polynomials
size(k);
// all should be in the ideal
size(reduce(k, i));
// calling with std(i) directly should give same size
size(finduni(std(i)));

tst_status(1);$
