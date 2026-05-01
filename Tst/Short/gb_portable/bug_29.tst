// bug_29 - invariant-checking version
// Original prints resolution; this checks betti numbers

LIB "tst.lib";
tst_init();

ring r=0,(x,y),dp;
ideal i=xy;
qring rq=std(i);
ideal j=x;
resolution mr=mres(j,0);
// check resolution length
size(list(mr));
// check betti numbers (invariant)
print(betti(mr), "betti");

tst_status(1);$
