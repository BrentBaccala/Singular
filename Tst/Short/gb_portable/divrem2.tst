// divrem2 - invariant-checking version
// Original prints remainder; this checks the division identity

LIB "tst.lib";
tst_init();

ring r=QQ,(x,y),dp;
ideal a=x^3+y^3+x*y;
ideal b=x;
list L=system("DivRemIdU",a,b);
// division identity: a * u == b * q + r
matrix(a)*matrix(L[3])==matrix(b)*matrix(L[2])+matrix(L[1]);
// remainder should match NF
ideal rest=NF(a,std(b));
matrix(rest)==matrix(L[1]);

tst_status(1);$
