// Modified from Short/bug_ring_quot.tst
// Changed to check invariants instead of raw GB output.
// Tests enterOnePairSpecial and quotient over integer ring.

LIB "tst.lib";
tst_init();

// wrong enterOnePairSpecial
ring R=integer,(x,y,z),dp;
ideal I=6xz+yz,xz+yz+z;
ideal sI = std(I);
// Check I reduces to 0 in its own GB
size(reduce(I, sI));
// Check size
size(sI);
// Check quotient
ideal Q = quotient(I,5);
size(Q);
dim(std(Q));

// tr. # 490
ring RR=integer,(x,y),dp;
module N=[x,1],[x,y];
module M=freemodule(2);
ideal Q2 = quotient(N,M);
size(Q2);

tst_status(1);$
