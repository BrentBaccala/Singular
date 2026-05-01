// bug_division_2 - invariant-checking version
// Original prints matrix results; this checks the division identity

LIB "tst.lib";
tst_init();

ring r= 0,(x),ds;
matrix M[2][2]=[1,0,0,0];
matrix NN[2][1]=[1,0];
module N=std(NN);
list l=division(M,N);
matrix T=l[1];
matrix R=l[2];
matrix U = l[3];
// division identity: M * U == N * T + R
matrix(M)*U==matrix(N)*T+matrix(R);

tst_status(1);$
