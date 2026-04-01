// bug_tr380 - invariant-checking version
// Original prints division result; this checks the division identity

LIB "tst.lib";
tst_init();

// error from idLift not tested for
ring r=(16,a),(e2,e1,y(2..1)),(wp(1+30,1,1,1));
qring rq=std(ideal(e1^2, e2^2, e1*e2));
list L=division(e2, std(0));
// division should return a list of size 3
size(L);
// check the division identity: e2 * L[3] == std(0) * L[1] + L[2]
// since we divide by std(0), L[2] should equal e2 * L[3]
matrix(e2)*matrix(L[3])==matrix(L[2]);

tst_status(1);$
