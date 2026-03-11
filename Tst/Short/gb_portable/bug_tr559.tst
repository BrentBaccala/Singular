// Modified from Short/bug_tr559.tst
// Changed to check invariants instead of raw GB output.
// Tests reduce/std over qring with ring-cf.

LIB "tst.lib";
tst_init();

// reduce/std over qring with ring-cf:
ring rng = integer,(x,y,z),dp;
ideal iq = 2*y+1, 4*x*z+3*y, y^2-2*x*z-y;
iq = std(iq);

qring rngQ = iq;

ideal J = 3, y-4, x*z;

ideal stdJ = groebner(J);
// Check that J reduces to 0 mod its own GB
size(reduce(J, stdJ));
// Check GB size
size(stdJ);

tst_status(1);$
