// Modified from Short/bug_tr288.tst
// Changed to check invariants instead of raw GB output.
// Tests cancelunit for 1+y (ecart==0) with a(1,-1) ordering.

LIB "tst.lib";
tst_init();

// missed cancelunit for 1+y (ecart==0)

ring r = 0, (x, y), (a(1, -1), dp); x+y+1;
size(std(ideal(1+y)));
size(std(ideal(1+y+x)));
option(redSB);
size(std(ideal(1+y)));
size(std(ideal(1+y+x)));

tst_status(1);$
