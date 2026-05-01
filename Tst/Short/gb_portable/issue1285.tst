// issue1285 - invariant-checking version
// Original prints reduce() result; this checks the reduction property

LIB "tst.lib";

ring R = (0,v1),z(1..3),ws(1..3);
ideal I = z(1)^2-z(2)+(-v1)*z(3)+(2*v1^2)*z(1)*z(3)+(-v1^3)*z(2)*z(3),
    z(2)^2-2*z(1)*z(3)+(v1)*z(2)*z(3),
    z(3)^2;
I = std(I);
poly p = 2*z(1)^2*z(3)+(v1)*z(1)*z(2)*z(3);
poly r = reduce(p, I);
// p - r should be in the ideal
size(reduce(p - r, I));
// r should be fully reduced (reducing again gives the same result)
reduce(r, I) == r;

tst_status(1);$
