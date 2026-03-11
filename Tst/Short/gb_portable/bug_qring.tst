// Modified from Short/bug_qring.tst
// Changed to check invariants instead of raw GB output.
// Tests normal form and complete NF in qrings.

LIB "tst.lib";
tst_init();

// normal form and complete NF in qrings:

ring p=0,(y,z,u,v),dp;
option(redSB);
qring q=std(ideal(y+u2+uv3,z+uv3));
q;
ideal i=y,z;
ideal si = std(i);
size(si);
dim(si);
// NF and interred should reduce to equivalent results
size(NF(si,0));
size(interred(si));
setring p;
ideal i=y,z,y+u2+uv3,z+uv3;
option(redSB);
ideal si = std(i);
size(si);
dim(si);
size(NF(si,std(ideal(y+u2+uv3,z+uv3))));
size(interred(si));
size(interred(NF(si,std(ideal(y+u2+uv3,z+uv3)))));

ring r = 2,(a, b),dp;
option(redSB);
size(interred(ideal(a^2,a*b+a^2)));
qring qq = std(ideal(b^2+a*b+a^2,a^3));
qq;
option(redSB);
size(interred(ideal(a^2,a*b+a^2)));

tst_status(1);$
