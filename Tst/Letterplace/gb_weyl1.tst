LIB "tst.lib"; tst_init();
LIB "freegb.lib";
ring r = 0,(d,x),Dp;
int upToDeg = 4;
def R = makeLetterplaceRing(upToDeg);
setring(R);
ideal Id = d*x-x*d-1;
option(redTail);
option(redSB);
std(Id);
tst_status(1);$
