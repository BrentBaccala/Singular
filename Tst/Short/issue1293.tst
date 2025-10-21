LIB "tst.lib";
tst_init();

ring sage12=(ZZ,27),(x, y),dp;
poly sage13=3;
def sage14=3;
poly sage15=sage13 ^ sage14;
poly sage16=0;
sage15 == sage16; // should return 1
// binomials:
sage13=3+x;
sage15=sage13 ^ sage14;
sage16=x3+9x2;
sage15==sage16;
// more terms
sage13=3+x+y;
sage15=sage13 ^ sage14;
sage16=x3+3x2y+3xy2+y3+9x2+18xy+9y2;
sage15==sage16;

tst_status(1);$
