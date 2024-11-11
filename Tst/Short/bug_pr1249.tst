LIB "tst.lib";
tst_init();

option(prot);
printlevel=16;

ring r = 0,(x,y,z),dp;
map f = r,x,xyz,z;
poly p = y;
ideal i = y;
f(i);

tst_status(1);$
