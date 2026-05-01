// bug_24 - invariant-checking version
// Original prints ideals; this checks preimage round-trip via reduce

LIB "tst.lib";
tst_init();

ring r=0,(x,y,z),ds;
map m=r,x+y+z+y2+x3,x-y-z,z;
ideal id=2y+2z+x^3+y^2,x-y-z,y2+2z3;
ideal pid=preimage(basering,m,id);
ideal imid=m(pid);
// id and m(preimage(m, id)) should generate the same ideal
size(reduce(id,std(imid)));
size(reduce(imid,std(id)));

tst_status(1);$
