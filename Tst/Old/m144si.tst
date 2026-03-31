  monitor("doe.tmp");
  ring r;
  poly f=x+y+z;
  int i=7;
  ideal I=f,x,y;
  monitor("");
  system("sh","rm -f doe.tmp");
LIB "tst.lib";tst_status(1);$
