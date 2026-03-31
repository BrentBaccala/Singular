LIB "tst.lib";
tst_init();

// dump of list(1) was list()
list T0;
list T1=1;
list T2=1,2;
list T3="a","b",3;
dump(":w dump-list");
read("dump-list");
getdump("dump-list");
T0;T1;
system("sh","rm -f dump-list");

tst_status(1);$
