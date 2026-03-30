STDIN   1> // Modified from Short/bug_tr874.tst
STDIN   2. // Output unchanged — codim() is an invariant check, already algorithm-independent.
STDIN   3. 
STDIN   4. LIB "tst.lib";
STDIN   5> tst_init();
init >> bug_tr874.new.stat
STDIN   6> 
STDIN   7. LIB"sing.lib";
STDIN   8> ring r = 0,(x,y),dp;
STDIN   9> ideal i1 = x;
STDIN  10> ideal i2 = x2, xy;
STDIN  11> codim(std(i1),std(i2));
1
STDIN  12> 
STDIN  13. tst_status(1);$
