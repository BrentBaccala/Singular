STDIN   1> // Modified from Short/bug_orderM.tst
STDIN   2. // Changed to check invariants instead of raw GB output.
STDIN   3. // Tests std with matrix ordering.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> bug_orderM.new.stat
STDIN   7> 
STDIN   8. ring r=0,(x,y,z),(M(1,1,1,0,1,0,1,0,0));
STDIN   9> ideal i=x + y + z, x^2 + y^2 + z^3;
STDIN  10> ideal g = std(i);
STDIN  11> // Check generators reduce to 0 in GB
STDIN  12. size(reduce(i, g));
0
STDIN  13> // Check GB size and dimension
STDIN  14. size(g);
2
STDIN  15> dim(g);
1
STDIN  16> 
STDIN  17. tst_status(1);$
