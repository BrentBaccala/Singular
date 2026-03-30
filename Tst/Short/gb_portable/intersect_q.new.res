STDIN   1> // Modified from Short/intersect_q.tst
STDIN   2. // Removed option(prot) to make output algorithm-independent.
STDIN   3. // Oscar #4249: intersect in quotient ring.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> intersect_q.new.stat
STDIN   7> 
STDIN   8. ring R=QQ,(x,y),dp;
STDIN   9> qring A=std(ideal(x^2-y^3, x-y));
STDIN  10> 
STDIN  11. intersect(ideal(y2),ideal(x));
_[1]=xy
STDIN  12> 
STDIN  13. tst_status(1);$
