STDIN   1> // Modified from Short/bug_orderM.tst
STDIN   2. // Changed to check invariants instead of raw GB output.
STDIN   3. // Tests std with matrix ordering.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> USER    :claude
init >> HOSTNAME:samsung
init >> uname -a:Linux samsung 6.8.1-1039-realtime #40-Ubuntu SMP PREEMPT_RT Mon Nov 24 21:35:45 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux
init >> date    :Tue Mar 10 08:55:39 PM EDT 2026
init >> version :44105
init >> ticks   :100
init >> 
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
1 >> tst_memory_0 :: samsung:86472
1 >> tst_memory_1 :: samsung:2207744
1 >> tst_memory_2 :: samsung:2232320
1 >> tst_timer_1 :: samsung:17
