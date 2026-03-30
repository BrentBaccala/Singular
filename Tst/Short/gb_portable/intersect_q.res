STDIN   1> // Modified from Short/intersect_q.tst
STDIN   2. // Removed option(prot) to make output algorithm-independent.
STDIN   3. // Oscar #4249: intersect in quotient ring.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> USER    :claude
init >> HOSTNAME:samsung
init >> uname -a:Linux samsung 6.8.1-1039-realtime #40-Ubuntu SMP PREEMPT_RT Mon Nov 24 21:35:45 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux
init >> date    :Tue Mar 10 08:55:42 PM EDT 2026
init >> version :44105
init >> ticks   :100
init >> 
STDIN   7> 
STDIN   8. ring R=QQ,(x,y),dp;
STDIN   9> qring A=std(ideal(x^2-y^3, x-y));
STDIN  10> 
STDIN  11. intersect(ideal(y2),ideal(x));
_[1]=xy
STDIN  12> 
STDIN  13. tst_status(1);$
1 >> tst_memory_0 :: samsung:86264
1 >> tst_memory_1 :: samsung:2207744
1 >> tst_memory_2 :: samsung:2207744
1 >> tst_timer_1 :: samsung:18
