STDIN   1> // Modified from Short/bug_tr874.tst
STDIN   2. // Output unchanged — codim() is an invariant check, already algorithm-independent.
STDIN   3. 
STDIN   4. LIB "tst.lib";
STDIN   5> tst_init();
init >> USER    :claude
init >> HOSTNAME:samsung
init >> uname -a:Linux samsung 6.8.1-1039-realtime #40-Ubuntu SMP PREEMPT_RT Mon Nov 24 21:35:45 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux
init >> date    :Tue Mar 10 08:55:41 PM EDT 2026
init >> version :44105
init >> ticks   :100
init >> 
STDIN   6> 
STDIN   7. LIB"sing.lib";
STDIN   8> ring r = 0,(x,y),dp;
STDIN   9> ideal i1 = x;
STDIN  10> ideal i2 = x2, xy;
STDIN  11> codim(std(i1),std(i2));
1
STDIN  12> 
STDIN  13. tst_status(1);$
1 >> tst_memory_0 :: samsung:197384
1 >> tst_memory_1 :: samsung:2207744
1 >> tst_memory_2 :: samsung:2232320
1 >> tst_timer_1 :: samsung:19
