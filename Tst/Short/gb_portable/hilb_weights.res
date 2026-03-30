STDIN   1> // Modified from Short/hilb_weights.tst
STDIN   2. // Removed option(prot) to make output algorithm-independent.
STDIN   3. // The test checks that hilbert-driven std produces the same result as plain std.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> USER    :claude
init >> HOSTNAME:samsung
init >> uname -a:Linux samsung 6.8.1-1039-realtime #40-Ubuntu SMP PREEMPT_RT Mon Nov 24 21:35:45 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux
init >> date    :Tue Mar 10 08:55:41 PM EDT 2026
init >> version :44105
init >> ticks   :100
init >> 
STDIN   7> 
STDIN   8. intvec w=23,11,11*23,11*23;
STDIN   9> ring s2=32003,(t,x,y,z),(wp(w));
STDIN  10> ideal i=
STDIN  11. 9x8+y7t+5x4y2t2+2xy2z3t2,
STDIN  12. 9y8+7xy6t+2x5yt2+2x2yz3t2,
STDIN  13. 9z8+3x2y2z2t2;
STDIN  14> ideal j = subst(i,t, t^11);
STDIN  15> j=subst(j,x,x23);
STDIN  16> ideal std_j = std(j);
STDIN  17> bigintvec v_j = hilb(std_j, 1, w);
STDIN  18> ideal k=std(j,v_j,w);
STDIN  19> 
STDIN  20. matrix(k) == matrix(std_j);
1
STDIN  21> 
STDIN  22. tst_status(1);$
1 >> tst_memory_0 :: samsung:347240
1 >> tst_memory_1 :: samsung:2383872
1 >> tst_memory_2 :: samsung:2408448
1 >> tst_timer_1 :: samsung:19
