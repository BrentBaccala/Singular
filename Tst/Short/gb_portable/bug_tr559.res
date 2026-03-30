STDIN   1> // Modified from Short/bug_tr559.tst
STDIN   2. // Changed to check invariants instead of raw GB output.
STDIN   3. // Tests reduce/std over qring with ring-cf.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> USER    :claude
init >> HOSTNAME:samsung
init >> uname -a:Linux samsung 6.8.1-1039-realtime #40-Ubuntu SMP PREEMPT_RT Mon Nov 24 21:35:45 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux
init >> date    :Tue Mar 10 08:55:40 PM EDT 2026
init >> version :44105
init >> ticks   :100
init >> 
STDIN   7> 
STDIN   8. // reduce/std over qring with ring-cf:
STDIN   9. ring rng = integer,(x,y,z),dp;
STDIN  10> ideal iq = 2*y+1, 4*x*z+3*y, y^2-2*x*z-y;
STDIN  11> iq = std(iq);
STDIN  12> 
STDIN  13. qring rngQ = iq;
STDIN  14> 
STDIN  15. ideal J = 3, y-4, x*z;
STDIN  16> 
STDIN  17. ideal stdJ = groebner(J);
STDIN  18> // Check that J reduces to 0 mod its own GB
STDIN  19. size(reduce(J, stdJ));
0
STDIN  20> // Check GB size
STDIN  21. size(stdJ);
3
STDIN  22> 
STDIN  23. tst_status(1);$
1 >> tst_memory_0 :: samsung:101472
1 >> tst_memory_1 :: samsung:2207744
1 >> tst_memory_2 :: samsung:2248704
1 >> tst_timer_1 :: samsung:18
