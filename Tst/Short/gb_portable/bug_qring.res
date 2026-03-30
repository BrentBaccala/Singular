STDIN   1> // Modified from Short/bug_qring.tst
STDIN   2. // Changed to check invariants instead of raw GB output.
STDIN   3. // Tests normal form and complete NF in qrings.
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
STDIN   8. // normal form and complete NF in qrings:
STDIN   9. 
STDIN  10. ring p=0,(y,z,u,v),dp;
STDIN  11> option(redSB);
STDIN  12> qring q=std(ideal(y+u2+uv3,z+uv3));
STDIN  13> q;
// coefficients: QQ considered as a field
// number of vars : 4
//        block   1 : ordering dp
//                  : names    y z u v
//        block   2 : ordering C
// quotient ring from ideal
_[1]=u2+y-z
_[2]=uv3+z
_[3]=yv3-zv3-zu
STDIN  14> ideal i=y,z;
STDIN  15> ideal si = std(i);
STDIN  16> size(si);
2
STDIN  17> dim(si);
1
STDIN  18> // NF and interred should reduce to equivalent results
STDIN  19. size(NF(si,0));
// ** _ is no standard basis
2
STDIN  20> size(interred(si));
2
STDIN  21> setring p;
STDIN  22> ideal i=y,z,y+u2+uv3,z+uv3;
STDIN  23> option(redSB);
STDIN  24> ideal si = std(i);
STDIN  25> size(si);
4
STDIN  26> dim(si);
1
STDIN  27> size(NF(si,std(ideal(y+u2+uv3,z+uv3))));
4
STDIN  28> size(interred(si));
4
STDIN  29> size(interred(NF(si,std(ideal(y+u2+uv3,z+uv3)))));
2
STDIN  30> 
STDIN  31. ring r = 2,(a, b),dp;
STDIN  32> option(redSB);
STDIN  33> size(interred(ideal(a^2,a*b+a^2)));
2
STDIN  34> qring qq = std(ideal(b^2+a*b+a^2,a^3));
STDIN  35> qq;
// coefficients: ZZ/2 considered as a field
// number of vars : 2
//        block   1 : ordering dp
//                  : names    a b
//        block   2 : ordering C
// quotient ring from ideal
_[1]=a2+ab+b2
_[2]=b3
STDIN  36> option(redSB);
STDIN  37> size(interred(ideal(a^2,a*b+a^2)));
2
STDIN  38> 
STDIN  39. tst_status(1);$
1 >> tst_memory_0 :: samsung:90944
1 >> tst_memory_1 :: samsung:2207744
1 >> tst_memory_2 :: samsung:2207744
1 >> tst_timer_1 :: samsung:17
