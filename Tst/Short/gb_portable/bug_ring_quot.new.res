STDIN   1> // Modified from Short/bug_ring_quot.tst
STDIN   2. // Changed to check invariants instead of raw GB output.
STDIN   3. // Tests enterOnePairSpecial and quotient over integer ring.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> bug_ring_quot.new.stat
STDIN   7> 
STDIN   8. // wrong enterOnePairSpecial
STDIN   9. ring R=integer,(x,y,z),dp;
STDIN  10> ideal I=6xz+yz,xz+yz+z;
STDIN  11> ideal sI = std(I);
STDIN  12> // Check I reduces to 0 in its own GB
STDIN  13. size(reduce(I, sI));
0
STDIN  14> // Check size
STDIN  15. size(sI);
2
STDIN  16> // Check quotient
STDIN  17. ideal Q = quotient(I,5);
STDIN  18> size(Q);
2
STDIN  19> dim(std(Q));
3
STDIN  20> 
STDIN  21. // tr. # 490
STDIN  22. ring RR=integer,(x,y),dp;
STDIN  23> module N=[x,1],[x,y];
STDIN  24> module M=freemodule(2);
STDIN  25> ideal Q2 = quotient(N,M);
STDIN  26> size(Q2);
1
STDIN  27> 
STDIN  28. tst_status(1);$
