STDIN   1> option(prot,sugarCrit);
STDIN   2> " ============= cyclic_roots_5(isol) + ==========================";
 ============= cyclic_roots_5(isol) + ==========================
STDIN   3> ring r4 = 0,(a,b,c,d,e),dp;
STDIN   4> r4;
// coefficients: QQ considered as a field
// number of vars : 5
//        block   1 : ordering dp
//                  : names    a b c d e
//        block   2 : ordering C
STDIN   5> poly s1=a+b+c+d+e;
STDIN   6> poly s2=de+1cd+1bc+1ae+1ab;
STDIN   7> poly s3=cde+1bcd+1ade+1abe+1abc;
STDIN   8> poly s4=bcde+1acde+1abde+1abce+1abcd;
STDIN   9> poly s5=abcde+1;
STDIN  10> ideal i=s1,s2,s3,s4,s5;
STDIN  11> ideal j=std(i);
[1048575:3]1(4)s2(3)s3(2)s4ss5(3)s(5)s(6)s(8)6-ss(9)s(10)s(12)s(14)s(16)-7-s(15)s(17)s(19)s(21)s(23)----s(21)----8-s(18)s(21)s(23)s(25)s(27)s(30)----------9-s(22)s(25)s(28)s(31)------------10-----------11-s(10)s(13)s(16)--12-ss(19)---------13-s(11)s(13)s(16)------14---------15-
product criterion:81 chain criterion:364
STDIN  12> "dim: "+ string(dim(j)) +",  mult: "+ string(mult(j)) +",  elem: "
STDIN  13.                                                         + string(size(j));
dim: 0,  mult: 70,  elem: 20
STDIN  14> kill r4;
STDIN  15> "====================  standard char0  =============================";
====================  standard char0  =============================
STDIN  16> //      H7 l, char 0, test0,11,1: 61/31 ohne vollst; Reduktion
STDIN  17. ring r=
STDIN  18. 0,(x,y),lp;
STDIN  19> poly f=x5+y11+xy9+x3y9;
STDIN  20> ideal i=jacob(f);
STDIN  21> ideal j=std( i);
std in char. 32003, homogenized ------------------
[1048575:1]11ss12s14s
product criterion:0 chain criterion:3
stdhilb in basering, homogenized ------------------
[1048575:1]11ss19s20s21sh23sh27s28sh29sh30sh31sh32sh33sh34sh35sh36sh37sh38sh39sh40shh
product criterion:0 chain criterion:156
hilbert series criterion:16
de-homogenize, interred ------------------
STDIN  22> size(j);
3
STDIN  23> kill r;
STDIN  24> " ============= cyclic_roots_6(homog) dp ==========================";
 ============= cyclic_roots_6(homog) dp ==========================
STDIN  25> ring r6 = 0,(a,b,c,d,e,f),dp;
STDIN  26> poly s1=a+b+c+d+e+f;
STDIN  27> poly s2=ab+bc+cd+de+ef+fa;
STDIN  28> poly s3=abc+bcd+cde+edf+efa+fab;
STDIN  29> poly s4=abcd+bcde+cdef+defa+efab+fabc;
STDIN  30> poly s5=abcde+bcdef+cdefa+defab+efabc+fabcd;
STDIN  31> poly s6=abcdef+1;
STDIN  32> ideal i=s1,s2,s3,s4,s5,s6;
STDIN  33> ideal j=std(i);
[1048575:3]1(5)s2(4)s3(3)s4ss5(4)s(6)s(7)s(9)6-s(11)s(12)s(14)s(16)s(19)s(21)s(24)-7-s(23)s(24)s(27)s(29)s(31)s(32)s(35)-s(37)s(40)s(42)s(44)s(45)--s(46)s(48)-----8-s(44)s(47)s(50)s(52)s(55)s(57)s(59)s(61)-s(63)----s(62)----s(61)s(64)-s(66)-----------s(58)-------9-s(53)s(56)s(59)s(62)s(65)s(68)s(71)s(74)s(77)s(80)s(83)s(86)s(90)s(95)s(102)s(108)--------(100)----------------------s(81)---10-s(83)s(88)s(90)s(94)s(99)s(104)s(109)s(114)-s(116)s(121)s(126)s(128)s(132)--------------------------------(100)--------------11-s(87)---------------------------------------12-s(50)--------13-s(44)s(47)s(51)-s(54)------------14-s(45)s(48)s(51)s(55)s(58)s(61)s(64)s(67)s(70)--------------------15-s(52)s(55)s(58)s(61)s(64)s(67)s(70)s(73)s(76)s(79)s(82)-------------------------------------16--------------------------------------------17-
product criterion:233 chain criterion:2353
STDIN  34> "dim: "+ string(dim(j)) +",  mult: "+ string(mult(j)) +",  elem: "
STDIN  35.                                                         + string(size(j));
dim: 0,  mult: 156,  elem: 45
STDIN  36> kill r6;
STDIN  37> LIB "tst.lib";tst_status(1);$
