STDIN   1> // Modified from Short/bug_tr159.tst
STDIN   2. // Removed option(prot) to make output algorithm-independent.
STDIN   3. // The test checks that NF is 0 as expected (bug was in kFind...InS).
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> bug_tr159.new.stat
STDIN   7> 
STDIN   8. // NF was not 0 as expected - not enough searching in kFind...InS
STDIN   9. intmat m1[8][8] =
STDIN  10.  2, 2, 2, 4, 1, 1, 3, 3,
STDIN  11.  0, 0, 0,-1, 0, 0, 0, 0,
STDIN  12.  0, 0, 0, 0,-1, 0, 0, 0,
STDIN  13. -1, 0, 0, 0, 0, 0, 0, 0,
STDIN  14.  0,-1, 0, 0, 0, 0, 0, 0,
STDIN  15.  0, 0,-1, 0, 0, 0, 0, 0,
STDIN  16.  0, 0, 0, 0, 0,-1, 0, 0,
STDIN  17.  0, 0, 0, 0, 0, 0,-1, 0;
STDIN  18> ring R1 = 2, (b_2_1,b_2_2,b_2_3,c_4_8,a_1_0,b_1_1,b_3_4,b_3_5), M(m1);
STDIN  19> ideal I;
STDIN  20> I[1]=a_1_0^2;
STDIN  21> I[2]=a_1_0*b_1_1;
STDIN  22> I[3]=b_2_3*a_1_0+b_2_2*a_1_0;
STDIN  23> I[4]=b_2_1*b_1_1+b_2_2*a_1_0;
STDIN  24> I[5]=b_2_2*b_1_1+b_2_2*a_1_0;
STDIN  25> I[6]=a_1_0*b_3_4;
STDIN  26> I[7]=a_1_0*b_3_5;
STDIN  27> I[8]=b_2_2^2+b_2_1*b_2_3;
STDIN  28> I[9]=b_1_1*b_3_4;
STDIN  29> I[10]=b_2_2*b_3_4+b_2_1*b_3_5+b_2_1*b_2_2*a_1_0;
STDIN  30> I[11]=b_2_3*b_3_4+b_2_2*b_3_5+b_2_1*b_2_2*a_1_0;
STDIN  31> I[12]=b_3_4^2+b_2_1*b_2_3^2+b_2_1^2*b_2_3;
STDIN  32> I[13]=b_3_4*b_3_5+b_2_2*b_2_3^2+b_2_1*b_2_2*b_2_3;
STDIN  33> I[14]=b_3_5^2+b_2_3*b_1_1*b_3_5+b_2_3^3+b_2_1*b_2_3^2+c_4_8*b_1_1^2;
STDIN  34> qring Q1 = groebner(I);
STDIN  35> 
STDIN  36. 
STDIN  37. intmat m2[5][5] =
STDIN  38.  2, 2, 2, 1, 1,
STDIN  39.  0,-1,-1, 0, 0,
STDIN  40.  0, 0, 0,-1, 0,
STDIN  41. -1, 0, 0, 0, 0,
STDIN  42.  0,-1, 0, 0, 0;
STDIN  43> ring R2 = 2, (@b_2_1,@c_2_2,@c_2_3,@a_1_0,@b_1_1), M(m2);
STDIN  44> ideal I;
STDIN  45> I[1]=@a_1_0^2;
STDIN  46> I[2]=@a_1_0*@b_1_1;
STDIN  47> I[3]=@b_2_1*@a_1_0;
STDIN  48> I[4]=@b_2_1^2+@c_2_2*@b_1_1^2;
STDIN  49> qring Q2 = groebner(I);
STDIN  50> 
STDIN  51. def Q = Q1+Q2;
STDIN  52> setring Q;
STDIN  53> Q;
// coefficients: ZZ/2 considered as a field
// number of vars : 13
//        block   1 : ordering M
//                  : names    b_2_1 b_2_2 b_2_3 c_4_8 a_1_0 b_1_1 b_3_4 b_3_5
//                  : weights      2     2     2     4     1     1     3     3
//                  : weights      0     0     0    -1     0     0     0     0
//                  : weights      0     0     0     0    -1     0     0     0
//                  : weights     -1     0     0     0     0     0     0     0
//                  : weights      0    -1     0     0     0     0     0     0
//                  : weights      0     0    -1     0     0     0     0     0
//                  : weights      0     0     0     0     0    -1     0     0
//                  : weights      0     0     0     0     0     0    -1     0
//        block   2 : ordering M
//                  : names    @b_2_1 @c_2_2 @c_2_3 @a_1_0 @b_1_1
//                  : weights       2      2      2      1      1
//                  : weights       0     -1     -1      0      0
//                  : weights       0      0      0     -1      0
//                  : weights      -1      0      0      0      0
//                  : weights       0     -1      0      0      0
//        block   3 : ordering C
// quotient ring from ideal
_[1]=a_1_0^2
_[2]=a_1_0*b_1_1
_[3]=b_2_3*a_1_0+b_2_2*a_1_0
_[4]=b_2_1*b_1_1+b_2_2*a_1_0
_[5]=b_2_2*b_1_1+b_2_2*a_1_0
_[6]=a_1_0*b_3_4
_[7]=a_1_0*b_3_5
_[8]=b_2_2^2+b_2_1*b_2_3
_[9]=b_1_1*b_3_4
_[10]=b_2_2*b_3_4+b_2_1*b_3_5+b_2_1*b_2_2*a_1_0
_[11]=b_2_3*b_3_4+b_2_2*b_3_5+b_2_1*b_2_2*a_1_0
_[12]=b_3_4^2+b_2_1*b_2_3^2+b_2_1^2*b_2_3
_[13]=b_3_4*b_3_5+b_2_2*b_2_3^2+b_2_1*b_2_2*b_2_3
_[14]=b_3_5^2+b_2_3*b_1_1*b_3_5+b_2_3^3+b_2_1*b_2_3^2+c_4_8*b_1_1^2
_[15]=@a_1_0^2
_[16]=@a_1_0*@b_1_1
_[17]=@b_2_1*@a_1_0
_[18]=@b_2_1^2+@c_2_2*@b_1_1^2
STDIN  54> 
STDIN  55. ideal t = a_1_0, b_1_1, b_2_1, b_2_2, b_2_3, b_3_4, b_3_5, c_4_8;
STDIN  56> NF(b_2_1,t);
// ** t is no standard basis
0
STDIN  57> 
STDIN  58. tst_status(1);$
