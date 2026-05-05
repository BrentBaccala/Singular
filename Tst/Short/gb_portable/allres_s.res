STDIN   1> // allres_s - invariant-checking version
STDIN   2. // Original prints betti tables from tst_test_res; this checks return values
STDIN   3. 
STDIN   4. LIB "tst.lib";
STDIN   5> tst_init();
init >> allres_s.new.stat
STDIN   6> 
STDIN   7. ring an=32003,(w,x,y,z),(dp,C);
STDIN   8> ideal i=
STDIN   9. 1w2xy+1w2xz+1w2yz+1wxyz+1x2yz+1xy2z+1xyz2,
STDIN  10. 1w4x+1w4z+1w3yz+1w2xyz+1wx2yz+1x2y2z+1xy2z2,
STDIN  11. 1w6+1w5z+1w4xz+1w3xyz+1w2xy2z+1wx2y2z+1x2y2z2;
STDIN  12> tst_test_res(i);
1
STDIN  13> kill an;
STDIN  14> 
STDIN  15. ring an=32003,(w,x,y,z),(dp,c);
STDIN  16> ideal i=
STDIN  17. wx2+y3,
STDIN  18. xy2+z3,
STDIN  19. yz2+w3,
STDIN  20. zw2+x3,
STDIN  21. xyz+yzw+zwx+wxy;
STDIN  22> tst_test_res(i);
1
STDIN  23> kill an;
STDIN  24> 
STDIN  25. ring an=32003,(w,x,y,z),(dp,c);
STDIN  26> ideal i=
STDIN  27. wx+y2,
STDIN  28. xy+z2,
STDIN  29. yz+w2,
STDIN  30. zw+x2,
STDIN  31. xy+yz+zw+wx;
STDIN  32> tst_test_res(i);
1
STDIN  33> kill an;
STDIN  34> 
STDIN  35. ring an=32003,(x,y,z),(dp,c);
STDIN  36> ideal i=
STDIN  37. zx2+y3,
STDIN  38. xy2+z3,
STDIN  39. yz2+x3,
STDIN  40. x2y+y2z+z2x;
STDIN  41> tst_test_res(i);
1
STDIN  42> kill an;
STDIN  43> 
STDIN  44. ring an=0,(w,x,y,z),(dp,C);
STDIN  45> ideal i=
STDIN  46. 1w2xy+1w2xz+1w2yz+1wxyz+1x2yz+1xy2z+1xyz2,
STDIN  47. 1w4x+1w4z+1w3yz+1w2xyz+1wx2yz+1x2y2z+1xy2z2,
STDIN  48. 1w6+1w5z+1w4xz+1w3xyz+1w2xy2z+1wx2y2z+1x2y2z2;
STDIN  49> tst_test_res(i);
1
STDIN  50> kill an;
STDIN  51> 
STDIN  52. ring an=0,(w,x,y,z),(dp,c);
STDIN  53> ideal i=
STDIN  54. wx2+y3,
STDIN  55. xy2+z3,
STDIN  56. yz2+w3,
STDIN  57. zw2+x3,
STDIN  58. xyz+yzw+zwx+wxy;
STDIN  59> tst_test_res(i, 1);
1
STDIN  60> kill an;
STDIN  61> 
STDIN  62. ring an=0,(w,x,y,z),(dp,c);
STDIN  63> ideal i=
STDIN  64. wx+y2,
STDIN  65. xy+z2,
STDIN  66. yz+w2,
STDIN  67. zw+x2,
STDIN  68. xy+yz+zw+wx;
STDIN  69> tst_test_res(i);
1
STDIN  70> kill an;
STDIN  71> 
STDIN  72. ring an=0,(x,y,z),(dp,c);
STDIN  73> ideal i=
STDIN  74. zx2+y3,
STDIN  75. xy2+z3,
STDIN  76. yz2+x3,
STDIN  77. x2y+y2z+z2x;
STDIN  78> tst_test_res(i);
1
STDIN  79> kill an;
STDIN  80> 
STDIN  81. tst_status(1);$
