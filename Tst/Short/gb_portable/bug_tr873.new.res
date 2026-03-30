STDIN   1> // Modified from Short/bug_tr873.tst
STDIN   2. // Changed to check size instead of raw GB output.
STDIN   3. // Tests ring creation in a loop with qring coefficients.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> bug_tr873.new.stat
STDIN   7> 
STDIN   8. ring Q1=0,i,lp;
STDIN   9> qring Q2=i;
STDIN  10> for(int i=1;i<5;i++)
STDIN  11. {
STDIN  12. ring F = 0,(a,b),dp;
STDIN  14. ring R = F,(c),dp;
STDIN  16. ideal I = c;
STDIN  18. size(std(I));
STDIN  20. ring RR = Q2,(c),dp;
STDIN  22. ideal I = c;
STDIN  24. size(std(I));
STDIN  26. 
STDIN  28. }
1
1
// ** redefining F (ring F = 0,(a,b),dp;)
// ** redefining R (ring R = F,(c),dp;)
1
// ** redefining RR (ring RR = Q2,(c),dp;)
1
// ** redefining F (ring F = 0,(a,b),dp;)
// ** redefining R (ring R = F,(c),dp;)
1
// ** redefining RR (ring RR = Q2,(c),dp;)
1
// ** redefining F (ring F = 0,(a,b),dp;)
// ** redefining R (ring R = F,(c),dp;)
1
// ** redefining RR (ring RR = Q2,(c),dp;)
1
STDIN  29> 
STDIN  30. tst_status(1);$
