// Modified from Short/bug_tr873.tst
// Changed to check size instead of raw GB output.
// Tests ring creation in a loop with qring coefficients.

LIB "tst.lib";
tst_init();

ring Q1=0,i,lp;
qring Q2=i;
for(int i=1;i<5;i++)
{
ring F = 0,(a,b),dp;
ring R = F,(c),dp;
ideal I = c;
size(std(I));
ring RR = Q2,(c),dp;
ideal I = c;
size(std(I));

}

tst_status(1);$
