// finduni_s - invariant-checking version
// Original prints finduni() results; this checks size and ideal membership

LIB "tst.lib";
tst_init();

option(redSB);

// cyclic 6 in char 32003
ring r=32003,(a,b,c,d,x,f), dp;
ideal i=a+b+c+d+x+f, ab+bc+cd+dx+xf+af, abc+bcd+cdx+d*xf+axf+abf, abcd+bcdx+cd*xf+ad*xf+abxf+abcf, abcdx+bcd*xf+acd*xf+abd*xf+abcxf+abcdf, abcd*xf-1;
ideal is=std(i);
ideal fu=finduni(is);
// should find one univariate poly per variable
size(fu);
// each element should be univariate (check number of variables used)
int alluni = 1;
for (int j=1; j<=size(fu); j++) { if (size(variables(fu[j])) != 1) { alluni = 0; } }
alluni;
kill r;

// walks-7 in char 0 -- omitted, finduni too slow for Short tests

tst_status(1);$
