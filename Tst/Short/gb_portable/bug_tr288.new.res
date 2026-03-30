STDIN   1> // Modified from Short/bug_tr288.tst
STDIN   2. // Changed to check invariants instead of raw GB output.
STDIN   3. // Tests cancelunit for 1+y (ecart==0) with a(1,-1) ordering.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> bug_tr288.new.stat
STDIN   7> 
STDIN   8. // missed cancelunit for 1+y (ecart==0)
STDIN   9. 
STDIN  10. ring r = 0, (x, y), (a(1, -1), dp); x+y+1;
x+1+y
STDIN  11> size(std(ideal(1+y)));
1
STDIN  12> size(std(ideal(1+y+x)));
1
STDIN  13> option(redSB);
STDIN  14> size(std(ideal(1+y)));
1
STDIN  15> size(std(ideal(1+y+x)));
1
STDIN  16> 
STDIN  17. tst_status(1);$
