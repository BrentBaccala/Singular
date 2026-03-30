STDIN   1> // Modified from Short/bug_tr132.tst
STDIN   2. // The reduce and std results should be algorithm-independent since they
STDIN   3. // operate in a quotient ring. Output unchanged since it checks reduce results.
STDIN   4. 
STDIN   5. LIB "tst.lib";
STDIN   6> tst_init();
init >> bug_tr132.new.stat
STDIN   7> 
STDIN   8. // reduce was not complete (strat->ak was 0)
STDIN   9. ring r = 0,(x,y,z), (c, dp);
STDIN  10> qring Q = std(ideal(var(1)**2, var(2)**2,
STDIN  11.  var(3)**2));
STDIN  12> reduce( maxideal(2) * gen(1), std(0));
_[1]=0
_[2]=[yz]
_[3]=0
_[4]=[xz]
_[5]=[xy]
_[6]=0
STDIN  13> 
STDIN  14. // std (compleReduce) was not complete (index bounds to small)
STDIN  15. vector v = var(1)**2 + var(2)**2 + var(1)*var(2); v;
[x2+xy+y2]
STDIN  16> option(redTail); option(redSB);
STDIN  17> std(v);
_[1]=[xy]
STDIN  18> 
STDIN  19. tst_status(1);$
