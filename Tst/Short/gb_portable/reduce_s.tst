// reduce_s - invariant-checking version
// Original prints reduce() results; this checks consistency properties

ring r = 32003, (x,y), ds;
poly p = y2+y3+x2y2;
ideal i = y2+x10y, x22;
ideal j = std(i);

// reduce(p,j,2) is tail-reduce: p - reduce(p,j,2) should be in ideal
size(reduce(p - reduce(p, j, 2), j));
// reduce(p,j) is full reduce: p - reduce(p,j) should be in ideal
size(reduce(p - reduce(p, j), j));
// full reduce should give zero remainder when applied twice
reduce(reduce(p,j),j) == reduce(p,j);

LIB "tst.lib"; tst_status(1);$
