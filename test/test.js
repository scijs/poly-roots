var test = require('tape'),
    cpoly = require('../index.js');

// Search epected real/imag for observed real/imag pairs:
function assertContainsCloseTo( t, er, ei, or, oi, tol ) {
  var hits = [];
  var i,j,idx;
  for(i=0; i<er.length; i++) {
    idx = null;
    for(j=0; j<oi.length; j++) {
      if( i in hits ) continue;
      var err = Math.sqrt(Math.pow(er[i]-or[j],2) + Math.pow(ei[i]-oi[j],2));
      if( err < tol ) {
        hits.push(j);
        idx = j;
        break;
      }
    }
    if( idx === null ) {
      t.assert(false,'did not contain root ' + er[i] + ' + ' + ei[i] + 'i');
    }
  }
  t.assert(hits.length = er.length);
}

test("factor z - 6", function(t) {
  var roots = cpoly([1,-6]);
  assertContainsCloseTo( t, [6],[0], roots[0],roots[1], 1e-13 );
  t.end();
});

test("factor 2z - 12", function(t) {
  var roots = cpoly([2,-12]);
  assertContainsCloseTo( t, [6],[0], roots[0],roots[1], 1e-13 );
  t.end();
});

test("factor z - 6i", function(t) {
  var roots = cpoly([1,0],[0,-6]);
  assertContainsCloseTo( t, [0],[6], roots[0],roots[1], 1e-13 );
  t.end();
});

test("factor z^2 + 2z - 3", function(t) {
  var roots = cpoly([1,2,-3]);
  assertContainsCloseTo( t, [-3,1],[0,0], roots[0],roots[1], 1e-13 );
  t.end();
});

test("factor z^2 + 1", function(t) {
  var roots = cpoly([1,0,1]);
  assertContainsCloseTo( t, [0,0],[1,-1], roots[0],roots[1], 1e-13 );
  t.end();
});

test("factor z^3 - 4z^2 + z + 6", function(t) {
  // Calculate roots of
  //                                3      2
  //   (z - 3) (z - 2) (z + 1)  =  z  - 4 z  + z + 6
  //
  var a_r = new Float64Array([1,-4, 1, 6]);

  var roots = cpoly( a_r );

  assertContainsCloseTo( t, [3,2,-1],[0,0,0], roots[0],roots[1], 1e-13 );

  t.end();
});

test("factor (z-1) * (z+1) * (z + 1 + 1e-4i) * (z + 1 - 1e-4i)",function(t) {
  // Calculate roots of:
  //
  // (z - 1) * (z + 1 + 1e-4i) * (z + 1 - 1e-4i) * (z + 1)  =
  //
  //     4      3                        2
  //    x  + 2 x  + 9.99999993922529e-9 x  - 2 x - 1.00000001
  //

  var a_r = new Float64Array([1, 2, 9.99999993922529e-9, -2, -1.00000001 ]);

  var roots = cpoly( a_r );

  assertContainsCloseTo( t, [1, -1, -1, -1],[0, 1e-4, -1e-4, 0], roots[0],roots[1], 1e-7 );

  t.end();
});

test("factor z^3 - 4z^2 + z + 6", function(t) {
  //
  // Calculate roots of
  //                                3      2
  //   (z - 3) (z - 2) (z + 1)  =  z  - 4 z  + z + 6
  //
  var a_r = new Float64Array([1,-4, 1, 6]);

  var roots = cpoly( a_r );

  assertContainsCloseTo( t, [3,2,-1],[0,0,0], roots[0],roots[1], 1e-13 );

  t.end();
});

test("factor z^3 - (4 + i)z^2 + (1 + i)z + (6 + 2i)", function(t) {
  //
  // Calculate roots of
  //                                    3            2
  //   (z - 3 - i) (z - 2) (z + 1)  =  z  - (4 + i) z  + (1 + i) z + (6 + 2i)
  //
  // roots:
  //  3 + i
  //  2
  //  -1
  //
  var a_r = new Float64Array([1,-4,  1,  6]);
  var a_i = new Float64Array([0,-1,  1,  2]);

  var roots = cpoly( a_r, a_i );

  assertContainsCloseTo( t, [3,2,-1],[1,0,0], roots[0], roots[1], 1e-13 );

  t.end();
});

test("factor (z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)*(z-7)*(z-8)*(z-9)*(z-10)",function(t) {
  var a_r = new Float64Array([1,-55,1320,-18150,157773,-902055,3416930,-8409500,12753576,-10628640,3628800]);

  var roots = cpoly( a_r );
  
  var rr = [1,2,3,4,5,6,7,8,9,10];
  var ri = [0,0,0,0,0,0,0,0,0,0];

  assertContainsCloseTo( t, rr, ri, roots[0], roots[1], 1e-9 );

  t.end();
});


test("factor z^20 + i * z^10 + 1",function(t) {
  var a_r = new Float64Array([1,-55,1320,-18150,157773,-902055,3416930,-8409500,12753576,-10628640,3628800]);

  var roots = cpoly(new Float64Array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]),
                    new Float64Array([0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]) );
  
  t.end();
});
