var ndarray = require('ndarray'),
    test = require('tape'),
    cpoly = require('../index.js');

test("derivative", function(t) {

  //t.equals(deriv([1,1,1,1]).join(), [1,2,3].join())
  //
  // Calculate roots of
  //                                3      2
  //   (x - 3) (x - 2) (x + 1)  =  x  - 4 x  + x + 6
  //
  var a_r = new Float64Array([1,-4,  1,  6]);
  //var a_i = new Float64Array([0, 1, -1, -2]);
  var a_i = new Float64Array([0, 0, 0, 0]);

  cpoly( a_r, a_i, 3 );

  t.end()
});
