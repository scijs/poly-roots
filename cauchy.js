module.exports = function cauchy (nn, pt, q ) {
  // Cauchy computs a lower bound on the moduli ofthe zeros of a polynomial.
  // pt is the modulus of the coefficients.
  //
  // The lower bound of the modulus of the zeros is given by the roots of the polynomial:
  //
  //   n         n-1
  //  z  + |a | z    + ... + |a   | z - |a |
  //         1                 n-1        n
  //
  var x, xmf, dx, df, i;

  // The sign of the last coefficient is reversed:
  pt[nn-1] = -pt[nn-1]

  // Compute upper estimate of bound:
  n = nn - 1;
  x = Math.exp((Math.log(-pt[nn-1]) - Math.log(pt[0])) / n);

  if( pt[n-1] !== 0 ) {
    xm = - pt[nn-1] / pt[n-1];
    if( xm < x ) {
      x = xm;
    }
  }

  // Chop the interval (0,x) until f <= 0
  while( true ) {
    xm = x * 0.1;
    f = pt[0];
    for(i=1; i<nn; i++) {
      f = f * xm + pt[i];
    }
    if( f > 0 ) {
      x = xm;
    } else {
      break;
    }
  }
  dx = x;

  // Do newton iteration until x converges to two decimal places
  while( Math.abs(dx/x) > 0.005 ) {
    q[0] = pt[0];
    for(i=1; i<nn; i++) {
      q[i] = q[i-1] * x + pt[i];
    }
    f = q[nn-1];
    df = q[0];
    for(i=1; i<n; i++) {
      df = df * x + q[i];
    }
    dx = f / df;
    x -= dx;
  }

  return x;
}
