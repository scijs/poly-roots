'use strict';

/*function show(a,b) {
  var i,nums;
  nums = [];
  for(i=0; i<a.length; i++) {
    nums.push('('+a[i]+'+'+b[i]+'i)');
  }
  return '['+nums.join(', ')+']';
}

function showPoly(n,a,b) {
  var i, terms;
  if( b===undefined ) {
    terms = [];
    for(i=0; i<n; i++) {
      terms.push(''+a[i]+' * z^' + (n-i-1));
    }
    return terms.join(' + ');
  } else {
    terms = [];
    for(i=0; i<n; i++) {
      terms.push('( '+a[i]+' + '+b[i]+'*%i ) * z^' + (n-i-1));
    }
    return terms.join(' + ');
  }
}
function showPolyPython(n,a,b) {
  var i, terms;
  if( b===undefined ) {
    terms = [];
    for(i=0; i<n; i++) {
      terms.push(''+a[i]+' * z^' + (n-i-1));
    }
    return terms.join(' + ');
  } else {
    terms = [];
    for(i=0; i<n; i++) {
      terms.push(a[i]+' + '+b[i]+'j');
    }
    return '[' + terms.join(', ') + ']';
  }
}*/

function polyev ( nn, sr, si, pr, pi, qr, qi, pvri ) {
  //
  // Evaluates a polynomial P at s by the Horner recurrence, placing the
  // partial sums in Q and the computed value in pvri, which for JS is
  // obviously an array since we can't pass by ref.
  //
  var i,t,pvr,pvi;

  pvr = qr[0] = pr[0];
  pvi = qi[0] = pi[0];

  for(i=1; i<nn; i++) {
    t = pvr * sr - pvi * si + pr[i];
    pvi = pvr * si + pvi * sr + pi[i];
    pvr = t;
    qr[i] = pvr;
    qi[i] = pvi;
  }
  pvri[0] = pvr;
  pvri[1] = pvi;
}

function calct( nn, sr, si, hr, hi, qhr, qhi, pvri, tri ) {
  // Computes t = -P(s)/H(s)
  // Returns 

  var bool, tmp;
  var n = nn - 1;
  var hvri = new Float64Array(2);

  polyev( n, sr, si, hr, hi, qhr, qhi, hvri );

  var m1 = Math.sqrt(hvri[0]*hvri[0]+hvri[0]*hvri[0]);
  var m2 = Math.sqrt(hr[n-1]*hr[n-1]+hi[n-1]*hi[n-1]);
  bool = (m1 <= 1e-15 * m2);

  if( ! bool ) {
    tmp = hvri[0]*hvri[0] + hvri[1]*hvri[1];
    tri[0] = - ( pvri[0]*hvri[0] + pvri[1]*hvri[1] ) / tmp;
    tri[1] = - ( pvri[1]*hvri[0] - pvri[0]*hvri[1] ) / tmp;
    return bool;
  }

  tri[0] = 0;
  tri[1] = 0;

  return bool;
}

function nexth( nn, bool, qhr, qhi, qpr, qpi, hr, hi, tri ) {
  // Calculates the next shifted H polynomial
  // return true if H(s) is essentially zero
  var t1, t2, j;
  var n = nn - 1;
  
  if( ! bool ) {
    for(j=1; j<n; j++) {
      t1 = qhr[j-1];
      t2 = qhi[j-1];
      hr[j] = tri[0] * t1 - tri[1] * t2 + qpr[j];
      hi[j] = tri[0] * t2 + tri[1] * t1 + qpi[j];
    }
    hr[0] = qpr[0];
    hi[0] = qpi[0];
    return;
  }

  // If H(s) is zero, replace h with qh:
  for(j=1; j<n; j++) {
    hr[j] = qhr[j-1];
    hi[j] = qhi[j-1];
  }
  hr[0] = 0;
  hi[0] = 0;
}

function errev ( nn, qr, qi, ms, mp ) {
  // Bounds the error in evaluating the polynomial by Horner recurrence
  // qr, qi: the partial sums
  // ms: modulus of the point
  // mp: modulus of the polynomial value
  // are, mre: error bounds on complex addition and multiplication

  var e, i;
  var are = 1.1e-16;
  var mre = 3.11e-16;

  e = Math.sqrt(qr[0]*qr[0]+qi[0]*qi[0]) * mre / (are + mre);
  for(i=0;i<nn;i++) {
    e = e * ms + Math.sqrt(qr[i]*qr[i]+qi[i]*qi[i]);
  }
  return e * (are + mre) - mp * mre;
}

function cauchy (nn, pt, q ) {
  // Cauchy computs a lower bound on the moduli ofthe zeros of a polynomial.
  // pt is the modulus of the coefficients.
  //
  // The lower bound of the modulus of the zeros is given by the roots of the polynomial:
  //
  //   n         n-1
  //  z  + |a | z    + ... + |a   | z - |a |
  //         1                 n-1        n
  //
  var x, dx, df, i, n, xm, f;

  // The sign of the last coefficient is reversed:
  pt[nn-1] = -pt[nn-1];

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

function noshft (l1, nn, hr, hi, pr, pi) {
  var i, j, jj, nm1, n, tr, ti;
  var xni, t1, t2, tmp;

  //console.log('begin stage 1');

  n = nn - 1;
  nm1 = n - 1;

  // From Eqn. 9.3 for the 'scaled recurrence', calculate
  //
  //  _ (0)        1
  //  H    (z)  =  - P' (z)
  //               n
  //
  // In my copy of the paper, the 'P' appears to be missing a prime. Since the
  // leading coefficient is constrained to be 1, this makes H a monic polynomial.
  for(i=0; i<n; i++) {
    xni = nn - i - 1;
    hr[i] = xni * pr[i] / n;
    hi[i] = xni * pi[i] / n;
  }

  for(jj=0; jj<l1; jj++) {
    //console.log('No shift iteration #',jj+1,'H =',showPolyPython(n,hr,hi));
    //console.log('No shift iteration #',jj+1,'P =',showPolyPython(n,pr,pi));

    // Compare the trailing coefficient of H to that of P:
    //console.log( Math.sqrt(hr[nm1]*hr[nm1]+hi[nm1]*hi[nm1]));
    //console.log( Math.sqrt(pr[n]*pr[n]+pi[n]*pi[n]));
    if( (hr[nm1]*hr[nm1]+hi[nm1]*hi[nm1]) > 1e-30 * (pr[n]*pr[n]+pi[n]*pi[n]) ) {

      // t = - p[n] / h[nm1]
      tmp = - (hr[nm1]*hr[nm1]+hi[nm1]*hi[nm1]);
      tr = ( pr[n]*hr[nm1]+pi[n]*hi[nm1] ) / tmp;
      ti = ( pi[n]*hr[nm1]-pr[n]*hi[nm1] ) / tmp;

      // Calculate the new polynomial using the horner recurrence:
      for(j=nm1; j>0; j--) {
        t1 = hr[j-1];
        t2 = hi[j-1];
        hr[j] = tr * t1 - ti * t2 + pr[j];
        hi[j] = tr * t2 + ti * t1 + pi[j];
      }
      hr[0] = pr[0];
      hi[0] = pi[0];

    } else {
      //console.log('edge case');

      // If the constant term is essentially zero, shift h coefficients
      for(i=n-1; i>0; i--) {
        //console.log("shifting",i);
        hr[i] = hr[i-1];
        hi[i] = hi[i-1];
      }
      hr[0] = 0;
      hi[0] = 0;
      //console.log('shifted H =',showPolyPython(nm1,hr,hi));
    }
  }
}

function vrshft (l3, nn, zri, sri, hr, hi, pr, pi, qpr, qpi, qhr, qhi, shr, shi, pvri, tri) {

  var mp, ms, omp, relstp, r1, r2, tp, bool, i, j;
  
  //console.log('begin stage 3');

  var conv = false;
  var b = false;

  sri[0] = zri[0];
  sri[1] = zri[1];
  var eta = 1.1e-16;


  // Main loop for stage three
  for(i=0; i<l3; i++) {
    //console.log('variable shift iteration #',i+1, ', H =',showPolyPython(nn-1,hr,hi));

    // Evaluate P at s and test for convergence
    polyev( nn, sri[0], sri[1], pr, pi, qpr, qpi, pvri );
    mp = Math.sqrt( pvri[0]*pvri[0] + pvri[1]*pvri[1] );
    ms = Math.sqrt( sri[0]*sri[0] + sri[1]*sri[1] );
    var err = errev(nn, qpr, qpi, ms, mp);
    //console.log('compare mp=',mp,' to err=',err);
    if( mp < 20 * err ) {
      // Polynomial value is smaller in value than a bound on the error in evaluating P,
      // terminate the iteration
      conv = true;
      zri[0] = sri[0];
      zri[1] = sri[1];
      //console.log('converged to',zri[0],'+',zri[1]+'i');
      return conv;
    }

    if( i !== 0 ) {
      if( ! ( b || mp < omp || relstp >= 0.5 ) ) {
        // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed
        // shift steps into the cluster to force one zero to dominate.
        tp = relstp;
        b = true;
        if( relstp < eta ) {
          tp = eta;
        }
        r1 = Math.sqrt(tp);
        r2 = sri[0] * (1+r1) - sri[1] * r1;
        sri[1] = sri[0] * r1 + sri[1] * (1+r1);
        sri[0] = r2;
        polyev( nn, sri[0], sri[1], pr, pi, qpr, qpi, pvri );
        for(j=1; j<5; j++) {
          calct( nn, sri[0], sri[1], hr, hi, qhr, qhi, pvri, tri );
          nexth( nn, bool, qhr, qhi, qpr, qpi, hr, hi, tri );
        }
        omp = Infinity;
      }
      if( mp * 0.1 > omp ) {
        return conv;
      }
    }
    omp = mp;

    bool = calct( nn, sri[0], sri[1], hr, hi, qhr, qhi, pvri, tri );
    bool = nexth( nn, bool, qhr, qhi, qpr, qpi, hr, hi, tri );
    bool = calct( nn, sri[0], sri[1], hr, hi, qhr, qhi, pvri, tri );
    if( ! bool ) {
      relstp = Math.sqrt(tri[0]*tri[0]+tri[1]*tri[1]) / Math.sqrt(sri[0]*sri[0]+sri[1]*sri[1]);
      sri[0] += tri[0];
      sri[1] += tri[1];
    }
  }
  return conv;
}

function fxshft (l2, nn, zri, sri, hr, hi, pr, pi, qpr, qpi, qhr, qhi, shr, shi, pvri, tri) {
  var i, j, n, conv;
  var otr, oti, bool, svsr, svsi;

  n = nn - 1;

  //console.log('begin stage 2');

  //console.log('np.polyval(',showPolyPython(nn,pr,pi)+', '+sr+'+'+si+'j)');
  //console.log('np.polyval(',showPolyPython(n,hr,hi)+', '+sr+'+'+si+'j)');

  // Evaluate P at s:
  polyev( nn, sri[0], sri[1], pr, pi, qpr, qpi, pvri );
  var test = true;
  var pasd = false;

  // Calculate first t = -P(s)/H(s):
  bool = calct( nn, sri[0], sri[1], hr, hi, qhr, qhi, pvri, tri );
  zri[0] = sri[0] + tri[0];
  zri[1] = sri[1] + tri[1];

  // Main loop for one second stage step
  for(j=0; j<l2; j++) {
    //console.log('fixed shift iteration #',j+1,' H =',showPolyPython(n,hr,hi));
    otr = tri[0];
    oti = tri[1];

    // Compute next h polynomial and new t:
    nexth( nn, bool, qhr, qhi, qpr, qpi, hr, hi, tri );
    bool = calct( nn, sri[0], sri[1], hr, hi, qhr, qhi, pvri, tri );
    zri[0] = sri[0] + tri[0];
    zri[1] = sri[1] + tri[1];


    // Test for convergence unless stage 3 has failed once or this is the last H polynomial
    if( ! ( bool || !test || j === l2-1) ) {
      var tmpi = tri[0] - otr;
      var tmpr = tri[1] - oti;
      var m1 = Math.sqrt(tmpr*tmpr + tmpi*tmpi);
      var m2 = Math.sqrt(zri[0]*zri[0] + zri[1]*zri[1]);
      if( m1 < 0.5 * m2 ) {
        if( pasd ) {
          //console.log('weak convergence passed twice');
          // The weak convergence test has been passed twice, start the third stage
          // iteration, after saving the current H polynomial and shift
          for(i=0; i<n; i++) {
            shr[i] = hr[i];
            shi[i] = hi[i];
          }
          svsr = sri[0];
          svsi = sri[1];

          conv = vrshft( 10, nn, zri, sri, hr, hi, pr, pi, qpr, qpi, qhr, qhi, shr, shi, pvri, tri );
          if( conv ) {
            return conv;
          }

          // The iteration failed to converge. Turn off testing and restore H, S, PV, and T.
          //console.log('iteration failed to converge. Turn off testing & restore H, S, PV, T');
          test = false;
          for(i=0; i<n; i++) {
            hr[i] = shr[i];
            hi[i] = shi[i];
          }
          sri[0] = svsr;
          sri[1] = svsi;
          polyev( nn, sri[0], sri[1], pr, pi, qpr, qpi, pvri );
          bool = calct( nn, sri[0], sri[1], hr, hi, qhr, qhi, pvri, tri );
          continue;
        } 
        pasd = true;
      } else {
        pasd = false;
      }
    }
  }

  conv = vrshft( 10, nn, zri, sri, hr, hi, pr, pi, qpr, qpi, qhr, qhi, shr, shi, pvri, tri );

  return conv;
}

var cpoly = function cpoly ( opr, opi ) {
  var i, bound, xxx, conv, cnt1, cnt2, tmp, idnn2;

  if( opi === undefined ) {
    opi = new Float64Array( opr.length );
  }

  if( opr.length !== opi.length ) {
    throw new TypeError('cpoly: error: real/complex polynomial coefficient input length mismatch');
  }

  var degree = opr.length - 1;

  // Initialization of constants
  var xx = 0.7071067811865475,
      yy = -xx,
      cosr = -0.06975647374412534,
      sinr = 0.9975640502598242,
      nn = degree + 1;

  // Output:
  var zeror = new Float64Array(degree),
      zeroi = new Float64Array(degree);

  // Allocate arrays
  var pr = new Float64Array(nn),
      pi = new Float64Array(nn),
      hr = new Float64Array(degree),
      hi = new Float64Array(degree),
      qpr = new Float64Array(nn),
      qpi = new Float64Array(nn),
      qhr = new Float64Array(degree),
      qhi = new Float64Array(degree),
      shr = new Float64Array(nn),
      shi = new Float64Array(nn),
      zri = new Float64Array(2), // for pasing around zeros since js doesn't do by reference;
      pvri = new Float64Array(2),  // for passing around the polynomial value, reason = ditto
      tri = new Float64Array(2),   // for passing around T = -P(S)/H(S)
      sri = new Float64Array(2);   // for passing around current s
  
  //console.log('degree =',degree);
  //console.log('nn =',nn);
  // Remove the zeros at the origin if any
  while( opr[nn] === 0 && opi[nn] === 0 ) {
    idnn2 = degree - nn + 1;
    zeror[idnn2] = 0;
    zeroi[idnn2] = 0;
    nn--;
  }

  // Make a copy of the coefficients
  for(i=0; i<nn; i++) {
    pr[i] = opr[i];
    pi[i] = opi[i];
    shr[i] = Math.sqrt(pr[i]*pr[i]+pi[i]*pi[i]);
  }

  // Skip scaling the polynomial for unusually large
  // or small coefficients. Caveat emptor.

  // Start the algorithm for one zero
  while( nn > 2 ) {

    // Calculate bound, a lower bound on the modulus of the zeros:
    for(i=0; i<nn; i++) {
      shr[i] = Math.sqrt(pr[i]*pr[i] + pi[i]*pi[i]);
    }
    bound = cauchy( nn, shr, shi);

    // Outer loop to control two major passes with different sequences of shifts
    for(cnt1=0; cnt1<2; cnt1++) {
      //console.log("BEGIN OUTER LOOP");

      noshft(5, nn, hr, hi, pr, pi);

      // Inner loop to select a shift
      for(cnt2=0; cnt2<9; cnt2++) {
        //console.log("BEGIN INNER LOOP");

        // rotate shift angle xx and yy:
        xxx = cosr * xx - sinr * yy;
        yy = sinr * xx + cosr * yy;
        xx = xxx;

        // Calculate the new shift:
        sri[0] = bound * xx;
        sri[1] = bound * yy;

        conv = fxshft( 10*cnt2, nn, zri, sri, hr, hi, pr, pi, qpr, qpi, qhr, qhi, shr, shi, pvri, tri );

        if( conv ) {
          //console.log('MAIN LOOP CONVERGENCE. STORE ZERO');
          idnn2 = degree - nn + 1;
          zeror[idnn2] = zri[0];
          zeroi[idnn2] = zri[1];
          nn--;
          //console.log('nn-- = ',nn);

          for(i=0; i<nn; i++) {
            pr[i] = qpr[i];
            pi[i] = qpi[i];
          }
          cnt1 = 3; // exit from the cnt2 loop also
          break;
        }
      }
    }
  }

  if( nn <= 2 ) {
    //console.log('END LOOP CONVERGENCE. STORE ZERO');
    // Calculate the final zero and return ( - p[1] / p[0] ):
    tmp = pr[0]*pr[0] + pi[0]*pi[0];
    zeror[degree-1] = - ( pr[1]*pr[0] + pi[1]*pi[0] ) / tmp;
    zeroi[degree-1] = - ( pi[1]*pr[0] - pr[1]*pi[0] ) / tmp;
  }

  return [zeror, zeroi];
};

module.exports = cpoly;
