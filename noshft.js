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
}

module.exports = function noshft (l1, nn, hr, hi, pr, pi) {
  var i, j, jj, nm1, n;
  var xni, t1, t2, tmp;

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

  console.log('P =',showPolyPython(nn, pr, pi));
  console.log('H =',showPolyPython(nn-1, hr, hi));

  for(jj=0; jj<l1; jj++) {

    // Compare the trailing coefficient of H to that of P:
    if( Math.sqrt(hr[nm1]*hr[nm1]+hi[nm1]*hi[nm1]) > 1e-15 * Math.sqrt(pr[nm1]*pr[nm1]+pi[nm1]*pi[nm1]) ) {

      // t = - p[n] / h[nm1]
      tmp = hr[nm1]*hr[nm1]+hi[nm1]*hi[nm1];
      tr = - ( pr[n]*hr[nm1]+pi[n]*hi[nm1] ) / tmp;
      ti = - ( pi[n]*hr[nm1]-pr[n]*hi[nm1] ) / tmp;

      for(i=0; i<nm1; i++) {
        j = nn - i - 2;
        t1 = hr[j-1];
        t2 = hi[j-1];
        hr[j] = tr * t1 - ti * t2 + pr[j];
        hi[j] = tr * t2 + ti * t1 + pi[j];
      }
      hr[0] = pr[0];
      hi[0] = pi[0];

      console.log('H =',showPolyPython(n, hr, hi));

    } else {

      // If the constant term is essentially zero, shift h coefficients
      for(i=0; i<nm1; i++) {
        j = nn - i;
        hr[j-1] = hr[j-2];
        hi[j-1] = hi[j-2];
      }
      hr[0] = 0;
      hi[0] = 0;
    }
  }

  return;

}
