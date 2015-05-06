"use strict"

function show(a,b) {
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
}


// Differentiates in place. Understood that the last entry is
// cut off (but still exists in memory) and that the order of
// the polynomial is reduced by 1.
function differentiate( n, Pr, Pi, Ppr, Ppi ) {
  var p,i;
  for(i=n-2; i>=0; i--) {
    p = n-i-1;
    Ppr[i] = Pr[i] * p;
    Ppi[i] = Pi[i] * p;
  }
}


// The lower bound of the modulus of the zeros is given by the roots of the polynomial
//
//   n         n-1
//  z  + |a | z    + ... + |a   | z - |a |
//         1                 n-1        n
//
function lowerRootBound(n, Pr, Pi, w) {
  var tol = 1e-8;
  var maxiter = 15;
  var i=0, dx=Infinity, iter=0, f, fp, betaAccum, xn, dx;
  // An initial guess (will always work, I'm pretty sure, since the poly is
  // strictly increasing in the right half-plane):
  var x = 2.0;

  // Set up the coefficients so that the first coefficient is, by definition, one:
  var w0 = Math.sqrt(Pr[0]*Pr[0] + Pi[0]*Pi[0]);
  w[0] = 1;
  for(var i=1; i<n; i++) {
    w[i] = Math.sqrt(Pr[i]*Pr[i] + Pi[i]*Pi[i]) / w0;
  }

  while(Math.abs(dx) > tol) {
    // Calculate the function and its derivative in one go so that we don't
    // have to calculate powers of x twice for no reason:
    f = - w[n-1]; // <-- trailing term is inverted
    fp =  w[n-2];
    xn = x;
    for(i=n-2; i>=1; i--) {
      f += w[i] * xn;
      fp += (n-i) * w[i-1] * xn;
      xn *= x;
    }
    f += w[i] * xn;
    dx = -f/fp;
    x += dx;

    if( ++iter >= maxiter ) {
      return 0;
    }
  }
  return x;
}


// Compute the no-shift iteration in-place wrt H:
//
//                                   (l)
//    (l+1)        1  /  (l)        H   (0)      \
//   H     (z)  =  - (  H   (z)  -  ------- P(z)  )
//                 z  \              P(0)        /
//
// Last coefficient vanishes, so we're left with a polynomial
// of order n-2, the same as H.
//
// P is the polynomial order n+1. H is the (pre-initialized) polynomial order n. Jenkins and
// Traub state that convergence criteria exist, but in general five iterations is sufficient
// to make the small roots stand out a bit more.
function noShiftIteration( n, Pr, Pi, Hr, Hi ) {
  var Pr0, Pi0, Hr0, Hi0, tmp1, tmp2, HPr, HPi, i, H0m;

  // Pull P(0):
  Pr0 = Pr[n-1];
  Pi0 = Pi[n-1];

  // Calculate coefficient H^(l)(0) / P(0):
  tmp1 = Pr0*Pr0 + Pi0*Pi0;
  HPr = (Hr[n-2] * Pr0 + Hi[n-2] * Pi0) / tmp1;
  HPi = (Hi[n-2] * Pr0 - Hr[n-2] * Pi0) / tmp1;

  // Calculate the new polynomial H^(l+1) in place:
  for(i=n-2; i>=1; i--) {
    Hr[i] = Hr[i-1] - HPr * Pr[i] + HPi * Pi[i];
    Hi[i] = Hi[i-1] - HPr * Pi[i] - HPi * Pr[i];
  }
  Hr[0] = - (HPr * Pr[0] - HPi * Pi[0]);
  Hi[0] = - (HPr * Pi[0] + HPi * Pr[0]);
  
  // For production: Normalize the answer so that it doesn't diverge
  H0m = Math.sqrt(Hr[0]*Hr[0] + Hi[0]*Hi[0]);
  for(i=0; i<n-1; i++) {
    Hr[i] /= H0m;
    Hi[i] /= H0m;
  }

  // For debugging: Normalize the polynomial to H0=1 so that it's easy to see convergence
  Hr0 = Hr[0];
  Hi0 = Hi[0];
  H0m = Hi0*Hi0 + Hr0*Hr0;
  for(i=0; i<n; i++) {
    tmp1 = Hr[i];
    Hr[i] = (Hr0*tmp1 + Hi0*Hi[i]) / H0m;
    Hi[i] = (Hr0*Hi[i] - Hi0*tmp1) / H0m;
  }
}


// Compute a fixed-shift iteration in-place wrt H:
//
//                                        (l)
//    (l+1)          1     /  (l)        H   (s)      \
//   H     (z)  =  -----  (  H   (z)  -  ------- P(z)  )
//                 z - s   \              P(s)        /
//
function fixedShiftIteration( n, Pr, Pi, Hr, Hi, sr, si ) {
  var Prs, Pis, Hrs, His, tmp1, tmp2, HPr, HPi, i, zr, zi, ar, ai, H0m, Hr0, Hi0;

  // Calculate H^(l)(s) and P(s) in tandem:
  Prs = Pr[n-1];
  Pis = Pi[n-1];
  Hrs = Hr[n-2];
  His = Hi[n-2];
  zr = sr;
  zi = si;
  for(i=n-2; i>0; i--) {
    // P += a_i * z^(n-i):
    Prs += Pr[i]*zr - Pi[i]*zi;
    Pis += Pr[i]*zi + Pi[i]*zr;

    // H += a_(i-1) * z^(n-i):
    Hrs += Hr[i-1]*zr - Hi[i-1]*zi;
    His += Hr[i-1]*zi + Hi[i-1]*zr;

    // Multiply by z:
    tmp1 = zr;
    zr = zr*sr - zi*si;
    zi = zi*sr + tmp1*si;
  }

  // One leftover term for P since it's a degree larger:
  Prs += Pr[0]*zr - Pi[0]*zi;
  Pis += Pr[0]*zi + Pi[0]*zr;

  // Calculate coefficient H^(l)(s) / P(s):
  tmp1 = Prs*Prs + Pis*Pis;
  HPr = (Hrs * Prs + His * Pis) / tmp1;
  HPi = (His * Prs - Hrs * Pis) / tmp1;

  // Calculate the new polynomial H^(l+1) * (z-s) in place:
  for(i=n-2; i>=1; i--) {
    Hr[i] = Hr[i-1] - HPr * Pr[i] + HPi * Pi[i];
    Hi[i] = Hi[i-1] - HPr * Pi[i] - HPi * Pr[i];
  }
  Hr[0] = - (HPr * Pr[0] - HPi * Pi[0]);
  Hi[0] = - (HPr * Pi[0] + HPi * Pr[0]);

  // Deflate by the root z-s, taking into account the fact that, by definition,
  // there's no remainder. Note that there's no change in the leading coefficient
  // whether or not it's equal to 1 because it just loses a z.
  for(i=1; i<n-1; i++) {
    Hr[i] += sr*Hr[i-1] - si*Hi[i-1];
    Hi[i] += sr*Hi[i-1] + si*Hr[i-1];
  }

  // For production: Normalize the answer so that it doesn't diverge
  H0m = Math.sqrt(Hr[0]*Hr[0] + Hi[0]*Hi[0]);
  for(i=0; i<n-1; i++) {
    Hr[i] /= H0m;
    Hi[i] /= H0m;
  }

  // For debugging: Normalize the polynomial to H0=1 so that it's easy to see convergence
  Hr0 = Hr[0];
  Hi0 = Hi[0];
  H0m = Hi0*Hi0 + Hr0*Hr0;
  for(i=0; i<n; i++) {
    tmp1 = Hr[i];
    Hr[i] = (Hr0*tmp1 + Hi0*Hi[i]) / H0m;
    Hi[i] = (Hr0*Hi[i] - Hi0*tmp1) / H0m;
  }
}

function evalPH( n, Pr, Pi, Hr, Hi, HP, sr, si ) {
  var zr, zi, i;
  // Calculate H^(l)(s) and P(s) in tandem:
  Prs = Pr[n-1];
  Pis = Pi[n-1];
  Hrs = Hr[n-2];
  His = Hi[n-2];
  zr = sr;
  zi = si;
  for(i=n-2; i>0; i--) {
    // P += a_i * z^(n-i):
    Prs += Pr[i]*zr - Pi[i]*zi;
    Pis += Pr[i]*zi + Pi[i]*zr;

    // H += a_(i-1) * z^(n-i):
    Hrs += Hr[i-1]*zr - Hi[i-1]*zi;
    His += Hr[i-1]*zi + Hi[i-1]*zr;

    // Multiply by z:
    tmp1 = zr;
    zr = zr*sr - zi*si;
    zi = zi*sr + tmp1*si;
  }
}

// Takes polynomial coefficients of format:
//
//       n         n+1
// a  * z   +  a  z    + ... + a   z  +  a
//  0           1               n-1       n
//
var polynomialRoots = function polynomialRoots( Pr, Pi ) {
  var i,j;

  if( Pi === undefined ) {
    Pi = new Float64Array(Pr.length);
  }

  if( Pr.length !== Pi.length ) {
    throw new TypeError("polynomialRoots: error: real and complex coefficient array lengths must be equal");
  }

  var n = Pr.length;
  
  // Iteration arrays and work arrays:
  var Hr = new Float64Array(n-1);
  var Hi = new Float64Array(n-1);
  var W = new Float64Array(n);

  // Initialize the stage one iteration with H(z) := P'(z)
  differentiate( n, Pr, Pi, Hr, Hi );

  var M = 3;
  var L = 5;
  var K = 30;
  // Five iterations of stage 1:
  for(i=0; i<M; i++) {
    noShiftIteration( n, Pr, Pi, Hr, Hi );
    console.log('H^(' + (i+1) + ') =',showPolyPython(n-1,Hr,Hi));
  }

  // Switch from stage 1 to stage 2. That means calculating a shift, which means taking the
  // lower bound on the zero and rotating it randomly about a circle in the complex plane:
  var lowerBound = lowerRootBound( n, Pr, Pi, W );
  var th = Math.random() * 2 * Math.PI;
  var sr = lowerBound * Math.cos(th);
  var si = lowerBound * Math.sin(th);

  console.log('initial shift =',sr,si);

  for(i=M; i<L; i++) {
    fixedShiftIteration( n, Pr, Pi, Hr, Hi, sr, si );
    console.log('H^(' + (i+1) + ') =',showPolyPython(n-1,Hr,Hi));
  }

  // Switch from stage 2 to stage 3. That means calculating a shift, which means taking the

  var zr, zi, Prs, Pis, Hrs, His, tmp1;
  for(i=L; i<K; i++) {

    // Calculate H^(l)(s) and P(s) in tandem:
    Prs = Pr[n-1];
    Pis = Pi[n-1];
    Hrs = Hr[n-2];
    His = Hi[n-2];
    zr = sr;
    zi = si;
    for(j=n-2; j>0; j--) {
      // P += a_i * z^(n-j):
      Prs += Pr[j]*zr - Pi[j]*zi;
      Pis += Pr[j]*zi + Pi[j]*zr;

      // H += a_(j-1) * z^(n-j):
      Hrs += Hr[j-1]*zr - Hi[j-1]*zi;
      His += Hr[j-1]*zi + Hi[j-1]*zr;

      // Multiply by z:
      tmp1 = zr;
      zr = zr*sr - zi*si;
      zi = zi*sr + tmp1*si;
    }

    // Update the shift for s:
    tmp1 = Hrs*Hrs + His*His;
    sr -= (Hrs*Prs + His*Pis) / tmp1;
    si -= (Hrs*Pis - His*Prs) / tmp1;

    console.log('updated shift =',sr,si);

    fixedShiftIteration( n, Pr, Pi, Hr, Hi, sr, si );
    console.log('H^(' + (i+1) + ') =',showPolyPython(n-1,Hr,Hi));
  }

};

module.exports = polynomialRoots;
