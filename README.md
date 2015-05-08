# poly-roots
[![Build Status](https://travis-ci.org/scijs/poly-roots.svg?branch=master)](https://travis-ci.org/scijs/poly-roots) [![npm version](https://badge.fury.io/js/poly-roots.svg)](http://badge.fury.io/js/poly-roots)  [![Dependency Status](https://david-dm.org/scijs/poly-roots.svg)](https://david-dm.org/scijs/poly-roots)

Find all [roots](http://en.wikipedia.org/wiki/Root_of_a_function) of a [polynomial](http://en.wikipedia.org/wiki/Polynomial) using the [Jenkins-Traub method](http://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm). In other words, it [factorizes the polynomial](http://en.wikipedia.org/wiki/Factorization_of_polynomials) over the complex numbers.

**N.B.: I fear I strayed too far toward translating cpoly while trying to understand the algorithm. It's similar enough to likely be covered under the original [ACM Software License Agreement](http://www.acm.org/publications/policies/softwarecrnotice). Sorry.**

## Introduction

This module factors a polynomial of the form 

![a0 * z^n + a1 * z^(n-1) + ... + a\_n-1 z + a\_n](docs/images/poly.png).

It uses the [Jenkins-Traub method](http://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm), and more specifically it's very nearly a line-by-line translation of the tried and true [cpoly.f90](http://jblevins.org/mirror/amiller/cpoly.f90). No really, it's almost a direct translation, taking some leeway in reworking goto statements into javascript. I started off with a pretty naive implementation of the original paper, [A three-stage variable shift iteration for polynomial zeros and its relation to generalized Ralyeigh iteration](http://octopus.library.cmu.edu/Collections/traub62/box00027/fld00056/bdl0004/doc0001/doc_27b56f4b1.pdf) by M. A. Jenkins and J. F. Traub, 1968, but there are some serious shortcuts and simplifications you can take if you stop and think about what you're doing. So I gave up cleaning up and refactoring my own version and reworked an existing implementation into JavaScript.

**The good**:

- It's reasonably fast
- It's numerically stable
- Memory usage is linear
- It benefits from the experimentation of the people who originally sat down and came up with a great implementation
- No dependencies

**The bad**:
- It's been translated by hand.
- The convergence criteria need a bit of work. I glossed over a couple subroutines that juggle some operations in order to prevent underflow errors, so I suspect the error estimates relative to machine epsilon aren't stricly accurate.
- It can maybe be translated better and more effectively via f2c + emscripten.
- The speed can be cut in half for polynomials with real coefficients by using the [rpoly.f90](http://jblevins.org/mirror/amiller/rpoly.f90) variant

You can go do some research about good root-finders, but for a quick rundown of what you have to work with if you want to stick with JavaScript, see a [quick benchmark](https://github.com/rreusser/scijs-polynomial-roots-comparison).

## Usage

#### `require('poly-roots')( real_coeffs [, imag_coeffs] )`

Computes the roots of a polynomial given the coefficients in descending order.

- `real_coeffs` the real coefficients of the polynomial arranged in order of decreasing degree
- `imag_coeffs` (optional) the imaginary coefficients of the polynomial. If not specified, assumed to be zero

**Returns**:  A pair of vectors representing the real and imaginary parts of the roots of the polynomial


## Example

```javascript
var roots = require('poly-roots');

// Roots of x^2 + 2x - 3:
var r1 = roots([1,2,-3]);

// Roots of z^3 - (4 + i)z^2 + (1 + i)z + (6 + 2i):
var r2 = roots([1,-4,1,6],[0,-1,1,2]);
```

## See also
For the companion roots version that determines roots by solution of an eigenvalue problem (via numeric.js), see [companion-roots](https://github.com/scijs/companion-roots).
For a blazing fast variant that might struggle in corner cases (like closely-spaced roots), see [durand-kerner](https://github.com/scijs/durand-kerner).

## Credits

Since this inherits a lot from `cpoly.f90` and `cpoly.f90` in turn is an update of the original code from CACM 419, I'm afraid that it may be subject to the [ACM Software License Agreement](http://www.acm.org/publications/policies/softwarecrnotice) which, in short, grants to you a royalty-free, nonexclusive right to execute, copy, modify and distribute both the binary and source code solely for academic, research and other similar noncommercial uses. :(
