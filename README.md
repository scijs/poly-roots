# poly-roots
=====================

Find all [roots](http://en.wikipedia.org/wiki/Root_of_a_function) of a [polynomial](http://en.wikipedia.org/wiki/Polynomial) using the [Jenkins-Traub method](http://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm). In other words, it [factorizes the polynomial](http://en.wikipedia.org/wiki/Factorization_of_polynomials) over the complex numbers.

## Introduction

This module factors a polynomial of the form a0 * z^n + a1 * z^(n-1)

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

# Credits
(c) 2015 Ricky Reusser. MIT License
