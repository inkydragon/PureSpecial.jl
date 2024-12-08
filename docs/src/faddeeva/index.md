# Faddeeva

> http://ab-initio.mit.edu/faddeeva/

Steven G. Johnson has written free/open-source C++ code
(with wrappers for C, Matlab, GNU Octave, Python, R, Scilab, and Julia)
to compute the various error functions of arbitrary complex arguments.

```cpp
// SPDX-License-Identifier: MIT

// compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
extern double complex Faddeeva_w(double complex z,double relerr);
extern double Faddeeva_w_im(double x); // special-case code for Im[w(x)] of real x

// Various functions that we can compute with the help of w(z)

// compute erfcx(z) = exp(z^2) erfc(z)
extern double complex Faddeeva_erfcx(double complex z, double relerr);
extern double Faddeeva_erfcx_re(double x);

// compute erf(z), the error function of complex arguments
extern double complex Faddeeva_erf(double complex z, double relerr);
extern double Faddeeva_erf_re(double x);

// compute erfi(z) = -i erf(iz), the imaginary error function
extern double complex Faddeeva_erfi(double complex z, double relerr);
extern double Faddeeva_erfi_re(double x);

// compute erfc(z) = 1 - erf(z), the complementary error function
extern double complex Faddeeva_erfc(double complex z, double relerr);
extern double Faddeeva_erfc_re(double x);

// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
extern double complex Faddeeva_Dawson(double complex z, double relerr);
extern double Faddeeva_Dawson_re(double x);
```
