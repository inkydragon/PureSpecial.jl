# Implementations

- [DLMF: Software Index](https://dlmf.nist.gov/software/)

## Opensource

- C++ [scipy/xsf: Special function implementations](https://github.com/scipy/xsf)
- F77/C++/Python [Special functions (scipy.special) — SciPy dev Manual](https://scipy.github.io/devdocs/reference/special.html)
- Python [Mathematical functions — mpmath documentation](https://mpmath.org/doc/current/functions/index.html)
- C++ [Boost - Chapter 8. Special Functions](https://www.boost.org/doc/libs/release/libs/math/doc/html/special.html)
- [Functions - SageMath Documentation](https://doc.sagemath.org/html/en/reference/functions/index.html)
- C [Special Functions — GSL documentation](https://www.gnu.org/software/gsl/doc/html/specfunc.html)
- C/F77 [JuliaMath/openspecfun: A collection of special mathematical functions](https://github.com/JuliaMath/openspecfun):
    Faddeeva (C);  AMOS (f77)
- Rust [Axect/puruspe: PURe RUSt SPEcial function library](https://github.com/Axect/puruspe)
- Rust [cpmech/russell: Rust Scientific Libary. ODE and DAE (Runge-Kutta) solvers. Special functions (Bessel, Elliptic, Beta, Gamma, Erf). Linear algebra. Sparse solvers (MUMPS, UMFPACK). Probability distributions. Tensor calculus.](https://github.com/cpmech/russell)

### Netlib Code

- [f77, 1975] FUNPACK:
    W. J. Cody, 1975, The FUNPACK Package of Special Function Subroutines
- [f77, 1977/1992], [netlib/fn](https://www.netlib.org/fn/), [netlib/slatec/fnlib](https://www.netlib.org/slatec/fnlib):
    Fullerton, W., (LANL), 1977~1992

  - function index: [netlib.org/fn/slatec_index](https://www.netlib.org/fn/slatec_index)

- [f77, 1983~1989], [netlib/amos](https://www.netlib.org/amos/):
    D. E. Amos, 1986, Sandia National Laboratories.  
    Algorithm 644: A portable package for Bessel functions of a complex argument and nonnegative order

    ```fortran
    ! Airy functions Ai(z), Ai′(z), Bi(z), Bi′(z)
    cairy, cbiry; 
    ! Bessel functions Hv(1)(z), Hv(2)(z), Iv(z), Jv(z), Kv(z), Yv(z) 
    cbesh;  cbesi, cbesj, cbesk, cbesy
    ```

- [f77, 1988], [netlib/specfun](https://www.netlib.org/specfun/):
    SPECFUN 2.5, by W. J. Cody, Argonne National Laboratory.  
    Algorithm 715: SPECFUN–a portable FORTRAN package of special function routines and test drivers

    ```fortran
    gamma, loggamma; psi; dawson; ei, e1, eix; erf, erfc, erfcx;
    besi0, besei0, besi1, besei1; besj0, besj1; besy0, besy1; besk0, besek0, besk1, besek1;
    I_{N+ALPHA}(X), J_{N+ALPHA}(X), K_{N+ALPHA}(X), Y_{N+ALPHA}(X);
    ```
- [C, 1989], [netlib/cephes](https://www.netlib.org/cephes):
    Stephen L. Moshier, 1989, Cephes Mathematical Library

- [f77, 1993], [netlib/vfnlib](https://www.netlib.org/vfnlib/):
    Ronald F. Boisvert and Bonita V. Saunders, 1993,
    Portable Vectorized Software for Bessel Function Evaluation

    VFNLIB is a suite of Fortran subprograms for the evaluation of Bessel
    functions and modified Bessel functions of orders zero and one for a
    vector of real arguments.
    ```vi0, vi1; vj0, vj1; vk0, vk1; vy0, vy1```

- [f77, 1996] specfun.f,:
    Shanjie Zhang and Jianming Jin, 1996.  
    "Computation of Special Functions", 1996, John Wiley & Sons, Inc.

- [f77, 1998] [netlib/math (math77,mathc90)](https://www.netlib.org/math/):
    Fred T. Krogh, Charles L. Lawson, W. Van Snyder, 1992~1998, California Institute of Technology.

    MATH77 and mathc90, Release 6.0  
    Libraries of Mathematical Subprograms in Fortran 77 and C

## Commercial

- [Special Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/SpecialFunctions.html)
- [IMSL® Fortran Math Special Functions Library Special Functions User Guide](https://help.imsl.com/fortran/fnlspecfunc/current/index.htm)
- [S (Specfun) : NAG Library Manual, Latest](https://support.nag.com/numeric/nl/nagdoc_latest/flhtml/s/sconts.html)
- [Special Functions - MATLAB & Simulink](https://www.mathworks.com/help/matlab/special-functions.html)
- MATLAB [Multiprecision Computing Toolbox for MATLAB](https://www.advanpix.com/documentation/function-reference/#Special_Functions)
- [Special Functions - Intel® oneAPI MKL](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-0/special-functions.html)


## Research

- [Rational Approximations for the Modified Bessel Function of the First Kind – I1(x) for Computations with Double Precision](https://www.advanpix.com/2015/11/12/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i1-for-computations-with-double-precision/)
