# Fortran Impl Ref

## specfun

last commit: https://github.com/scipy/scipy/blob/23ed5fed8f9450446fc0aad0acb5002d4d7f84f7/scipy/special/specfun/specfun.f

```fortran
C       COMPUTATION OF SPECIAL FUNCTIONS
C
C          Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs.
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
C
C       Scipy changes:
C       - Compiled into a single source file and changed REAL To DBLE throughout.
C       - Changed according to ERRATA.
C       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
C       - Made functions return sf_error codes in ISFER variables instead
C         of printing warnings. The codes are
C         - SF_ERROR_OK        = 0: no error
C         - SF_ERROR_SINGULAR  = 1: singularity encountered
C         - SF_ERROR_UNDERFLOW = 2: floating point underflow
C         - SF_ERROR_OVERFLOW  = 3: floating point overflow
C         - SF_ERROR_SLOW      = 4: too many iterations required
C         - SF_ERROR_LOSS      = 5: loss of precision
C         - SF_ERROR_NO_RESULT = 6: no result obtained
C         - SF_ERROR_DOMAIN    = 7: out of domain
C         - SF_ERROR_ARG       = 8: invalid input parameter
C         - SF_ERROR_OTHER     = 9: unclassified error
C       - Improved initial guesses for roots in JYZO.
```
