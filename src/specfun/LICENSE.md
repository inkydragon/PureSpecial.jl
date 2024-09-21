The code implementation in this folder is translated from the C code in [`scipy/special/specfun.c`](https://github.com/scipy/scipy/blob/2b8417343ebe919599e6391ac9787001505f9d34/scipy/special/specfun.c) to julia.

The license for the SciPy project is [BSD 3-Clause "New" or "Revised" License](https://github.com/scipy/scipy/blob/main/LICENSE.txt).

See below for the license.
```
/*
 *
 * This file accompanied with the header file specfun.h is a partial
 * C translation of the Fortran code by Zhang and Jin following
 * original description:
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *       COMPUTATION OF SPECIAL FUNCTIONS
 *
 *          Shanjie Zhang and Jianming Jin
 *
 *       Copyrighted but permission granted to use code in programs.
 *       Buy their book:
 *
 *          Shanjie Zhang, Jianming Jin,
 *          Computation of Special Functions,
 *          Wiley, 1996,
 *          ISBN: 0-471-11963-6,
 *          LC: QA351.C45.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *       Scipy changes:
 *       - Compiled into a single source file and changed REAL To DBLE throughout.
 *       - Changed according to ERRATA.
 *       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
 *       - Made functions return sf_error codes in ISFER variables instead
 *         of printing warnings. The codes are
 *         - SF_ERROR_OK        = 0: no error
 *         - SF_ERROR_SINGULAR  = 1: singularity encountered
 *         - SF_ERROR_UNDERFLOW = 2: floating point underflow
 *         - SF_ERROR_OVERFLOW  = 3: floating point overflow
 *         - SF_ERROR_SLOW      = 4: too many iterations required
 *         - SF_ERROR_LOSS      = 5: loss of precision
 *         - SF_ERROR_NO_RESULT = 6: no result obtained
 *         - SF_ERROR_DOMAIN    = 7: out of domain
 *         - SF_ERROR_ARG       = 8: invalid input parameter
 *         - SF_ERROR_OTHER     = 9: unclassified error
 *       - Improved initial guesses for roots in JYZO.
 *
 *
 */

/*
 * Copyright (C) 2024 SciPy developers
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Names of the SciPy Developers may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

```
