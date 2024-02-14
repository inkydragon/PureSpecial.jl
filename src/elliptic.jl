# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Elliptic functions and integrals

## cephes
ellipj(u, m)
    Jacobian elliptic functions
    + cephes_ellpj

ellipk(m)
    Complete elliptic integral of the first kind.
    + ellipk -> ellpk -> cephes_ellpk

ellipkm1(p)
    Complete elliptic integral of the first kind around m = 1
    + cephes_ellpk

ellipkinc(phi, m)
    Incomplete elliptic integral of the first kind
    + cephes_ellik

ellipe(m)
    Complete elliptic integral of the second kind
    + cephes_ellpe

ellipeinc(phi, m)
    Incomplete elliptic integral of the second kind
    + cephes_ellie

    
## ellint_carlson
elliprc(x, y)
    Degenerate symmetric elliptic integral.
    + fellint_RC,cellint_RC

elliprd(x, y, z)
    Symmetric elliptic integral of the second kind.
    + fellint_RD,cellint_RD

elliprf(x, y, z)
    Completely-symmetric elliptic integral of the first kind.
    + fellint_RF,cellint_RF

elliprg(x, y, z)
    Completely-symmetric elliptic integral of the second kind.
    + fellint_RG,cellint_RG

elliprj(x, y, z, p)
    Symmetric elliptic integral of the third kind.
    + fellint_RJ,cellint_RJ
"""
