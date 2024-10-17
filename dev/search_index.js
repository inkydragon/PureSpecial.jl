var documenterSearchIndex = {"docs":
[{"location":"special/#Special-functions","page":"Special functions","title":"Special functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"ref: scipy.special ## SciPy Manual   based commit: 2b84173","category":"page"},{"location":"special/","page":"Special functions","title":"Special functions","text":"cephes.h: 100\nspecfun_wrappers.h: 49\namos_wrappers.h: 14\nfaddeeva.h++: 9\nboost_special_functions.h++: 7","category":"page"},{"location":"special/#Airy-functions","page":"Special functions","title":"Airy functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] airy(z): cephesairy; amosairy,amos_biry\n[ ] airye(z): amosairy,amosbiry\n[ ] ai_zeros(nt): specfun.airyzo -> specfunairyzo\n[ ] bi_zeros(nt): specfun_airyzo\n[ ] itairy(x): specfun_itairy","category":"page"},{"location":"special/#Elliptic-functions","page":"Special functions","title":"Elliptic functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] ellipj(u, m): cephes_ellpj\n[ ] ellipk(m): cephes_ellpk\n[ ] ellipkm1(p): cephes_ellpk\n[ ] ellipkinc(phi, m): cephes_ellik\n[ ] ellipe(m): cephes_ellpe\n[ ] ellipeinc(phi, m): cephes_ellie\n[ ] elliprc(x, y): fellintRC,cellintRC\n[ ] elliprd(x, y, z): fellintRD,cellintRD\n[ ] elliprf(x, y, z): fellintRF,cellintRF\n[ ] elliprg(x, y, z): fellintRG,cellintRG\n[ ] elliprj(x, y, z, p): fellintRJ,cellintRJ","category":"page"},{"location":"special/#Bessel-functions","page":"Special functions","title":"Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] lmbda(v, x): specfun.lamv -> specfunlamv; specfun.lamn -> specfunlamn","category":"page"},{"location":"special/#Zeros-of-Bessel-functions","page":"Special functions","title":"Zeros of Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] jnjnp_zeros(nt): specfun.jdzo -> specfunjdzo\n[ ] jnyn_zeros(n, nt): specfun.jyzo -> specfunjyzo\njn_zeros(n, nt)\njnp_zeros(n, nt)\nyn_zeros(n, nt)\nynp_zeros(n, nt)\n[ ] specfun.cyzo -> specfuncyzo:\ny0_zeros(nt)\ny1_zeros(nt)\ny1p_zeros(nt)","category":"page"},{"location":"special/#Faster-Common-Bessel-functions","page":"Special functions","title":"Faster Common Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"Use cephes","category":"page"},{"location":"special/#Integrals-of-Bessel-functions","page":"Special functions","title":"Integrals of Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] itj0y0(x): specfun_itjya\n[ ] it2j0y0(x): specfun_ittjya\n[ ] iti0k0(x): specfun_itika\n[ ] it2i0k0(x): specfun_ittika\n[ ] besselpoly(a, lmb, nu): cephes_besselpoly","category":"page"},{"location":"special/#Derivatives-of-Bessel-functions","page":"Special functions","title":"Derivatives of Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"No 3rd deps.","category":"page"},{"location":"special/#Spherical-Bessel-functions","page":"Special functions","title":"Spherical Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"AMOS or pure cython.","category":"page"},{"location":"special/#Riccati-Bessel-functions","page":"Special functions","title":"Riccati-Bessel functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] riccati_jn: specfun.rctj -> specfunrctj\n[ ] riccati_yn: specfun.rcty -> specfunrcty","category":"page"},{"location":"special/#Struve-functions","page":"Special functions","title":"Struve functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"Use cephes, specfun","category":"page"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] struve: cephesstruveh\n[ ] modstruve: cephesstruvel\n[ ] itstruve0: specfun_itsh0\n[ ] it2struve0: specfun_itth0\n[ ] itmodstruve0: specfun_itsl0","category":"page"},{"location":"special/#Statistical-functions","page":"Special functions","title":"Statistical functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"Skip!","category":"page"},{"location":"special/#Information-Theory-functions","page":"Special functions","title":"Information Theory functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"Skip!","category":"page"},{"location":"special/#Gamma-functions","page":"Special functions","title":"Gamma functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"Use cephes","category":"page"},{"location":"special/#Error-function","page":"Special functions","title":"Error function","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"faddeeva,cephes,boost_special","category":"page"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] erf(z): faddeevaerf -> Faddeeva::erf; cepheserf\n[ ] erfc(x): faddeevaerfccomplex -> Faddeeva::erfc; cephes_erfc\n[ ] erfcx(x): faddeevaerfcx,faddeevaerfcx_complex -> Faddeeva::erfcx\n[ ] erfi(z): faddeevaerfi,faddeevaerfi_complex -> Faddeeva::erfi\n[ ] erfinv(y): \n[ ] erfcinv(y): \n[ ] wofz(z): \n[ ] dawsn(x): faddeevadawsn,faddeevadawsn_complex -> Faddeeva::Dawson\n[ ] fresnel(z): cephesfresnl; cfresnlwrap -> specfuncfs,specfuncfc\n[ ] fresnel_zeros(nt): \n[ ] modfresnelp(x): modifiedfresnelpluswrap -> specfunffk\n[ ] modfresnelm(x): modifiedfresnelminuswrap -> specfunffk\n[ ] voigt_profile(x, sigma, gamma): faddeevavoigtprofile -> Faddeeva::w\n[ ] erf_zeros: specfun.cerzo -> specfuncerzo\n[ ] fresnel_zeros: specfun.fcszo -> specfunfcszo\nfresnelc_zeros\nfresnels_zeros","category":"page"},{"location":"special/#Legendre-functions","page":"Special functions","title":"Legendre functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] sph_harm(m, n, theta, phi): cephes_poch\n[ ] lpmv(m, v, x): specfun_lpmv\n[ ] clpmn(m, n, z): specfun.clpmn -> specfunclpmn\n[ ] lpn(n, z): specfun.lpn -> specfunlpn; specfun.clpn -> specfunclpn\n[ ] lqn(n, z): specfun.lqnb -> specfunlqnb ; specfun.clqn -> specfunclqn\n[ ] lpmn(m, n, z): specfun.lpmn -> specfunlpmn\n[ ] lqmn(m, n, z): specfun.lqmn -> specfunlqmn; specfun.clqmn -> specfunclqmn","category":"page"},{"location":"special/#Ellipsoidal-harmonics","page":"Special functions","title":"Ellipsoidal harmonics","text":"","category":"section"},{"location":"special/#Orthogonal-polynomials","page":"Special functions","title":"Orthogonal polynomials","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"evaluate values\nroots","category":"page"},{"location":"special/#Hypergeometric-functions","page":"Special functions","title":"Hypergeometric functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] hyp2f1: cepheshyp2f1，hyp2f1complex\n[ ] hyp1f1: boost::hyp1f1double; chyp1f1wrap -> specfun_cchg\n[ ] hyperu: cephespoch,specfunchgu\n[ ] hyp0f1: cbesiwrap,cbesjwrap,hyp0f1asy","category":"page"},{"location":"special/#Parabolic-Cylinder-functions","page":"Special functions","title":"Parabolic Cylinder functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] pbdv: specfun_pbdv\n[ ] pbvv: specfun_pbvv\n[ ] pbwa: specfun_pbwa\n[ ] pbdv_seq: specfun.pbdv -> specfunpbdv\n[ ] pbvv_seq: specfun.pbvv -> specfunpbvv\n[ ] pbwa_seq: specfun.cpbdn -> specfuncpbdn","category":"page"},{"location":"special/#Mathieu-functions","page":"Special functions","title":"Mathieu functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] mathieu_a(m, q): cemcvawrap -> specfun_cva2\n[ ] mathieu_b(m, q): semcvawrap -> specfun_cva2\n[ ] specfun.fcoef -> specfunfcoef:\nmathieu_even_coef(m, q)\nmathieu_odd_coef(m, q)\n[ ] mathieu_cem(m, q, x): specfun_mtu0\n[ ] mathieu_sem(m, q, x): specfun_mtu0\n[ ] mathieu_modcem1(m, q, x): mcm1wrap -> specfunmtu12\n[ ] mathieu_modcem2(m, q, x): mcm2wrap -> specfunmtu12\n[ ] mathieu_modsem1(m, q, x): msm1wrap -> specfunmtu12\n[ ] mathieu_modsem2(m, q, x): msm2wrap -> specfunmtu12","category":"page"},{"location":"special/#Spheroidal-Wave-functions","page":"Special functions","title":"Spheroidal Wave functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] pro_ang1(m, n, c, x): prolateaswfanocvwrap -> specfunsegv,specfun_aswfa\n[ ] pro_rad1(m, n, c, x): prolateradial1nocvwrap -> specfunsegv(1),specfun_rswfp(1)\n[ ] pro_rad2(m, n, c, x): prolateradial2nocvwrap -> specfunsegv(1),specfun_rswfp(2)\n[ ] obl_ang1(m, n, c, x): oblateaswfanocvwrap -> specfunsegv,specfun_aswfa\n[ ] obl_rad1(m, n, c, x): oblateradial1nocvwrap -> specfunsegv,specfun_rswfo\n[ ] obl_rad2(m, n, c, x): oblateradial2nocvwrap -> specfunsegv,specfun_rswfo\n[ ] pro_cv(m, n, c): prolatesegvwrap -> specfun_segv\n[ ] obl_cv(m, n, c): oblatesegvwrap -> specfun_segv\n[ ] pro_cv_seq(m, n, c): specfun.segv -> specfunsegv\n[ ] obl_cv_seq(m, n, c): specfun.segv -> specfunsegv\n[ ] pro_ang1_cv(m, n, c, cv, x): prolateaswfawrap -> specfun_aswfa\n[ ] pro_rad1_cv(m, n, c, cv, x): prolateradial1wrap -> specfun_rswfp(1)\n[ ] pro_rad2_cv(m, n, c, cv, x): prolateradial2wrap -> specfun_rswfp(2)\n[ ] obl_ang1_cv(m, n, c, cv, x): oblateaswfawrap -> specfun_aswfa\n[ ] obl_rad1_cv(m, n, c, cv, x): oblateradial1wrap -> specfun_rswfo\n[ ] obl_rad2_cv(m, n, c, cv, x): oblateradial2wrap -> specfun_rswfo","category":"page"},{"location":"special/#Kelvin-functions","page":"Special functions","title":"Kelvin functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] kelvin(x): kelvinwrap -> specfunklvna\nber(x)\nbei(x)\nberp(x)\nbeip(x)\nker(x)\nkei(x)\nkerp(x)\nkeip(x)\n[ ] kelvin_zeros(nt): specfun_klvnzo","category":"page"},{"location":"special/#Combinatorics","page":"Special functions","title":"Combinatorics","text":"","category":"section"},{"location":"special/#Lambert-W-functions","page":"Special functions","title":"Lambert W functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"in cpp","category":"page"},{"location":"special/#Other-Special-functions","page":"Special functions","title":"Other Special functions","text":"","category":"section"},{"location":"special/","page":"Special functions","title":"Special functions","text":"[ ] agm(a, b): \n[ ] bernoulli(n): specfun.bernob -> specfunbernob\n[ ] binom(x, y): \n[ ] diric(x, n): \n[ ] euler(n): specfun.eulerb -> specfuneulerb\n[ ] expn(n, x): \n[ ] exp1(z): exp1wrap -> specfune1xb; cexp1wrap -> specfune1z\n[ ] expi(x): expiwrap -> specfuneix; cexpiwrap -> specfuneixz\n[ ] factorial(n): \n[ ] factorial2(n): \n[ ] factorialk(n, k): \n[ ] shichi(x): \n[ ] sici(x): \n[ ] softmax(x): \n[ ] log_softmax(x): \n[ ] spence(z): \n[ ] zeta(x): \n[ ] zetac(x): ","category":"page"},{"location":"special/#Convenience-functions","page":"Special functions","title":"Convenience functions","text":"","category":"section"},{"location":"specfun/#f77-specfun","page":"f77 specfun","title":"f77 specfun","text":"","category":"section"},{"location":"specfun/","page":"f77 specfun","title":"f77 specfun","text":"ref: test\\f77\\specfun.f","category":"page"},{"location":"specfun/","page":"f77 specfun","title":"f77 specfun","text":"155 SUBROUTINE\n5 FUNCTION","category":"page"},{"location":"specfun/#scipy-used-functions","page":"f77 specfun","title":"scipy used functions","text":"","category":"section"},{"location":"specfun/","page":"f77 specfun","title":"f77 specfun","text":"// https://github.com/scipy/scipy/blob/2b84173/scipy/special/specfun.h\n    specfun_airyzo\n    specfun_aswfa\n    specfun_bernob\n    specfun_cerror\n    specfun_cerzo\n    specfun_cchg\n    specfun_cfc\n    specfun_cfs\n    specfun_cgama\n    specfun_chgm\n    specfun_chgu\n    specfun_clpmn\n    specfun_clpn\n    specfun_clqmn\n    specfun_clqn\n    specfun_cpbdn\n    specfun_cva2\n    specfun_cyzo\n    specfun_eix\n    specfun_e1xb\n    specfun_eixz\n    specfun_e1z\n    specfun_eulerb\n    specfun_fcoef\n    specfun_fcszo\n    specfun_ffk\n    specfun_hygfz\n    specfun_itairy\n    specfun_itika\n    specfun_itjya\n    specfun_itsh0\n    specfun_itsl0\n    specfun_itth0\n    specfun_ittika\n    specfun_ittjya\n    specfun_jdzo\n    specfun_jyzo\n    specfun_klvna\n    specfun_klvnzo\n    specfun_lamn\n    specfun_lamv\n    specfun_lpmn\n    specfun_lpmv\n    specfun_lpn\n    specfun_lqmn\n    specfun_lqnb\n    specfun_mtu0\n    specfun_mtu12\n    specfun_pbdv\n    specfun_pbvv\n    specfun_pbwa\n    specfun_rctj\n    specfun_rcty\n    specfun_rswfp\n    specfun_rswfo\n    specfun_sdmn\n    specfun_segv","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Scipy4j","category":"page"},{"location":"#Scipy4j","page":"Home","title":"Scipy4j","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Scipy4j.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Scipy4j]","category":"page"},{"location":"#Scipy4j.ai_zeros-Tuple{Any}","page":"Home","title":"Scipy4j.ai_zeros","text":"ai_zeros(nt)\n\nCompute the first nt zeros of Airy functions Ai(x) and Ai'(x), a and a', and the associated values of Ai(a') and Ai'(a).\n\n\n\n\n\n","category":"method"},{"location":"#Scipy4j.bi_zeros-Tuple{Any}","page":"Home","title":"Scipy4j.bi_zeros","text":"bi_zeros(nt)\n\nCompute the first nt zeros of Airy functions Bi(x) and Bi'(x), b and b', and the associated values of Bi(b') and Bi'(b).\n\n\n\n\n\n","category":"method"},{"location":"#Scipy4j.itairy-Tuple{Any}","page":"Home","title":"Scipy4j.itairy","text":"itairy(x)\n\nCompute the integrals of Airy functions with respect to t from 0 and x ( x ≥ 0 )\n\n\n\n\n\n","category":"method"}]
}
