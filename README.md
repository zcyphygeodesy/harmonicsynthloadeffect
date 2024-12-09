## Fortran codes for spherical harmonic synthesis of all-element surface load effects
https://www.zcyphygeodesy.com/en/h-nd-130.html
## [Algorithm purpose]
    From the surface atmosphere, land water or sea level variation load spherical harmonic coefficient model (m), compute the non-tidal load effects on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (10μE) or horizontal gravity gradient (NW, to the north and to the west, 10μE) using spherical harmonic synthesis.
    The Algorithm can computate unified analytically various load effects on all-element geodetic variations in whole Earth space. The calculated surface load effect time time seriescan be employed to calibrate various parameters of the satellite's key geo-detic payloads, and effectively improve and check the quality, reliability and accuracy of the time-varying monitoring for satellite gravity field.
## [Computation Output]
    tdn(14): the non-tidal load effects on all-element geodetic variations.
    tdn(1:14) stores the non-tidal load effects on 10 kinds of geodetic variations, which are the load effects on height anomaly tdn(1) (mm), ground gravity #tdn(2) (μGal), gravity disturbance tdn(3) (μGal), ground tilt #tdn(4:5) (SW, to the south and to the west, mas), vertical deflection tdn(6:7) (SW, to the south and to the west, mas), horizontal displacement #tdn(8:9) (EN, to the east and to the north, mm), ground radial displacement #tdn(10) (mm), ground normal or orthometric height #tdn(11) (mm), radial gravity gradient tdn(12 )(10μE) or horizontal gravity gradient tdn(13:14) (NW, to the north and to the west, 10μE).
    The calculation point can be on the ground, low altitude, satellite, ocean or underwater space. The geodetic variations abvove marked with # are valid only when the site is fixed with the solid Earth.
## [Geophysical models]
    The Earth’s Load Love number file love_load_cm.dat from a Regional EIAstic Rebound calculator (REAR1.0, 2015).
## [Main program for test entrance]
    Harmsynthloadeffect.f90
    The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), height (m) relative to landsea surface......
    The record format of the output file reslt.txt: Behind the record of the calculation point file, appends 14 columns of load effects on the height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (10μE) or horizontal gravity gradient (NW, to the north and to the west, 10μE).
## (1) Algorithm module for spherical harmonic synthesis of all-element load effects
    Loadeffectpnm(rln,maxn,cnm,snm,flv,tdn,GRS,pnm,dpt1,dpt2,gr,hh)
    Input parameters: cnm, snm, maxn - the load direct influences of geopotential coefficients and maximum calculation degree
    Input parameters: pnm,dpt1,dpt2 - the normalized associative Legendre function and their first and second derivatives at the calculation point.
    Input parameters: rln(3), gr, hh – the spherical coordinates, normal gravity and height relative to landsea surface of the calculation point.
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value.
    When GRS(6) = -1, the atmospheric load effect is calculated. When GRS(6) > 0, the load effect of land water, sea level variation or the sum of the two is calculated.
## (2) Calculation module for normal Earth’s gravity field
    normdjn(GRS,djn); GNormalfd(BLH, NFD, GRS)
    Return parameters: NFD(5) - the normal geopotential (m²/s²), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass)..
## (3) Algorithm module for normalized associative Legendre function and its derivative
    BelPnmdt(pnm,dpt1,dpt2,maxn,t)
    The normalized associated Legendre functions and thier derivatives. Improved Belikov recursion algorithm for Pnm and Non-singular recursive algorithm for derivative of Pnm.
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm module for transforming ellipsoid geodetic coordinates into spherical coordinates
    BLH_RLAT(GRS,BLH,RLAT)
## (6) Other auxiliary modules
    PickRecord(str0, kln, rec, nn)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] ETideLoad4.5 User Reference https://www.zcyphygeodesy.com/en/
    8.2.2 The normalized spherical harmonic series expansion for surface load deformation field
    8.2.3 The normalized associated Legendre functions and thier derivatives
DOS executable test file, geophysical models and all input and output data.
