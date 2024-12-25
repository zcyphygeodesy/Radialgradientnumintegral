## Fortran codes for operation of radial gradient numerical integral on gravity field element
https://www.zcyphygeodesy.com/en/h-nd-148.html
## [Algorithm purpose]
    From the ellipsoidal height grid of the equipotential boundary surface and anomalous gravity field element grid on the surface, compute the radial gradient (/km) of the field element on the surface by the numerical integral.
    The radial gradient integral algorithm of the anomalous field element is derived from the solution of the Stokes boundary value problem, which requires the integrand field elements to be on the equipotential surface.
    The radial direction points from the center of the Earth to the outside of the Earth.
    The equipotential surface can be constructed from a global geopotential model (not greater than 360 degrees), which can also be represent by a normal (orthometric) equiheight surface with the altitude of not more than ten kilometers.
## [Main program for test entrance]
    Radialgradientnumintegral.f90
    Input parameters: dr - the integral radius (m).
    Input parameters: calcpntfl - the calculation point file name on the equipotential boundary surface. The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees),......
    Input parameters: dwmhgrdfl - the ellipsoidal height grid file name of the equipotential boundary surface. The grid will be employed to calculate the integral distance.
    Input parameters: gravgrdfl -  the residual field element grid file name on the equipotential surface.
## (1) Module for radial gradient operation on residual field element
    Radialgradientpnt(calcpntfl,dwmhgrdfl,gravgrdfl,dr)
    The output file reslt.txt, whose record format: Behind the input calculation point file record, appends a column of ellipsoidal height of the calculation point interpolated from the ellipsoidal height grid of the equipotential surface and a column of integral value of the radial gradient (in unit of /km).
## (2) Module for numerical integral of radial gradient of residual field element
    real*8 function RaGradientBLH(BLH,gra,dwm,nlat,nlon,hd,dr,GRS)
    Input parameters: BLH(3) - longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m) of the calculation point.
    Input parameters: dwm(nlat,nlon) - the ellipsoidal height grid of the equipotential boundary surface, which employed to calculate the integral distance.
    Input parameters: gra(nlat,nlon) - the residual field element grid on the equipotential surface.
    Input parameters: hd(6) - the grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid).
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Return - the calculated residual radial gradient (in unit of /m).
## (3) Calculation module for the normal gravity field
    GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Other auxiliary modules
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    CGrdPntD2(lon,lat,dt,row,col,hd); PickRecord(line, kln, rec, nn)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
The zip compression package in the attachment includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file and all input and output data.
