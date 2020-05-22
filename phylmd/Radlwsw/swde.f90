module swde_m

  IMPLICIT NONE

contains

  SUBROUTINE swde(pgg, pref, prmuz, pto1, pw, pre1, pre2, ptr1, ptr2)
    
    USE conf_phys_m, only: kdlon

    ! PURPOSE.
    ! --------
    ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
    ! LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.

    ! METHOD.
    ! -------

    ! STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.

    ! REFERENCE.
    ! ----------

    ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    ! AUTHOR.
    ! -------
    ! JEAN-JACQUES MORCRETTE  *ECMWF*

    ! MODIFICATIONS.
    ! --------------
    ! ORIGINAL : 88-12-15
    ! ------------------------------------------------------------------
    ! * ARGUMENTS:

    DOUBLE PRECISION pgg(kdlon) ! ASSYMETRY FACTOR
    DOUBLE PRECISION pref(kdlon) ! REFLECTIVITY OF THE UNDERLYING LAYER
    DOUBLE PRECISION prmuz(kdlon) ! COSINE OF SOLAR ZENITH ANGLE
    DOUBLE PRECISION pto1(kdlon) ! OPTICAL THICKNESS
    DOUBLE PRECISION pw(kdlon) ! SINGLE SCATTERING ALBEDO
    DOUBLE PRECISION pre1(kdlon) ! LAYER REFLECTIVITY (NO UNDERLYING-LAYER REFLECTION)
    DOUBLE PRECISION pre2(kdlon) ! LAYER REFLECTIVITY
    DOUBLE PRECISION ptr1(kdlon) ! LAYER TRANSMISSIVITY (NO UNDERLYING-LAYER REFLECTION)
    DOUBLE PRECISION ptr2(kdlon) ! LAYER TRANSMISSIVITY

    ! * LOCAL VARIABLES:

    INTEGER jl
    DOUBLE PRECISION zff, zgp, ztop, zwcp, zdt, zx1, zwm
    DOUBLE PRECISION zrm2, zrk, zx2, zrp, zalpha, zbeta, zarg
    DOUBLE PRECISION zexmu0, zarg2, zexkp, zexkm, zxp2p, zxm2p, zap2b
    DOUBLE PRECISION zam2b
    DOUBLE PRECISION za11, za12, za13, za21, za22, za23
    DOUBLE PRECISION zdena, zc1a, zc2a, zri0a, zri1a
    DOUBLE PRECISION zri0b, zri1b
    DOUBLE PRECISION zb21, zb22, zb23, zdenb, zc1b, zc2b
    DOUBLE PRECISION zri0c, zri1c, zri0d, zri1d
    ! ------------------------------------------------------------------

    ! *         1.      DELTA-EDDINGTON CALCULATIONS


    DO jl = 1, kdlon

       ! *         1.1     SET UP THE DELTA-MODIFIED PARAMETERS


       zff = pgg(jl)*pgg(jl)
       zgp = pgg(jl)/(1.+pgg(jl))
       ztop = (1.-pw(jl)*zff)*pto1(jl)
       zwcp = (1-zff)*pw(jl)/(1.-pw(jl)*zff)
       zdt = 2./3.
       zx1 = 1. - zwcp*zgp
       zwm = 1. - zwcp
       zrm2 = prmuz(jl)*prmuz(jl)
       zrk = sqrt(3.*zwm*zx1)
       zx2 = 4.*(1.-zrk*zrk*zrm2)
       zrp = zrk/zx1
       zalpha = 3.*zwcp*zrm2*(1.+zgp*zwm)/zx2
       zbeta = 3.*zwcp*prmuz(jl)*(1.+3.*zgp*zrm2*zwm)/zx2
       ! MAF      ZARG=MIN(ZTOP/PRMUZ(JL),200.)
       zarg = min(ztop/prmuz(jl), 2.0D+2)
       zexmu0 = exp(-zarg)
       ! MAF      ZARG2=MIN(ZRK*ZTOP,200.)
       zarg2 = min(zrk*ztop, 2.0D+2)
       zexkp = exp(zarg2)
       zexkm = 1./zexkp
       zxp2p = 1. + zdt*zrp
       zxm2p = 1. - zdt*zrp
       zap2b = zalpha + zdt*zbeta
       zam2b = zalpha - zdt*zbeta

       ! *         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER


       za11 = zxp2p
       za12 = zxm2p
       za13 = zap2b
       za22 = zxp2p*zexkp
       za21 = zxm2p*zexkm
       za23 = zam2b*zexmu0
       zdena = za11*za22 - za21*za12
       zc1a = (za22*za13-za12*za23)/zdena
       zc2a = (za11*za23-za21*za13)/zdena
       zri0a = zc1a + zc2a - zalpha
       zri1a = zrp*(zc1a-zc2a) - zbeta
       pre1(jl) = (zri0a-zdt*zri1a)/prmuz(jl)
       zri0b = zc1a*zexkm + zc2a*zexkp - zalpha*zexmu0
       zri1b = zrp*(zc1a*zexkm-zc2a*zexkp) - zbeta*zexmu0
       ptr1(jl) = zexmu0 + (zri0b+zdt*zri1b)/prmuz(jl)

       ! *         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER


       zb21 = za21 - pref(jl)*zxp2p*zexkm
       zb22 = za22 - pref(jl)*zxm2p*zexkp
       zb23 = za23 - pref(jl)*zexmu0*(zap2b-prmuz(jl))
       zdenb = za11*zb22 - zb21*za12
       zc1b = (zb22*za13-za12*zb23)/zdenb
       zc2b = (za11*zb23-zb21*za13)/zdenb
       zri0c = zc1b + zc2b - zalpha
       zri1c = zrp*(zc1b-zc2b) - zbeta
       pre2(jl) = (zri0c-zdt*zri1c)/prmuz(jl)
       zri0d = zc1b*zexkm + zc2b*zexkp - zalpha*zexmu0
       zri1d = zrp*(zc1b*zexkm-zc2b*zexkp) - zbeta*zexmu0
       ptr2(jl) = zexmu0 + (zri0d+zdt*zri1d)/prmuz(jl)
    END DO

  END SUBROUTINE swde

end module swde_m
