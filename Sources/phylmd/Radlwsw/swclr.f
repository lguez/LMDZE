module swclr_m

  IMPLICIT NONE

contains

  SUBROUTINE swclr(knu, flag_aer, palbp, pdsig, prayl, psec, pcgaz, ppizaz, &
       pray1, pray2, prefz, prj, prk, prmu0, ptauaz, ptra1, ptra2)
    
    USE raddim, only: kdlon, kflev
    USE radepsi, only: repsct, zepsec
    USE radopt, only: novlp

    ! ------------------------------------------------------------------
    ! PURPOSE.
    ! --------
    ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
    ! CLEAR-SKY COLUMN

    ! REFERENCE.
    ! ----------

    ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    ! AUTHOR.
    ! -------
    ! JEAN-JACQUES MORCRETTE  *ECMWF*

    ! MODIFICATIONS.
    ! --------------
    ! ORIGINAL : 94-11-15
    ! ------------------------------------------------------------------
    ! * ARGUMENTS:

    INTEGER knu
    ! -OB
    logical, intent(in):: flag_aer
    DOUBLE PRECISION palbp(kdlon, 2)
    DOUBLE PRECISION pdsig(kdlon, kflev)
    DOUBLE PRECISION prayl(kdlon)
    DOUBLE PRECISION psec(kdlon)

    DOUBLE PRECISION pcgaz(kdlon, kflev)
    DOUBLE PRECISION ppizaz(kdlon, kflev)
    DOUBLE PRECISION pray1(kdlon, kflev+1)
    DOUBLE PRECISION pray2(kdlon, kflev+1)
    DOUBLE PRECISION prefz(kdlon, 2, kflev+1)
    DOUBLE PRECISION prj(kdlon, 6, kflev+1)
    DOUBLE PRECISION prk(kdlon, 6, kflev+1)
    DOUBLE PRECISION prmu0(kdlon, kflev+1)
    DOUBLE PRECISION ptauaz(kdlon, kflev)
    DOUBLE PRECISION ptra1(kdlon, kflev+1)
    DOUBLE PRECISION ptra2(kdlon, kflev+1)

    ! * LOCAL VARIABLES:

    DOUBLE PRECISION zc0i(kdlon, kflev+1)
    DOUBLE PRECISION zclear(kdlon)
    DOUBLE PRECISION zr21(kdlon)
    DOUBLE PRECISION zss0(kdlon)
    DOUBLE PRECISION zscat(kdlon)
    DOUBLE PRECISION ztr(kdlon, 2, kflev+1)

    INTEGER jl, jk, ja, jkl, jklp1, jaj, jkm1
    DOUBLE PRECISION ztray, zgar, zratio, zff, zfacoa, zcorae
    DOUBLE PRECISION zmue, zgap, zww, zto, zden, zmu1, zden1
    DOUBLE PRECISION zbmu0, zbmu1, zre11

    ! ------------------------------------------------------------------

    ! *         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
    ! --------------------------------------------


    DO jk = 1, kflev + 1
       DO ja = 1, 6
          DO jl = 1, kdlon
             prj(jl, ja, jk) = 0.
             prk(jl, ja, jk) = 0.
          END DO
       END DO
    END DO

    DO jk = 1, kflev
       DO jl = 1, kdlon
          ptauaz(jl, jk) = 0d0
          ppizaz(jl, jk) = 0d0
          pcgaz(jl, jk) = 0d0
       END DO

       IF (flag_aer) THEN
          ! -OB
          DO jl = 1, kdlon
             ztray = prayl(jl)*pdsig(jl, jk)
             zratio = ztray/(ztray+ptauaz(jl,jk))
             zgar = pcgaz(jl, jk)
             zff = zgar*zgar
             ptauaz(jl, jk) = ztray + ptauaz(jl, jk)*(1.-ppizaz(jl,jk)*zff)
             pcgaz(jl, jk) = zgar*(1.-zratio)/(1.+zgar)
             ppizaz(jl, jk) = zratio + (1.-zratio)*ppizaz(jl, jk)*(1.-zff)/(1.- &
                  ppizaz(jl,jk)*zff)
          END DO
       ELSE
          DO jl = 1, kdlon
             ztray = prayl(jl)*pdsig(jl, jk)
             ptauaz(jl, jk) = ztray
             pcgaz(jl, jk) = 0.
             ppizaz(jl, jk) = 1. - repsct
          END DO
       END IF
    END DO

    ! ------------------------------------------------------------------

    ! *         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
    ! ----------------------------------------------


    DO jl = 1, kdlon
       zc0i(jl, kflev+1) = 0.
       zclear(jl) = 1.
       zscat(jl) = 0.
    END DO

    jk = 1
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
       zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
       zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
       zr21(jl) = exp(-zcorae)
       zss0(jl) = 1. - zr21(jl)

       IF (novlp==1) THEN
          ! * maximum-random
          zclear(jl) = zclear(jl)*(1.0-max(zss0(jl),zscat(jl)))/ &
               (1.0-min(zscat(jl),1.-zepsec))
          zc0i(jl, jkl) = 1.0 - zclear(jl)
          zscat(jl) = zss0(jl)
       ELSE IF (novlp==2) THEN
          ! * maximum
          zscat(jl) = max(zss0(jl), zscat(jl))
          zc0i(jl, jkl) = zscat(jl)
       ELSE IF (novlp==3) THEN
          ! * random
          zclear(jl) = zclear(jl)*(1.0-zss0(jl))
          zscat(jl) = 1.0 - zclear(jl)
          zc0i(jl, jkl) = zscat(jl)
       END IF
    END DO

    DO jk = 2, kflev
       jkl = kflev + 1 - jk
       jklp1 = jkl + 1
       DO jl = 1, kdlon
          zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
          zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
          zr21(jl) = exp(-zcorae)
          zss0(jl) = 1. - zr21(jl)

          IF (novlp==1) THEN
             ! * maximum-random
             zclear(jl) = zclear(jl)*(1.0-max(zss0(jl),zscat(jl)))/ &
                  (1.0-min(zscat(jl),1.-zepsec))
             zc0i(jl, jkl) = 1.0 - zclear(jl)
             zscat(jl) = zss0(jl)
          ELSE IF (novlp==2) THEN
             ! * maximum
             zscat(jl) = max(zss0(jl), zscat(jl))
             zc0i(jl, jkl) = zscat(jl)
          ELSE IF (novlp==3) THEN
             ! * random
             zclear(jl) = zclear(jl)*(1.0-zss0(jl))
             zscat(jl) = 1.0 - zclear(jl)
             zc0i(jl, jkl) = zscat(jl)
          END IF
       END DO
    END DO

    ! ------------------------------------------------------------------

    ! *         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
    ! -----------------------------------------------


    DO jl = 1, kdlon
       pray1(jl, kflev+1) = 0.
       pray2(jl, kflev+1) = 0.
       prefz(jl, 2, 1) = palbp(jl, knu)
       prefz(jl, 1, 1) = palbp(jl, knu)
       ptra1(jl, kflev+1) = 1.
       ptra2(jl, kflev+1) = 1.
    END DO

    DO jk = 2, kflev + 1
       jkm1 = jk - 1
       DO jl = 1, kdlon


          ! ------------------------------------------------------------------

          ! *         3.1  EQUIVALENT ZENITH ANGLE
          ! -----------------------


          zmue = (1.-zc0i(jl,jk))*psec(jl) + zc0i(jl, jk)*1.66
          prmu0(jl, jk) = 1./zmue


          ! ------------------------------------------------------------------

          ! *         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
          ! ----------------------------------------------------


          zgap = pcgaz(jl, jkm1)
          zbmu0 = 0.5 - 0.75*zgap/zmue
          zww = ppizaz(jl, jkm1)
          zto = ptauaz(jl, jkm1)
          zden = 1. + (1.-zww+zbmu0*zww)*zto*zmue + (1-zww)*(1.-zww+2.*zbmu0*zww) &
               *zto*zto*zmue*zmue
          pray1(jl, jkm1) = zbmu0*zww*zto*zmue/zden
          ptra1(jl, jkm1) = 1./zden

          zmu1 = 0.5
          zbmu1 = 0.5 - 0.75*zgap*zmu1
          zden1 = 1. + (1.-zww+zbmu1*zww)*zto/zmu1 + (1-zww)*(1.-zww+2.*zbmu1*zww &
               )*zto*zto/zmu1/zmu1
          pray2(jl, jkm1) = zbmu1*zww*zto/zmu1/zden1
          ptra2(jl, jkm1) = 1./zden1



          prefz(jl, 1, jk) = (pray1(jl,jkm1)+prefz(jl,1,jkm1)*ptra1(jl,jkm1)* &
               ptra2(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1,jkm1)))

          ztr(jl, 1, jkm1) = (ptra1(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1, &
               jkm1)))

          prefz(jl, 2, jk) = (pray1(jl,jkm1)+prefz(jl,2,jkm1)*ptra1(jl,jkm1)* &
               ptra2(jl,jkm1))

          ztr(jl, 2, jkm1) = ptra1(jl, jkm1)

       END DO
    END DO
    DO jl = 1, kdlon
       zmue = (1.-zc0i(jl,1))*psec(jl) + zc0i(jl, 1)*1.66
       prmu0(jl, 1) = 1./zmue
    END DO


    ! ------------------------------------------------------------------

    ! *         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
    ! -------------------------------------------------


    IF (knu==1) THEN
       jaj = 2
       DO jl = 1, kdlon
          prj(jl, jaj, kflev+1) = 1.
          prk(jl, jaj, kflev+1) = prefz(jl, 1, kflev+1)
       END DO

       DO jk = 1, kflev
          jkl = kflev + 1 - jk
          jklp1 = jkl + 1
          DO jl = 1, kdlon
             zre11 = prj(jl, jaj, jklp1)*ztr(jl, 1, jkl)
             prj(jl, jaj, jkl) = zre11
             prk(jl, jaj, jkl) = zre11*prefz(jl, 1, jkl)
          END DO
       END DO

    ELSE

       DO jaj = 1, 2
          DO jl = 1, kdlon
             prj(jl, jaj, kflev+1) = 1.
             prk(jl, jaj, kflev+1) = prefz(jl, jaj, kflev+1)
          END DO

          DO jk = 1, kflev
             jkl = kflev + 1 - jk
             jklp1 = jkl + 1
             DO jl = 1, kdlon
                zre11 = prj(jl, jaj, jklp1)*ztr(jl, jaj, jkl)
                prj(jl, jaj, jkl) = zre11
                prk(jl, jaj, jkl) = zre11*prefz(jl, jaj, jkl)
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE swclr

end module swclr_m
