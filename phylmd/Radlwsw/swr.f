module swr_m

  IMPLICIT NONE

contains

  SUBROUTINE swr(knu, palbd, pcg, pcld, pomega, psec, ptau, &
       pcgaz, ppizaz, pray1, pray2, prefz, prj, prk, prmue, ptauaz, ptra1, &
       ptra2)
    USE dimens_m
    USE dimphy
    USE raddim
    USE radepsi
    USE radopt
    use swde_m, only: swde

    ! ------------------------------------------------------------------
    ! PURPOSE.
    ! --------
    ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
    ! CONTINUUM SCATTERING

    ! METHOD.
    ! -------

    ! 1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
    ! OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)

    ! REFERENCE.
    ! ----------

    ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    ! AUTHOR.
    ! -------
    ! JEAN-JACQUES MORCRETTE  *ECMWF*

    ! MODIFICATIONS.
    ! --------------
    ! ORIGINAL : 89-07-14
    ! ------------------------------------------------------------------
    ! * ARGUMENTS:

    INTEGER knu
    DOUBLE PRECISION palbd(kdlon, 2)
    DOUBLE PRECISION pcg(kdlon, 2, kflev)
    DOUBLE PRECISION pcld(kdlon, kflev)
    DOUBLE PRECISION pomega(kdlon, 2, kflev)
    DOUBLE PRECISION psec(kdlon)
    DOUBLE PRECISION ptau(kdlon, 2, kflev)

    DOUBLE PRECISION pray1(kdlon, kflev+1)
    DOUBLE PRECISION pray2(kdlon, kflev+1)
    DOUBLE PRECISION prefz(kdlon, 2, kflev+1)
    DOUBLE PRECISION prj(kdlon, 6, kflev+1)
    DOUBLE PRECISION prk(kdlon, 6, kflev+1)
    DOUBLE PRECISION prmue(kdlon, kflev+1)
    DOUBLE PRECISION pcgaz(kdlon, kflev)
    DOUBLE PRECISION ppizaz(kdlon, kflev)
    DOUBLE PRECISION ptauaz(kdlon, kflev)
    DOUBLE PRECISION ptra1(kdlon, kflev+1)
    DOUBLE PRECISION ptra2(kdlon, kflev+1)

    ! * LOCAL VARIABLES:

    DOUBLE PRECISION zc1i(kdlon, kflev+1)
    DOUBLE PRECISION zclear(kdlon)
    DOUBLE PRECISION zcloud(kdlon)
    DOUBLE PRECISION zgg(kdlon)
    DOUBLE PRECISION zref(kdlon)
    DOUBLE PRECISION zre1(kdlon)
    DOUBLE PRECISION zre2(kdlon)
    DOUBLE PRECISION zrmuz(kdlon)
    DOUBLE PRECISION zrneb(kdlon)
    DOUBLE PRECISION zr21(kdlon)
    DOUBLE PRECISION zr22(kdlon)
    DOUBLE PRECISION zss1(kdlon)
    DOUBLE PRECISION zto1(kdlon)
    DOUBLE PRECISION ztr(kdlon, 2, kflev+1)
    DOUBLE PRECISION ztr1(kdlon)
    DOUBLE PRECISION ztr2(kdlon)
    DOUBLE PRECISION zw(kdlon)

    INTEGER jk, jl, ja, jkl, jklp1, jkm1, jaj
    DOUBLE PRECISION zfacoa, zfacoc, zcorae, zcorcd
    DOUBLE PRECISION zmue, zgap, zww, zto, zden, zden1
    DOUBLE PRECISION zmu1, zre11, zbmu0, zbmu1

    ! ------------------------------------------------------------------

    ! *         1.    INITIALIZATION
    ! --------------


    DO jk = 1, kflev + 1
       DO ja = 1, 6
          DO jl = 1, kdlon
             prj(jl, ja, jk) = 0.
             prk(jl, ja, jk) = 0.
          END DO
       END DO
    END DO


    ! ------------------------------------------------------------------

    ! *         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
    ! ----------------------------------------------


    DO jl = 1, kdlon
       zc1i(jl, kflev+1) = 0.
       zclear(jl) = 1.
       zcloud(jl) = 0.
    END DO

    jk = 1
    jkl = kflev + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, kdlon
       zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
       zfacoc = 1. - pomega(jl, knu, jkl)*pcg(jl, knu, jkl)*pcg(jl, knu, jkl)
       zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
       zcorcd = zfacoc*ptau(jl, knu, jkl)*psec(jl)
       zr21(jl) = exp(-zcorae)
       zr22(jl) = exp(-zcorcd)
       zss1(jl) = pcld(jl, jkl)*(1.0-zr21(jl)*zr22(jl)) + &
            (1.0-pcld(jl,jkl))*(1.0-zr21(jl))

       IF (novlp==1) THEN
          ! * maximum-random
          zclear(jl) = zclear(jl)*(1.0-max(zss1(jl),zcloud(jl)))/ &
               (1.0-min(zcloud(jl),1.-zepsec))
          zc1i(jl, jkl) = 1.0 - zclear(jl)
          zcloud(jl) = zss1(jl)
       ELSE IF (novlp==2) THEN
          ! * maximum
          zcloud(jl) = max(zss1(jl), zcloud(jl))
          zc1i(jl, jkl) = zcloud(jl)
       ELSE IF (novlp==3) THEN
          ! * random
          zclear(jl) = zclear(jl)*(1.0-zss1(jl))
          zcloud(jl) = 1.0 - zclear(jl)
          zc1i(jl, jkl) = zcloud(jl)
       END IF
    END DO

    DO jk = 2, kflev
       jkl = kflev + 1 - jk
       jklp1 = jkl + 1
       DO jl = 1, kdlon
          zfacoa = 1. - ppizaz(jl, jkl)*pcgaz(jl, jkl)*pcgaz(jl, jkl)
          zfacoc = 1. - pomega(jl, knu, jkl)*pcg(jl, knu, jkl)*pcg(jl, knu, jkl)
          zcorae = zfacoa*ptauaz(jl, jkl)*psec(jl)
          zcorcd = zfacoc*ptau(jl, knu, jkl)*psec(jl)
          zr21(jl) = exp(-zcorae)
          zr22(jl) = exp(-zcorcd)
          zss1(jl) = pcld(jl, jkl)*(1.0-zr21(jl)*zr22(jl)) + &
               (1.0-pcld(jl,jkl))*(1.0-zr21(jl))

          IF (novlp==1) THEN
             ! * maximum-random
             zclear(jl) = zclear(jl)*(1.0-max(zss1(jl),zcloud(jl)))/ &
                  (1.0-min(zcloud(jl),1.-zepsec))
             zc1i(jl, jkl) = 1.0 - zclear(jl)
             zcloud(jl) = zss1(jl)
          ELSE IF (novlp==2) THEN
             ! * maximum
             zcloud(jl) = max(zss1(jl), zcloud(jl))
             zc1i(jl, jkl) = zcloud(jl)
          ELSE IF (novlp==3) THEN
             ! * random
             zclear(jl) = zclear(jl)*(1.0-zss1(jl))
             zcloud(jl) = 1.0 - zclear(jl)
             zc1i(jl, jkl) = zcloud(jl)
          END IF
       END DO
    END DO

    ! ------------------------------------------------------------------

    ! *         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
    ! -----------------------------------------------


    DO jl = 1, kdlon
       pray1(jl, kflev+1) = 0.
       pray2(jl, kflev+1) = 0.
       prefz(jl, 2, 1) = palbd(jl, knu)
       prefz(jl, 1, 1) = palbd(jl, knu)
       ptra1(jl, kflev+1) = 1.
       ptra2(jl, kflev+1) = 1.
    END DO

    DO jk = 2, kflev + 1
       jkm1 = jk - 1
       DO jl = 1, kdlon
          zrneb(jl) = pcld(jl, jkm1)
          zre1(jl) = 0.
          ztr1(jl) = 0.
          zre2(jl) = 0.
          ztr2(jl) = 0.


          ! ------------------------------------------------------------------

          ! *         3.1  EQUIVALENT ZENITH ANGLE
          ! -----------------------


          zmue = (1.-zc1i(jl,jk))*psec(jl) + zc1i(jl, jk)*1.66
          prmue(jl, jk) = 1./zmue


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
          ! PRINT *,' LOOP 342 ** 3 ** JL=',JL,PRAY1(JL,JKM1),PTRA1(JL,JKM1)

          zmu1 = 0.5
          zbmu1 = 0.5 - 0.75*zgap*zmu1
          zden1 = 1. + (1.-zww+zbmu1*zww)*zto/zmu1 + (1-zww)*(1.-zww+2.*zbmu1*zww &
               )*zto*zto/zmu1/zmu1
          pray2(jl, jkm1) = zbmu1*zww*zto/zmu1/zden1
          ptra2(jl, jkm1) = 1./zden1


          ! ------------------------------------------------------------------

          ! *         3.3  EFFECT OF CLOUD LAYER
          ! ---------------------


          zw(jl) = pomega(jl, knu, jkm1)
          zto1(jl) = ptau(jl, knu, jkm1)/zw(jl) + ptauaz(jl, jkm1)/ppizaz(jl, &
               jkm1)
          zr21(jl) = ptau(jl, knu, jkm1) + ptauaz(jl, jkm1)
          zr22(jl) = ptau(jl, knu, jkm1)/zr21(jl)
          zgg(jl) = zr22(jl)*pcg(jl, knu, jkm1) + (1.-zr22(jl))*pcgaz(jl, jkm1)
          ! Modif PhD - JJM 19/03/96 pour erreurs arrondis
          ! machine
          ! PHD PROTECTION ZW(JL) = ZR21(JL) / ZTO1(JL)
          IF (zw(jl)==1. .AND. ppizaz(jl,jkm1)==1.) THEN
             zw(jl) = 1.
          ELSE
             zw(jl) = zr21(jl)/zto1(jl)
          END IF
          zref(jl) = prefz(jl, 1, jkm1)
          zrmuz(jl) = prmue(jl, jk)
       END DO

       CALL swde(zgg, zref, zrmuz, zto1, zw, zre1, zre2, ztr1, ztr2)

       DO jl = 1, kdlon

          prefz(jl, 1, jk) = (1.-zrneb(jl))*(pray1(jl,jkm1)+prefz(jl,1,jkm1)* &
               ptra1(jl,jkm1)*ptra2(jl,jkm1)/(1.-pray2(jl,jkm1)*prefz(jl,1, &
               jkm1))) + zrneb(jl)*zre2(jl)

          ztr(jl, 1, jkm1) = zrneb(jl)*ztr2(jl) + (ptra1(jl,jkm1)/(1.-pray2(jl, &
               jkm1)*prefz(jl,1,jkm1)))*(1.-zrneb(jl))

          prefz(jl, 2, jk) = (1.-zrneb(jl))*(pray1(jl,jkm1)+prefz(jl,2,jkm1)* &
               ptra1(jl,jkm1)*ptra2(jl,jkm1)) + zrneb(jl)*zre1(jl)

          ztr(jl, 2, jkm1) = zrneb(jl)*ztr1(jl) + ptra1(jl, jkm1)*(1.-zrneb(jl))

       END DO
    END DO
    DO jl = 1, kdlon
       zmue = (1.-zc1i(jl,1))*psec(jl) + zc1i(jl, 1)*1.66
       prmue(jl, 1) = 1./zmue
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

  END SUBROUTINE swr

end module swr_m
