module sw2s_m

  IMPLICIT NONE

contains

  SUBROUTINE sw2s(knu, paki, palbd, palbp, pcg, pcld, pclear, pdsig, pomega, &
       poz, prmu, psec, ptau, pud, pwv, pqs, pfdown, pfup)
    
    USE dimensions
    USE dimphy
    USE raddim
    USE radepsi
    use swclr_m, only: swclr
    use swde_m, only: swde
    use swr_m, only: swr

    ! ------------------------------------------------------------------
    ! PURPOSE.
    ! --------

    ! THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
    ! SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).

    ! METHOD.
    ! -------

    ! 1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
    ! CONTINUUM SCATTERING
    ! 2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
    ! A GREY MOLECULAR ABSORPTION
    ! 3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
    ! OF ABSORBERS
    ! 4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
    ! 5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

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
    ! 94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
    ! ------------------------------------------------------------------
    ! * ARGUMENTS:

    INTEGER knu
    DOUBLE PRECISION paki(kdlon, 2)
    DOUBLE PRECISION palbd(kdlon, 2)
    DOUBLE PRECISION palbp(kdlon, 2)
    DOUBLE PRECISION pcg(kdlon, 2, kflev)
    DOUBLE PRECISION pcld(kdlon, kflev)
    DOUBLE PRECISION pclear(kdlon)
    DOUBLE PRECISION pdsig(kdlon, kflev)
    DOUBLE PRECISION pomega(kdlon, 2, kflev)
    DOUBLE PRECISION poz(kdlon, kflev)
    DOUBLE PRECISION pqs(kdlon, kflev)
    DOUBLE PRECISION prmu(kdlon)
    DOUBLE PRECISION psec(kdlon)
    DOUBLE PRECISION ptau(kdlon, 2, kflev)
    DOUBLE PRECISION pud(kdlon, 5, kflev+1)
    DOUBLE PRECISION pwv(kdlon, kflev)

    DOUBLE PRECISION pfdown(kdlon, kflev+1)
    DOUBLE PRECISION pfup(kdlon, kflev+1)

    ! * LOCAL VARIABLES:

    INTEGER iind2(2), iind3(3)
    DOUBLE PRECISION zcgaz(kdlon, kflev)
    DOUBLE PRECISION zfd(kdlon, kflev+1)
    DOUBLE PRECISION zfu(kdlon, kflev+1)
    DOUBLE PRECISION zg(kdlon)
    DOUBLE PRECISION zgg(kdlon)
    DOUBLE PRECISION zpizaz(kdlon, kflev)
    DOUBLE PRECISION zrayl(kdlon)
    DOUBLE PRECISION zray1(kdlon, kflev+1)
    DOUBLE PRECISION zray2(kdlon, kflev+1)
    DOUBLE PRECISION zref(kdlon)
    DOUBLE PRECISION zrefz(kdlon, 2, kflev+1)
    DOUBLE PRECISION zre1(kdlon)
    DOUBLE PRECISION zre2(kdlon)
    DOUBLE PRECISION zrj(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrj0(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrk(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrk0(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrl(kdlon, 8)
    DOUBLE PRECISION zrmue(kdlon, kflev+1)
    DOUBLE PRECISION zrmu0(kdlon, kflev+1)
    DOUBLE PRECISION zrmuz(kdlon)
    DOUBLE PRECISION zrneb(kdlon)
    DOUBLE PRECISION zr1(kdlon)
    DOUBLE PRECISION zr2(kdlon, 2)
    DOUBLE PRECISION zr3(kdlon, 3)
    DOUBLE PRECISION zr4(kdlon)
    DOUBLE PRECISION zr21(kdlon)
    DOUBLE PRECISION zr22(kdlon)
    DOUBLE PRECISION zs(kdlon)
    DOUBLE PRECISION ztauaz(kdlon, kflev)
    DOUBLE PRECISION zto1(kdlon)
    DOUBLE PRECISION ztr(kdlon, 2, kflev+1)
    DOUBLE PRECISION ztra1(kdlon, kflev+1)
    DOUBLE PRECISION ztra2(kdlon, kflev+1)
    DOUBLE PRECISION ztr1(kdlon)
    DOUBLE PRECISION ztr2(kdlon)
    DOUBLE PRECISION zw(kdlon)
    DOUBLE PRECISION zw1(kdlon)
    DOUBLE PRECISION zw2(kdlon, 2)
    DOUBLE PRECISION zw3(kdlon, 3)
    DOUBLE PRECISION zw4(kdlon)
    DOUBLE PRECISION zw5(kdlon)

    INTEGER jl, jk, k, jaj, ikm1, ikl, jn, jabs, jkm1
    INTEGER jref, jkl, jklp1, jajp, jkki, jkkp4, jn2j, iabs
    DOUBLE PRECISION zrmum1, zwh2o, zcneb, zaa, zbb, zrki, zre11

    ! * Prescribed Data:

    DOUBLE PRECISION rsun(2)
    SAVE rsun
    DOUBLE PRECISION rray(2, 6)
    SAVE rray
    DATA rsun(1)/0.441676d0/
    DATA rsun(2)/0.558324d0/
    DATA (rray(1,k), k=1, 6)/.428937d-01, .890743d+00, -.288555d+01, &
         .522744d+01, -.469173d+01, .161645d+01/
    DATA (rray(2,k), k=1, 6)/.697200d-02, .173297d-01, -.850903d-01, &
         .248261d+00, -.302031d+00, .129662d+00/

    ! ------------------------------------------------------------------

    ! *         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
    ! -------------------------------------------



    ! *         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
    ! -----------------------------------------


    DO jl = 1, kdlon
       zrmum1 = 1. - prmu(jl)
       zrayl(jl) = rray(knu, 1) + zrmum1*(rray(knu,2)+zrmum1*(rray(knu, &
            3)+zrmum1*(rray(knu,4)+zrmum1*(rray(knu,5)+zrmum1*rray(knu,6)))))
    END DO


    ! ------------------------------------------------------------------

    ! *         2.    CONTINUUM SCATTERING CALCULATIONS
    ! ---------------------------------


    ! *         2.1   CLEAR-SKY FRACTION OF THE COLUMN
    ! --------------------------------


    CALL swclr(knu, palbp, pdsig, zrayl, psec, zpizaz, zray1, zray2, zrefz, &
         zrj0, zrk0, zrmu0, ztauaz, ztra1, ztra2)


    ! *         2.2   CLOUDY FRACTION OF THE COLUMN
    ! -----------------------------


    zcgaz = 0d0
    CALL swr(knu, palbd, pcg, pcld, pomega, psec, ptau, zcgaz, &
         zpizaz, zray1, zray2, zrefz, zrj, zrk, zrmue, ztauaz, ztra1, ztra2)


    ! ------------------------------------------------------------------

    ! *         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
    ! ------------------------------------------------------


    jn = 2

    DO jabs = 1, 2


       ! *         3.1  SURFACE CONDITIONS
       ! ------------------


       DO jl = 1, kdlon
          zrefz(jl, 2, 1) = palbd(jl, knu)
          zrefz(jl, 1, 1) = palbd(jl, knu)
       END DO


       ! *         3.2  INTRODUCING CLOUD EFFECTS
       ! -------------------------


       DO jk = 2, kflev + 1
          jkm1 = jk - 1
          ikl = kflev + 1 - jkm1
          DO jl = 1, kdlon
             zrneb(jl) = pcld(jl, jkm1)
             IF (jabs==1 .AND. zrneb(jl)>2.*zeelog) THEN
                zwh2o = max(pwv(jl,jkm1), zeelog)
                zcneb = max(zeelog, min(zrneb(jl),1.-zeelog))
                zbb = pud(jl, jabs, jkm1)*pqs(jl, jkm1)/zwh2o
                zaa = max((pud(jl,jabs,jkm1)-zcneb*zbb)/(1.-zcneb), zeelog)
             ELSE
                zaa = pud(jl, jabs, jkm1)
                zbb = zaa
             END IF
             zrki = paki(jl, jabs)
             zs(jl) = exp(-zrki*zaa*1.66)
             zg(jl) = exp(-zrki*zaa/zrmue(jl,jk))
             ztr1(jl) = 0.
             zre1(jl) = 0.
             ztr2(jl) = 0.
             zre2(jl) = 0.

             zw(jl) = pomega(jl, knu, jkm1)
             zto1(jl) = ptau(jl, knu, jkm1)/zw(jl) + ztauaz(jl, jkm1)/zpizaz(jl, &
                  jkm1) + zbb*zrki

             zr21(jl) = ptau(jl, knu, jkm1) + ztauaz(jl, jkm1)
             zr22(jl) = ptau(jl, knu, jkm1)/zr21(jl)
             zgg(jl) = zr22(jl)*pcg(jl, knu, jkm1) + (1.-zr22(jl))*zcgaz(jl, jkm1)
             zw(jl) = zr21(jl)/zto1(jl)
             zref(jl) = zrefz(jl, 1, jkm1)
             zrmuz(jl) = zrmue(jl, jk)
          END DO

          CALL swde(zgg, zref, zrmuz, zto1, zw, zre1, zre2, ztr1, ztr2)

          DO jl = 1, kdlon

             zrefz(jl, 2, jk) = (1.-zrneb(jl))*(zray1(jl,jkm1)+zrefz(jl,2,jkm1)* &
                  ztra1(jl,jkm1)*ztra2(jl,jkm1))*zg(jl)*zs(jl) + zrneb(jl)*zre1(jl)

             ztr(jl, 2, jkm1) = zrneb(jl)*ztr1(jl) + (ztra1(jl,jkm1))*zg(jl)*(1.- &
                  zrneb(jl))

             zrefz(jl, 1, jk) = (1.-zrneb(jl))*(zray1(jl,jkm1)+zrefz(jl,1,jkm1)* &
                  ztra1(jl,jkm1)*ztra2(jl,jkm1)/(1.-zray2(jl,jkm1)*zrefz(jl,1, &
                  jkm1)))*zg(jl)*zs(jl) + zrneb(jl)*zre2(jl)

             ztr(jl, 1, jkm1) = zrneb(jl)*ztr2(jl) + (ztra1(jl,jkm1)/(1.-zray2(jl, &
                  jkm1)*zrefz(jl,1,jkm1)))*zg(jl)*(1.-zrneb(jl))

          END DO
       END DO

       ! *         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
       ! -------------------------------------------------


       DO jref = 1, 2

          jn = jn + 1

          DO jl = 1, kdlon
             zrj(jl, jn, kflev+1) = 1.
             zrk(jl, jn, kflev+1) = zrefz(jl, jref, kflev+1)
          END DO

          DO jk = 1, kflev
             jkl = kflev + 1 - jk
             jklp1 = jkl + 1
             DO jl = 1, kdlon
                zre11 = zrj(jl, jn, jklp1)*ztr(jl, jref, jkl)
                zrj(jl, jn, jkl) = zre11
                zrk(jl, jn, jkl) = zre11*zrefz(jl, jref, jkl)
             END DO
          END DO
       END DO
    END DO


    ! ------------------------------------------------------------------

    ! *         4.    INVERT GREY AND CONTINUUM FLUXES
    ! --------------------------------



    ! *         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
    ! ---------------------------------------------


    DO jk = 1, kflev + 1
       DO jaj = 1, 5, 2
          jajp = jaj + 1
          DO jl = 1, kdlon
             zrj(jl, jaj, jk) = zrj(jl, jaj, jk) - zrj(jl, jajp, jk)
             zrk(jl, jaj, jk) = zrk(jl, jaj, jk) - zrk(jl, jajp, jk)
             zrj(jl, jaj, jk) = max(zrj(jl,jaj,jk), zeelog)
             zrk(jl, jaj, jk) = max(zrk(jl,jaj,jk), zeelog)
          END DO
       END DO
    END DO

    DO jk = 1, kflev + 1
       DO jaj = 2, 6, 2
          DO jl = 1, kdlon
             zrj(jl, jaj, jk) = max(zrj(jl,jaj,jk), zeelog)
             zrk(jl, jaj, jk) = max(zrk(jl,jaj,jk), zeelog)
          END DO
       END DO
    END DO

    ! *         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
    ! ---------------------------------------------


    DO jk = 1, kflev + 1
       jkki = 1
       DO jaj = 1, 2
          iind2(1) = jaj
          iind2(2) = jaj
          DO jn = 1, 2
             jn2j = jn + 2*jaj
             jkkp4 = jkki + 4

             ! *         4.2.1  EFFECTIVE ABSORBER AMOUNTS
             ! --------------------------


             DO jl = 1, kdlon
                zw2(jl, 1) = log(zrj(jl,jn,jk)/zrj(jl,jn2j,jk))/paki(jl, jaj)
                zw2(jl, 2) = log(zrk(jl,jn,jk)/zrk(jl,jn2j,jk))/paki(jl, jaj)
             END DO

             ! *         4.2.2  TRANSMISSION FUNCTION
             ! ---------------------


             CALL swtt1(knu, 2, iind2, zw2, zr2)

             DO jl = 1, kdlon
                zrl(jl, jkki) = zr2(jl, 1)
                zrl(jl, jkkp4) = zr2(jl, 2)
             END DO

             jkki = jkki + 1
          END DO
       END DO

       ! *         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
       ! ------------------------------------------------------


       DO jl = 1, kdlon
          pfdown(jl, jk) = zrj(jl, 1, jk)*zrl(jl, 1)*zrl(jl, 3) + &
               zrj(jl, 2, jk)*zrl(jl, 2)*zrl(jl, 4)
          pfup(jl, jk) = zrk(jl, 1, jk)*zrl(jl, 5)*zrl(jl, 7) + &
               zrk(jl, 2, jk)*zrl(jl, 6)*zrl(jl, 8)
       END DO
    END DO


    ! ------------------------------------------------------------------

    ! *         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
    ! ----------------------------------------



    ! *         5.1   DOWNWARD FLUXES
    ! ---------------


    jaj = 2
    iind3(1) = 1
    iind3(2) = 2
    iind3(3) = 3

    DO jl = 1, kdlon
       zw3(jl, 1) = 0.
       zw3(jl, 2) = 0.
       zw3(jl, 3) = 0.
       zw4(jl) = 0.
       zw5(jl) = 0.
       zr4(jl) = 1.
       zfd(jl, kflev+1) = zrj0(jl, jaj, kflev+1)
    END DO
    DO jk = 1, kflev
       ikl = kflev + 1 - jk
       DO jl = 1, kdlon
          zw3(jl, 1) = zw3(jl, 1) + pud(jl, 1, ikl)/zrmu0(jl, ikl)
          zw3(jl, 2) = zw3(jl, 2) + pud(jl, 2, ikl)/zrmu0(jl, ikl)
          zw3(jl, 3) = zw3(jl, 3) + poz(jl, ikl)/zrmu0(jl, ikl)
          zw4(jl) = zw4(jl) + pud(jl, 4, ikl)/zrmu0(jl, ikl)
          zw5(jl) = zw5(jl) + pud(jl, 5, ikl)/zrmu0(jl, ikl)
       END DO

       CALL swtt1(knu, 3, iind3, zw3, zr3)

       DO jl = 1, kdlon
          ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
          zfd(jl, ikl) = zr3(jl, 1)*zr3(jl, 2)*zr3(jl, 3)*zr4(jl)* &
               zrj0(jl, jaj, ikl)
       END DO
    END DO


    ! *         5.2   UPWARD FLUXES
    ! -------------


    DO jl = 1, kdlon
       zfu(jl, 1) = zfd(jl, 1)*palbp(jl, knu)
    END DO

    DO jk = 2, kflev + 1
       ikm1 = jk - 1
       DO jl = 1, kdlon
          zw3(jl, 1) = zw3(jl, 1) + pud(jl, 1, ikm1)*1.66
          zw3(jl, 2) = zw3(jl, 2) + pud(jl, 2, ikm1)*1.66
          zw3(jl, 3) = zw3(jl, 3) + poz(jl, ikm1)*1.66
          zw4(jl) = zw4(jl) + pud(jl, 4, ikm1)*1.66
          zw5(jl) = zw5(jl) + pud(jl, 5, ikm1)*1.66
       END DO

       CALL swtt1(knu, 3, iind3, zw3, zr3)

       DO jl = 1, kdlon
          ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
          zfu(jl, jk) = zr3(jl, 1)*zr3(jl, 2)*zr3(jl, 3)*zr4(jl)* &
               zrk0(jl, jaj, jk)
       END DO
    END DO


    ! ------------------------------------------------------------------

    ! *         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
    ! --------------------------------------------------

    iabs = 3

    ! *         6.1    DOWNWARD FLUXES
    ! ---------------

    DO jl = 1, kdlon
       zw1(jl) = 0.
       zw4(jl) = 0.
       zw5(jl) = 0.
       zr1(jl) = 0.
       pfdown(jl, kflev+1) = ((1.-pclear(jl))*pfdown(jl,kflev+1)+pclear(jl)*zfd( &
            jl,kflev+1))*rsun(knu)
    END DO

    DO jk = 1, kflev
       ikl = kflev + 1 - jk
       DO jl = 1, kdlon
          zw1(jl) = zw1(jl) + poz(jl, ikl)/zrmue(jl, ikl)
          zw4(jl) = zw4(jl) + pud(jl, 4, ikl)/zrmue(jl, ikl)
          zw5(jl) = zw5(jl) + pud(jl, 5, ikl)/zrmue(jl, ikl)
          ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
       END DO

       CALL swtt(knu, iabs, zw1, zr1)

       DO jl = 1, kdlon
          pfdown(jl, ikl) = ((1.-pclear(jl))*zr1(jl)*zr4(jl)*pfdown(jl,ikl)+ &
               pclear(jl)*zfd(jl,ikl))*rsun(knu)
       END DO
    END DO


    ! *         6.2    UPWARD FLUXES
    ! -------------

    DO jl = 1, kdlon
       pfup(jl, 1) = ((1.-pclear(jl))*zr1(jl)*zr4(jl)*pfup(jl,1)+pclear(jl)*zfu( &
            jl,1))*rsun(knu)
    END DO

    DO jk = 2, kflev + 1
       ikm1 = jk - 1
       DO jl = 1, kdlon
          zw1(jl) = zw1(jl) + poz(jl, ikm1)*1.66
          zw4(jl) = zw4(jl) + pud(jl, 4, ikm1)*1.66
          zw5(jl) = zw5(jl) + pud(jl, 5, ikm1)*1.66
          ! ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
       END DO

       CALL swtt(knu, iabs, zw1, zr1)

       DO jl = 1, kdlon
          pfup(jl, jk) = ((1.-pclear(jl))*zr1(jl)*zr4(jl)*pfup(jl,jk)+pclear(jl)* &
               zfu(jl,jk))*rsun(knu)
       END DO
    END DO

  END SUBROUTINE sw2s

end module sw2s_m
