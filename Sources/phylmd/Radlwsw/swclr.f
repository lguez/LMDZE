SUBROUTINE swclr(knu, paer, flag_aer, tauae, pizae, cgae, palbp, pdsig, &
    prayl, psec, pcgaz, ppizaz, pray1, pray2, prefz, prj, prk, prmu0, ptauaz, &
    ptra1, ptra2)
  USE dimens_m
  USE dimphy
  USE raddim
  USE radepsi
  USE radopt
  IMPLICIT NONE

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
  DOUBLE PRECISION flag_aer
  DOUBLE PRECISION tauae(kdlon, kflev, 2)
  DOUBLE PRECISION pizae(kdlon, kflev, 2)
  DOUBLE PRECISION cgae(kdlon, kflev, 2)
  DOUBLE PRECISION paer(kdlon, kflev, 5)
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
  DOUBLE PRECISION zcle0(kdlon, kflev)
  DOUBLE PRECISION zclear(kdlon)
  DOUBLE PRECISION zr21(kdlon)
  DOUBLE PRECISION zr23(kdlon)
  DOUBLE PRECISION zss0(kdlon)
  DOUBLE PRECISION zscat(kdlon)
  DOUBLE PRECISION ztr(kdlon, 2, kflev+1)

  INTEGER jl, jk, ja, jkl, jklp1, jaj, jkm1, in
  DOUBLE PRECISION ztray, zgar, zratio, zff, zfacoa, zcorae
  DOUBLE PRECISION zmue, zgap, zww, zto, zden, zmu1, zden1
  DOUBLE PRECISION zbmu0, zbmu1, zre11

  ! * Prescribed Data for Aerosols:

  DOUBLE PRECISION taua(2, 5), rpiza(2, 5), rcga(2, 5)
  SAVE taua, rpiza, rcga
  DATA ((taua(in,ja),ja=1,5), in=1, 2)/.730719, .912819, .725059, .745405, &
    .682188, .730719, .912819, .725059, .745405, .682188/
  DATA ((rpiza(in,ja),ja=1,5), in=1, 2)/.872212, .982545, .623143, .944887, &
    .997975, .872212, .982545, .623143, .944887, .997975/
  DATA ((rcga(in,ja),ja=1,5), in=1, 2)/.647596, .739002, .580845, .662657, &
    .624246, .647596, .739002, .580845, .662657, .624246/
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
    ! -OB
    ! DO 104 JL = 1, KDLON
    ! PCGAZ(JL,JK) = 0.
    ! PPIZAZ(JL,JK) =  0.
    ! PTAUAZ(JL,JK) = 0.
    ! 104  CONTINUE
    ! -OB
    ! DO 106 JAE=1,5
    ! DO 105 JL = 1, KDLON
    ! PTAUAZ(JL,JK)=PTAUAZ(JL,JK)
    ! S        +PAER(JL,JK,JAE)*TAUA(KNU,JAE)
    ! PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JK,JAE)
    ! S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)
    ! PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JK,JAE)
    ! S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
    ! 105  CONTINUE
    ! 106  CONTINUE
    ! -OB
    DO jl = 1, kdlon
      ptauaz(jl, jk) = flag_aer*tauae(jl, jk, knu)
      ppizaz(jl, jk) = flag_aer*pizae(jl, jk, knu)
      pcgaz(jl, jk) = flag_aer*cgae(jl, jk, knu)
    END DO

    IF (flag_aer>0) THEN
      ! -OB
      DO jl = 1, kdlon
        ! PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
        ! PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
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
    END IF ! check flag_aer
    ! 107  CONTINUE
    ! PRINT 9107,JK,((PAER(JL,JK,JAE),JAE=1,5)
    ! $ ,PTAUAZ(JL,JK),PPIZAZ(JL,JK),PCGAZ(JL,JK),JL=1,KDLON)
    ! 9107 FORMAT(1X,'SWCLR_107',I3,8E12.5)

  END DO

  ! ------------------------------------------------------------------

  ! *         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
  ! ----------------------------------------------


  DO jl = 1, kdlon
    zr23(jl) = 0.
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
    zcle0(jl, jkl) = zss0(jl)

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
      zcle0(jl, jkl) = zss0(jl)

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

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE swclr
