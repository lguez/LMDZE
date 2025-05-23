module swclr_m

  IMPLICIT NONE

contains

  SUBROUTINE swclr(knu, albp, pdsig, prayl, psec, ppizaz, pray1, pray2, &
       prefz, prj, prk, prmu0, ptauaz, ptra1, ptra2)

    ! PURPOSE.
    ! COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
    ! CLEAR-SKY COLUMN

    ! REFERENCE.
    ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    ! AUTHOR.
    ! JEAN-JACQUES MORCRETTE *ECMWF*

    ! MODIFICATIONS.
    ! ORIGINAL : 94-11-15

    use dimensions, only: llm
    use dimphy, only: klon
    USE radepsi, only: zepsec
    USE radopt, only: novlp

    ! ARGUMENTS:

    INTEGER, intent(in):: knu
    DOUBLE PRECISION, intent(in):: albp(klon, 2)
    DOUBLE PRECISION, intent(in):: pdsig(klon, llm)
    DOUBLE PRECISION, intent(in):: prayl(klon)
    DOUBLE PRECISION psec(klon)

    DOUBLE PRECISION, intent(out):: ppizaz(klon, llm)
    DOUBLE PRECISION pray1(klon, llm + 1)
    DOUBLE PRECISION pray2(klon, llm + 1)
    DOUBLE PRECISION prefz(klon, 2, llm + 1)
    DOUBLE PRECISION prj(klon, 6, llm + 1)
    DOUBLE PRECISION prk(klon, 6, llm + 1)
    DOUBLE PRECISION prmu0(klon, llm + 1)
    DOUBLE PRECISION, intent(out):: ptauaz(klon, llm)
    DOUBLE PRECISION ptra1(klon, llm + 1)
    DOUBLE PRECISION ptra2(klon, llm + 1)

    ! LOCAL VARIABLES:
    DOUBLE PRECISION zc0i(klon, llm + 1)
    DOUBLE PRECISION zclear(klon)
    DOUBLE PRECISION zr21(klon)
    DOUBLE PRECISION zss0(klon)
    DOUBLE PRECISION zscat(klon)
    DOUBLE PRECISION ztr(klon, 2, llm + 1)
    INTEGER jl, jk, ja, jkl, jklp1, jaj, jkm1
    DOUBLE PRECISION zfacoa, zcorae
    DOUBLE PRECISION zmue, zgap, zww, zto, zden, zmu1, zden1
    DOUBLE PRECISION zbmu0, zbmu1, zre11
    double precision, parameter:: REPSCT = 1d-10

    !------------------------------------------------------------------
    
    ! 1. OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH

    DO jk = 1, llm + 1
       DO ja = 1, 6
          DO jl = 1, klon
             prj(jl, ja, jk) = 0d0
             prk(jl, ja, jk) = 0d0
          END DO
       END DO
    END DO

    DO jk = 1, llm
       DO jl = 1, klon
          ptauaz(jl, jk) = prayl(jl) * pdsig(jl, jk)
          ppizaz(jl, jk) = 1d0 - repsct
       END DO
    END DO

    ! 2. TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL

    DO jl = 1, klon
       zc0i(jl, llm + 1) = 0d0
       zclear(jl) = 1d0
       zscat(jl) = 0d0
    END DO

    jk = 1
    jkl = llm + 1 - jk
    jklp1 = jkl + 1
    DO jl = 1, klon
       zfacoa = 1d0
       zcorae = zfacoa * ptauaz(jl, jkl) * psec(jl)
       zr21(jl) = exp(- zcorae)
       zss0(jl) = 1d0 - zr21(jl)

       IF (novlp == 1) THEN
          ! maximum-random
          zclear(jl) = zclear(jl) * (1d0 - max(zss0(jl), zscat(jl))) / (1d0 &
               - min(zscat(jl), 1d0 - zepsec))
          zc0i(jl, jkl) = 1d0 - zclear(jl)
          zscat(jl) = zss0(jl)
       ELSE IF (novlp == 2) THEN
          ! maximum
          zscat(jl) = max(zss0(jl), zscat(jl))
          zc0i(jl, jkl) = zscat(jl)
       ELSE IF (novlp == 3) THEN
          ! random
          zclear(jl) = zclear(jl) * (1d0 - zss0(jl))
          zscat(jl) = 1d0 - zclear(jl)
          zc0i(jl, jkl) = zscat(jl)
       END IF
    END DO

    DO jk = 2, llm
       jkl = llm + 1 - jk
       jklp1 = jkl + 1
       DO jl = 1, klon
          zfacoa = 1d0
          zcorae = zfacoa * ptauaz(jl, jkl) * psec(jl)
          zr21(jl) = exp(- zcorae)
          zss0(jl) = 1d0 - zr21(jl)

          IF (novlp == 1) THEN
             ! maximum-random
             zclear(jl) = zclear(jl) * (1d0 - max(zss0(jl), zscat(jl))) &
                  / (1d0 - min(zscat(jl), 1d0 - zepsec))
             zc0i(jl, jkl) = 1d0 - zclear(jl)
             zscat(jl) = zss0(jl)
          ELSE IF (novlp == 2) THEN
             ! maximum
             zscat(jl) = max(zss0(jl), zscat(jl))
             zc0i(jl, jkl) = zscat(jl)
          ELSE IF (novlp == 3) THEN
             ! random
             zclear(jl) = zclear(jl) * (1d0 - zss0(jl))
             zscat(jl) = 1d0 - zclear(jl)
             zc0i(jl, jkl) = zscat(jl)
          END IF
       END DO
    END DO

    ! 3. REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING

    DO jl = 1, klon
       pray1(jl, llm + 1) = 0d0
       pray2(jl, llm + 1) = 0d0
       prefz(jl, 2, 1) = albp(jl, knu)
       prefz(jl, 1, 1) = albp(jl, knu)
       ptra1(jl, llm + 1) = 1d0
       ptra2(jl, llm + 1) = 1d0
    END DO

    DO jk = 2, llm + 1
       jkm1 = jk - 1
       DO jl = 1, klon

          ! 3.1 EQUIVALENT ZENITH ANGLE

          zmue = (1d0 - zc0i(jl, jk)) * psec(jl) + zc0i(jl, jk) * 1.66d0
          prmu0(jl, jk) = 1d0 / zmue

          ! 3.2 REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS

          zgap = 0d0
          zbmu0 = 0.5d0 - 0.75d0 * zgap / zmue
          zww = ppizaz(jl, jkm1)
          zto = ptauaz(jl, jkm1)
          zden = 1d0 + (1d0 - zww + zbmu0 * zww) * zto * zmue + (1d0 - zww) &
               * (1d0 - zww + 2d0 * zbmu0 * zww) * zto * zto * zmue * zmue
          pray1(jl, jkm1) = zbmu0 * zww * zto * zmue / zden
          ptra1(jl, jkm1) = 1d0 / zden

          zmu1 = 0.5d0
          zbmu1 = 0.5d0 - 0.75d0 * zgap * zmu1
          zden1 = 1d0 + (1d0 - zww + zbmu1 * zww) * zto / zmu1 + (1d0 - zww) &
               * (1d0 - zww + 2d0 * zbmu1 * zww &
               ) * zto * zto / zmu1 / zmu1
          pray2(jl, jkm1) = zbmu1 * zww * zto / zmu1 / zden1
          ptra2(jl, jkm1) = 1d0 / zden1

          prefz(jl, 1, jk) = (pray1(jl, jkm1) + prefz(jl, 1, jkm1) &
               * ptra1(jl, jkm1)* ptra2(jl, jkm1) / (1d0 - pray2(jl, jkm1) &
               * prefz(jl, 1, jkm1)))

          ztr(jl, 1, jkm1) = (ptra1(jl, jkm1) / (1d0 - pray2(jl, jkm1) &
               * prefz(jl, 1, jkm1)))

          prefz(jl, 2, jk) = (pray1(jl, jkm1) + prefz(jl, 2, jkm1) &
               * ptra1(jl, jkm1) * ptra2(jl, jkm1))

          ztr(jl, 2, jkm1) = ptra1(jl, jkm1)

       END DO
    END DO
    DO jl = 1, klon
       zmue = (1d0 - zc0i(jl, 1)) * psec(jl) + zc0i(jl, 1) * 1.66d0
       prmu0(jl, 1) = 1d0 / zmue
    END DO

    ! 3.5 REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL

    IF (knu == 1) THEN
       jaj = 2
       DO jl = 1, klon
          prj(jl, jaj, llm + 1) = 1d0
          prk(jl, jaj, llm + 1) = prefz(jl, 1, llm + 1)
       END DO

       DO jk = 1, llm
          jkl = llm + 1 - jk
          jklp1 = jkl + 1
          DO jl = 1, klon
             zre11 = prj(jl, jaj, jklp1) * ztr(jl, 1, jkl)
             prj(jl, jaj, jkl) = zre11
             prk(jl, jaj, jkl) = zre11 * prefz(jl, 1, jkl)
          END DO
       END DO
    ELSE
       DO jaj = 1, 2
          DO jl = 1, klon
             prj(jl, jaj, llm + 1) = 1d0
             prk(jl, jaj, llm + 1) = prefz(jl, jaj, llm + 1)
          END DO

          DO jk = 1, llm
             jkl = llm + 1 - jk
             jklp1 = jkl + 1
             DO jl = 1, klon
                zre11 = prj(jl, jaj, jklp1) * ztr(jl, jaj, jkl)
                prj(jl, jaj, jkl) = zre11
                prk(jl, jaj, jkl) = zre11 * prefz(jl, jaj, jkl)
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE swclr

end module swclr_m
