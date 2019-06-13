SUBROUTINE lwc(klim, pcldld, pcldlu, pemis, pfluc, pbint, pbsuin, pcts, &
    pcntrb, pflux)
  USE dimensions
  USE dimphy
  USE raddim
  USE radepsi
  USE radopt
  IMPLICIT NONE

  ! PURPOSE.
  ! --------
  ! INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
  ! RADIANCES

  ! EXPLICIT ARGUMENTS :
  ! --------------------
  ! ==== INPUTS ===
  ! PBINT  : (KDLON,0:KFLEV)     ; HALF LEVEL PLANCK FUNCTION
  ! PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
  ! PCLDLD : (KDLON,KFLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
  ! PCLDLU : (KDLON,KFLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
  ! PCNTRB : (KDLON,KFLEV+1,KFLEV+1); CLEAR-SKY ENERGY EXCHANGE
  ! PCTS   : (KDLON,KFLEV)       ; CLEAR-SKY LAYER COOLING-TO-SPACE
  ! PEMIS  : (KDLON)             ; SURFACE EMISSIVITY
  ! PFLUC
  ! ==== OUTPUTS ===
  ! PFLUX(KDLON,2,KFLEV)         ; RADIATIVE FLUXES :
  ! 1  ==>  UPWARD   FLUX TOTAL
  ! 2  ==>  DOWNWARD FLUX TOTAL

  ! METHOD.
  ! -------

  ! 1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
  ! 2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
  ! 3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
  ! CLOUDS

  ! REFERENCE.
  ! ----------

  ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
  ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  ! AUTHOR.
  ! -------
  ! JEAN-JACQUES MORCRETTE  *ECMWF*

  ! MODIFICATIONS.
  ! --------------
  ! ORIGINAL : 89-07-14
  ! Voigt lines (loop 231 to 233)  - JJM & PhD - 01/96
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:
  INTEGER klim
  DOUBLE PRECISION pfluc(kdlon, 2, kflev+1) ! CLEAR-SKY RADIATIVE FLUXES
  DOUBLE PRECISION pbint(kdlon, kflev+1) ! HALF LEVEL PLANCK FUNCTION
  DOUBLE PRECISION pbsuin(kdlon) ! SURFACE PLANCK FUNCTION
  DOUBLE PRECISION pcntrb(kdlon, kflev+1, kflev+1) !CLEAR-SKY ENERGY EXCHANGE
  DOUBLE PRECISION pcts(kdlon, kflev) ! CLEAR-SKY LAYER COOLING-TO-SPACE

  DOUBLE PRECISION pcldld(kdlon, kflev)
  DOUBLE PRECISION pcldlu(kdlon, kflev)
  DOUBLE PRECISION pemis(kdlon)

  DOUBLE PRECISION pflux(kdlon, 2, kflev+1)
  ! -----------------------------------------------------------------------
  ! * LOCAL VARIABLES:
  INTEGER imx(kdlon), imxp(kdlon)

  DOUBLE PRECISION zclear(kdlon), zcloud(kdlon)
  DOUBLE PRECISION zdnf(kdlon, kflev+1, kflev+1), zfd(kdlon), zfn10(kdlon), &
    zfu(kdlon), zupf(kdlon, kflev+1, kflev+1)
  DOUBLE PRECISION zclm(kdlon, kflev+1, kflev+1)

  INTEGER jk, jl, imaxc, imx1, imx2, jkj, jkp1, jkm1
  INTEGER jk1, jk2, jkc, jkcp1, jcloud
  DOUBLE PRECISION zcfrac
  ! ------------------------------------------------------------------

  ! *         1.     INITIALIZATION
  ! --------------


  imaxc = 0

  DO jl = 1, kdlon
    imx(jl) = 0
    imxp(jl) = 0
    zcloud(jl) = 0.
  END DO

  ! *         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
  ! -------------------------------------------


  DO jk = 1, kflev
    DO jl = 1, kdlon
      imx1 = imx(jl)
      imx2 = jk
      IF (pcldlu(jl,jk)>zepsc) THEN
        imxp(jl) = imx2
      ELSE
        imxp(jl) = imx1
      END IF
      imaxc = max(imxp(jl), imaxc)
      imx(jl) = imxp(jl)
    END DO
  END DO
  ! GM*******
  imaxc = kflev
  ! GM*******

  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pflux(jl, 1, jk) = pfluc(jl, 1, jk)
      pflux(jl, 2, jk) = pfluc(jl, 2, jk)
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
  ! ---------------------------------------

  IF (imaxc>0) THEN
    ! *         2.0     INITIALIZE TO CLEAR-SKY FLUXES
    ! ------------------------------


    DO jk1 = 1, kflev + 1
      DO jk2 = 1, kflev + 1
        DO jl = 1, kdlon
          zupf(jl, jk2, jk1) = pfluc(jl, 1, jk1)
          zdnf(jl, jk2, jk1) = pfluc(jl, 2, jk1)
        END DO
      END DO
    END DO

    ! *         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
    ! ----------------------------------------------


    DO jkc = 1, imaxc
      jcloud = jkc
      jkcp1 = jcloud + 1

      ! *         2.1.1   ABOVE THE CLOUD
      ! ---------------


      DO jk = jkcp1, kflev + 1
        jkm1 = jk - 1
        DO jl = 1, kdlon
          zfu(jl) = 0.
        END DO
        IF (jk>jkcp1) THEN
          DO jkj = jkcp1, jkm1
            DO jl = 1, kdlon
              zfu(jl) = zfu(jl) + pcntrb(jl, jk, jkj)
            END DO
          END DO
        END IF

        DO jl = 1, kdlon
          zupf(jl, jkcp1, jk) = pbint(jl, jk) - zfu(jl)
        END DO
      END DO

      ! *         2.1.2   BELOW THE CLOUD
      ! ---------------


      DO jk = 1, jcloud
        jkp1 = jk + 1
        DO jl = 1, kdlon
          zfd(jl) = 0.
        END DO

        IF (jk<jcloud) THEN
          DO jkj = jkp1, jcloud
            DO jl = 1, kdlon
              zfd(jl) = zfd(jl) + pcntrb(jl, jk, jkj)
            END DO
          END DO
        END IF
        DO jl = 1, kdlon
          zdnf(jl, jkcp1, jk) = -pbint(jl, jk) - zfd(jl)
        END DO
      END DO

    END DO


    ! *         2.2     CLOUD COVER MATRIX
    ! ------------------

    ! *    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
    ! HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1


    DO jk1 = 1, kflev + 1
      DO jk2 = 1, kflev + 1
        DO jl = 1, kdlon
          zclm(jl, jk1, jk2) = 0.
        END DO
      END DO
    END DO



    ! *         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
    ! ------------------------------------------


    DO jk1 = 2, kflev + 1
      DO jl = 1, kdlon
        zclear(jl) = 1.
        zcloud(jl) = 0.
      END DO
      DO jk = jk1 - 1, 1, -1
        DO jl = 1, kdlon
          IF (novlp==1) THEN
            ! * maximum-random
            zclear(jl) = zclear(jl)*(1.0-max(pcldlu(jl, &
              jk),zcloud(jl)))/(1.0-min(zcloud(jl),1.-zepsec))
            zclm(jl, jk1, jk) = 1.0 - zclear(jl)
            zcloud(jl) = pcldlu(jl, jk)
          ELSE IF (novlp==2) THEN
            ! * maximum
            zcloud(jl) = max(zcloud(jl), pcldlu(jl,jk))
            zclm(jl, jk1, jk) = zcloud(jl)
          ELSE IF (novlp==3) THEN
            ! * random
            zclear(jl) = zclear(jl)*(1.0-pcldlu(jl,jk))
            zcloud(jl) = 1.0 - zclear(jl)
            zclm(jl, jk1, jk) = zcloud(jl)
          END IF
        END DO
      END DO
    END DO


    ! *         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
    ! ------------------------------------------


    DO jk1 = 1, kflev
      DO jl = 1, kdlon
        zclear(jl) = 1.
        zcloud(jl) = 0.
      END DO
      DO jk = jk1, kflev
        DO jl = 1, kdlon
          IF (novlp==1) THEN
            ! * maximum-random
            zclear(jl) = zclear(jl)*(1.0-max(pcldld(jl, &
              jk),zcloud(jl)))/(1.0-min(zcloud(jl),1.-zepsec))
            zclm(jl, jk1, jk) = 1.0 - zclear(jl)
            zcloud(jl) = pcldld(jl, jk)
          ELSE IF (novlp==2) THEN
            ! * maximum
            zcloud(jl) = max(zcloud(jl), pcldld(jl,jk))
            zclm(jl, jk1, jk) = zcloud(jl)
          ELSE IF (novlp==3) THEN
            ! * random
            zclear(jl) = zclear(jl)*(1.0-pcldld(jl,jk))
            zcloud(jl) = 1.0 - zclear(jl)
            zclm(jl, jk1, jk) = zcloud(jl)
          END IF
        END DO
      END DO
    END DO



    ! *         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
    ! ----------------------------------------------


    ! *         3.1     DOWNWARD FLUXES
    ! ---------------


    DO jl = 1, kdlon
      pflux(jl, 2, kflev+1) = 0.
    END DO

    DO jk1 = kflev, 1, -1

      ! *                 CONTRIBUTION FROM CLEAR-SKY FRACTION

      DO jl = 1, kdlon
        zfd(jl) = (1.-zclm(jl,jk1,kflev))*zdnf(jl, 1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM ADJACENT CLOUD

      DO jl = 1, kdlon
        zfd(jl) = zfd(jl) + zclm(jl, jk1, jk1)*zdnf(jl, jk1+1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

      DO jk = kflev - 1, jk1, -1
        DO jl = 1, kdlon
          zcfrac = zclm(jl, jk1, jk+1) - zclm(jl, jk1, jk)
          zfd(jl) = zfd(jl) + zcfrac*zdnf(jl, jk+2, jk1)
        END DO
      END DO

      DO jl = 1, kdlon
        pflux(jl, 2, jk1) = zfd(jl)
      END DO

    END DO




    ! *         3.2     UPWARD FLUX AT THE SURFACE
    ! --------------------------


    DO jl = 1, kdlon
      pflux(jl, 1, 1) = pemis(jl)*pbsuin(jl) - (1.-pemis(jl))*pflux(jl, 2, 1)
    END DO



    ! *         3.3     UPWARD FLUXES
    ! -------------


    DO jk1 = 2, kflev + 1

      ! *                 CONTRIBUTION FROM CLEAR-SKY FRACTION

      DO jl = 1, kdlon
        zfu(jl) = (1.-zclm(jl,jk1,1))*zupf(jl, 1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM ADJACENT CLOUD

      DO jl = 1, kdlon
        zfu(jl) = zfu(jl) + zclm(jl, jk1, jk1-1)*zupf(jl, jk1, jk1)
      END DO

      ! *                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

      DO jk = 2, jk1 - 1
        DO jl = 1, kdlon
          zcfrac = zclm(jl, jk1, jk-1) - zclm(jl, jk1, jk)
          zfu(jl) = zfu(jl) + zcfrac*zupf(jl, jk, jk1)
        END DO
      END DO

      DO jl = 1, kdlon
        pflux(jl, 1, jk1) = zfu(jl)
      END DO

    END DO


  END IF


  ! *         2.3     END OF CLOUD EFFECT COMPUTATIONS


  IF (.NOT. levoigt) THEN
    DO jl = 1, kdlon
      zfn10(jl) = pflux(jl, 1, klim) + pflux(jl, 2, klim)
    END DO
    DO jk = klim + 1, kflev + 1
      DO jl = 1, kdlon
        zfn10(jl) = zfn10(jl) + pcts(jl, jk-1)
        pflux(jl, 1, jk) = zfn10(jl)
        pflux(jl, 2, jk) = 0.0
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE lwc
