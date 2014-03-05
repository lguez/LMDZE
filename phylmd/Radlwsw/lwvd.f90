SUBROUTINE lwvd(kuaer, ktraer, pabcu, pdbdt, pga, pgb, pcntrb, pdisd, pdisu)
  USE dimens_m
  USE dimphy
  USE raddim
  USE raddimlw
  IMPLICIT NONE

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS

  ! METHOD.
  ! -------

  ! 1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
  ! CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE

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
  ! -----------------------------------------------------------------------
  ! * ARGUMENTS:

  INTEGER kuaer, ktraer

  DOUBLE PRECISION pabcu(kdlon, nua, 3*kflev+1) ! ABSORBER AMOUNTS
  DOUBLE PRECISION pdbdt(kdlon, ninter, kflev) ! LAYER PLANCK FUNCTION GRADIENT
  DOUBLE PRECISION pga(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
  DOUBLE PRECISION pgb(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS

  DOUBLE PRECISION pcntrb(kdlon, kflev+1, kflev+1) ! ENERGY EXCHANGE MATRIX
  DOUBLE PRECISION pdisd(kdlon, kflev+1) !  CONTRIBUTION BY DISTANT LAYERS
  DOUBLE PRECISION pdisu(kdlon, kflev+1) !  CONTRIBUTION BY DISTANT LAYERS

  ! * LOCAL VARIABLES:

  DOUBLE PRECISION zglayd(kdlon)
  DOUBLE PRECISION zglayu(kdlon)
  DOUBLE PRECISION ztt(kdlon, ntra)
  DOUBLE PRECISION ztt1(kdlon, ntra)
  DOUBLE PRECISION ztt2(kdlon, ntra)

  INTEGER jl, jk, ja, ikp1, ikn, ikd1, jkj, ikd2
  INTEGER ikjp1, ikm1, ikj, jlk, iku1, ijkl, iku2
  INTEGER ind1, ind2, ind3, ind4, itt
  DOUBLE PRECISION zww, zdzxdg, zdzxmg

  ! *         1.    INITIALIZATION
  ! --------------


  ! *         1.1     INITIALIZE LAYER CONTRIBUTIONS
  ! ------------------------------


  DO jk = 1, kflev + 1
    DO jl = 1, kdlon
      pdisd(jl, jk) = 0.
      pdisu(jl, jk) = 0.
    END DO
  END DO

  ! *         1.2     INITIALIZE TRANSMISSION FUNCTIONS
  ! ---------------------------------



  DO ja = 1, ntra
    DO jl = 1, kdlon
      ztt(jl, ja) = 1.0
      ztt1(jl, ja) = 1.0
      ztt2(jl, ja) = 1.0
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      VERTICAL INTEGRATION
  ! --------------------


  ind1 = 0
  ind3 = 0
  ind4 = 1
  ind2 = 1


  ! *         2.2     CONTRIBUTION FROM DISTANT LAYERS
  ! ---------------------------------



  ! *         2.2.1   DISTANT AND ABOVE LAYERS
  ! ------------------------




  ! *         2.2.2   FIRST UPPER LEVEL
  ! -----------------


  DO jk = 1, kflev - 1
    ikp1 = jk + 1
    ikn = (jk-1)*ng1p1 + 1
    ikd1 = jk*ng1p1 + 1

    CALL lwttm(pga(1,1,1,jk), pgb(1,1,1,jk), pabcu(1,1,ikn), pabcu(1,1,ikd1), &
      ztt1)



    ! *         2.2.3   HIGHER UP
    ! ---------


    itt = 1
    DO jkj = ikp1, kflev
      IF (itt==1) THEN
        itt = 2
      ELSE
        itt = 1
      END IF
      ikjp1 = jkj + 1
      ikd2 = jkj*ng1p1 + 1

      IF (itt==1) THEN
        CALL lwttm(pga(1,1,1,jkj), pgb(1,1,1,jkj), pabcu(1,1,ikn), &
          pabcu(1,1,ikd2), ztt1)
      ELSE
        CALL lwttm(pga(1,1,1,jkj), pgb(1,1,1,jkj), pabcu(1,1,ikn), &
          pabcu(1,1,ikd2), ztt2)
      END IF

      DO ja = 1, ktraer
        DO jl = 1, kdlon
          ztt(jl, ja) = (ztt1(jl,ja)+ztt2(jl,ja))*0.5
        END DO
      END DO

      DO jl = 1, kdlon
        zww = pdbdt(jl, 1, jkj)*ztt(jl, 1)*ztt(jl, 10) + &
          pdbdt(jl, 2, jkj)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
          pdbdt(jl, 3, jkj)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
          pdbdt(jl, 4, jkj)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
          pdbdt(jl, 5, jkj)*ztt(jl, 3)*ztt(jl, 14) + &
          pdbdt(jl, 6, jkj)*ztt(jl, 6)*ztt(jl, 15)
        zglayd(jl) = zww
        zdzxdg = zglayd(jl)
        pdisd(jl, jk) = pdisd(jl, jk) + zdzxdg
        pcntrb(jl, jk, ikjp1) = zdzxdg
      END DO


    END DO
  END DO


  ! *         2.2.4   DISTANT AND BELOW LAYERS
  ! ------------------------




  ! *         2.2.5   FIRST LOWER LEVEL
  ! -----------------


  DO jk = 3, kflev + 1
    ikn = (jk-1)*ng1p1 + 1
    ikm1 = jk - 1
    ikj = jk - 2
    iku1 = ikj*ng1p1 + 1


    CALL lwttm(pga(1,1,1,ikj), pgb(1,1,1,ikj), pabcu(1,1,iku1), &
      pabcu(1,1,ikn), ztt1)



    ! *         2.2.6   DOWN BELOW
    ! ----------


    itt = 1
    DO jlk = 1, ikj
      IF (itt==1) THEN
        itt = 2
      ELSE
        itt = 1
      END IF
      ijkl = ikm1 - jlk
      iku2 = (ijkl-1)*ng1p1 + 1


      IF (itt==1) THEN
        CALL lwttm(pga(1,1,1,ijkl), pgb(1,1,1,ijkl), pabcu(1,1,iku2), &
          pabcu(1,1,ikn), ztt1)
      ELSE
        CALL lwttm(pga(1,1,1,ijkl), pgb(1,1,1,ijkl), pabcu(1,1,iku2), &
          pabcu(1,1,ikn), ztt2)
      END IF

      DO ja = 1, ktraer
        DO jl = 1, kdlon
          ztt(jl, ja) = (ztt1(jl,ja)+ztt2(jl,ja))*0.5
        END DO
      END DO

      DO jl = 1, kdlon
        zww = pdbdt(jl, 1, ijkl)*ztt(jl, 1)*ztt(jl, 10) + &
          pdbdt(jl, 2, ijkl)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
          pdbdt(jl, 3, ijkl)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
          pdbdt(jl, 4, ijkl)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
          pdbdt(jl, 5, ijkl)*ztt(jl, 3)*ztt(jl, 14) + &
          pdbdt(jl, 6, ijkl)*ztt(jl, 6)*ztt(jl, 15)
        zglayu(jl) = zww
        zdzxmg = zglayu(jl)
        pdisu(jl, jk) = pdisu(jl, jk) + zdzxmg
        pcntrb(jl, jk, ijkl) = zdzxmg
      END DO


    END DO
  END DO

  RETURN
END SUBROUTINE lwvd
