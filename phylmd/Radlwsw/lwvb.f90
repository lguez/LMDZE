SUBROUTINE lwvb(kuaer, ktraer, klim, pabcu, padjd, padju, pb, pbint, pbsui, &
    pbsur, pbtop, pdisd, pdisu, pemis, ppmb, pga, pgb, pgasur, pgbsur, &
    pgatop, pgbtop, pcts, pfluc)
  USE dimensions
  USE dimphy
  USE radopt
  USE raddimlw
  IMPLICIT NONE

  ! -----------------------------------------------------------------------
  ! PURPOSE.
  ! --------
  ! INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
  ! INTEGRATION

  ! METHOD.
  ! -------

  ! 1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
  ! ATMOSPHERE
  ! 2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
  ! TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
  ! 3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES

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
  ! Voigt lines (loop 2413 to 2427)  - JJM & PhD - 01/96
  ! -----------------------------------------------------------------------

  ! *       0.1   ARGUMENTS
  ! ---------

  INTEGER kuaer, ktraer, klim

  DOUBLE PRECISION pabcu(klon, nua, 3*llm+1) ! ABSORBER AMOUNTS
  DOUBLE PRECISION padjd(klon, llm+1) ! CONTRIBUTION BY ADJACENT LAYERS
  DOUBLE PRECISION padju(klon, llm+1) ! CONTRIBUTION BY ADJACENT LAYERS
  DOUBLE PRECISION pb(klon, ninter, llm+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
  DOUBLE PRECISION pbint(klon, llm+1) ! HALF-LEVEL PLANCK FUNCTIONS
  DOUBLE PRECISION pbsur(klon, ninter) ! SPECTRAL SURFACE PLANCK FUNCTION
  DOUBLE PRECISION pbsui(klon) ! SURFACE PLANCK FUNCTION
  DOUBLE PRECISION pbtop(klon, ninter) ! SPECTRAL T.O.A. PLANCK FUNCTION
  DOUBLE PRECISION pdisd(klon, llm+1) ! CONTRIBUTION BY DISTANT LAYERS
  DOUBLE PRECISION pdisu(klon, llm+1) ! CONTRIBUTION BY DISTANT LAYERS
  DOUBLE PRECISION pemis(klon) ! SURFACE EMISSIVITY
  DOUBLE PRECISION ppmb(klon, llm+1) ! PRESSURE MB
  DOUBLE PRECISION pga(klon, 8, 2, llm) ! PADE APPROXIMANTS
  DOUBLE PRECISION pgb(klon, 8, 2, llm) ! PADE APPROXIMANTS
  DOUBLE PRECISION pgasur(klon, 8, 2) ! SURFACE PADE APPROXIMANTS
  DOUBLE PRECISION pgbsur(klon, 8, 2) ! SURFACE PADE APPROXIMANTS
  DOUBLE PRECISION pgatop(klon, 8, 2) ! T.O.A. PADE APPROXIMANTS
  DOUBLE PRECISION pgbtop(klon, 8, 2) ! T.O.A. PADE APPROXIMANTS

  DOUBLE PRECISION pfluc(klon, 2, llm+1) ! CLEAR-SKY RADIATIVE FLUXES
  DOUBLE PRECISION pcts(klon, llm) ! COOLING-TO-SPACE TERM

  ! * LOCAL VARIABLES:

  DOUBLE PRECISION zbgnd(klon)
  DOUBLE PRECISION zfd(klon)
  DOUBLE PRECISION zfn10(klon)
  DOUBLE PRECISION zfu(klon)
  DOUBLE PRECISION ztt(klon, ntra)
  DOUBLE PRECISION ztt1(klon, ntra)
  DOUBLE PRECISION zuu(klon, nua)
  DOUBLE PRECISION zcnsol(klon)
  DOUBLE PRECISION zcntop(klon)

  INTEGER jk, jl, ja
  INTEGER jstra, jstru
  INTEGER in, jlim
  DOUBLE PRECISION zctstr
  ! -----------------------------------------------------------------------

  ! *         1.    INITIALIZATION
  ! --------------



  ! *         1.2     INITIALIZE TRANSMISSION FUNCTIONS
  ! ---------------------------------


  DO ja = 1, ntra
    DO jl = 1, klon
      ztt(jl, ja) = 1.0
      ztt1(jl, ja) = 1.0
    END DO
  END DO

  DO ja = 1, nua
    DO jl = 1, klon
      zuu(jl, ja) = 1.0
    END DO
  END DO

  ! ------------------------------------------------------------------

  ! *         2.      VERTICAL INTEGRATION
  ! --------------------

  ! *         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
  ! -----------------------------------


  DO jk = 1, llm
    in = (jk-1)*ng1p1 + 1

    DO ja = 1, kuaer
      DO jl = 1, klon
        zuu(jl, ja) = pabcu(jl, ja, in)
      END DO
    END DO


    CALL lwtt(pgatop(1,1,1), pgbtop(1,1,1), zuu, ztt)

    DO jl = 1, klon
      zcntop(jl) = pbtop(jl, 1)*ztt(jl, 1)*ztt(jl, 10) + &
        pbtop(jl, 2)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
        pbtop(jl, 3)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
        pbtop(jl, 4)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
        pbtop(jl, 5)*ztt(jl, 3)*ztt(jl, 14) + pbtop(jl, 6)*ztt(jl, 6)*ztt(jl, &
        15)
      zfd(jl) = zcntop(jl) - pbint(jl, jk) - pdisd(jl, jk) - padjd(jl, jk)
      pfluc(jl, 2, jk) = zfd(jl)
    END DO

  END DO

  jk = llm + 1
  in = (jk-1)*ng1p1 + 1

  DO jl = 1, klon
    zcntop(jl) = pbtop(jl, 1) + pbtop(jl, 2) + pbtop(jl, 3) + pbtop(jl, 4) + &
      pbtop(jl, 5) + pbtop(jl, 6)
    zfd(jl) = zcntop(jl) - pbint(jl, jk) - pdisd(jl, jk) - padjd(jl, jk)
    pfluc(jl, 2, jk) = zfd(jl)
  END DO

  ! *         2.4     COOLING-TO-SPACE OF LAYERS ABOVE 10 HPA
  ! ---------------------------------------



  ! *         2.4.1   INITIALIZATION
  ! --------------


  jlim = llm

  IF (.NOT. levoigt) THEN
    DO jk = llm, 1, -1
      IF (ppmb(1,jk)<10.0) THEN
        jlim = jk
      END IF
    END DO
  END IF
  klim = jlim

  IF (.NOT. levoigt) THEN
    DO ja = 1, ktraer
      DO jl = 1, klon
        ztt1(jl, ja) = 1.0
      END DO
    END DO

    ! *         2.4.2   LOOP OVER LAYERS ABOVE 10 HPA
    ! -----------------------------


    DO jstra = llm, jlim, -1
      jstru = (jstra-1)*ng1p1 + 1

      DO ja = 1, kuaer
        DO jl = 1, klon
          zuu(jl, ja) = pabcu(jl, ja, jstru)
        END DO
      END DO


      CALL lwtt(pga(1,1,1,jstra), pgb(1,1,1,jstra), zuu, ztt)

      DO jl = 1, klon
        zctstr = (pb(jl,1,jstra)+pb(jl,1,jstra+1))* &
          (ztt1(jl,1)*ztt1(jl,10)-ztt(jl,1)*ztt(jl,10)) + &
          (pb(jl,2,jstra)+pb(jl,2,jstra+1))*(ztt1(jl,2)*ztt1(jl,7)*ztt1(jl,11 &
          )-ztt(jl,2)*ztt(jl,7)*ztt(jl,11)) + (pb(jl,3,jstra)+pb(jl,3,jstra+1 &
          ))*(ztt1(jl,4)*ztt1(jl,8)*ztt1(jl,12)-ztt(jl,4)*ztt(jl,8)*ztt(jl,12 &
          )) + (pb(jl,4,jstra)+pb(jl,4,jstra+1))*(ztt1(jl,5)*ztt1(jl,9)*ztt1( &
          jl,13)-ztt(jl,5)*ztt(jl,9)*ztt(jl,13)) + (pb(jl,5,jstra)+pb(jl,5, &
          jstra+1))*(ztt1(jl,3)*ztt1(jl,14)-ztt(jl,3)*ztt(jl,14)) + &
          (pb(jl,6,jstra)+pb(jl,6,jstra+1))*(ztt1(jl,6)*ztt1(jl,15)-ztt(jl,6) &
          *ztt(jl,15))
        pcts(jl, jstra) = zctstr*0.5
      END DO
      DO ja = 1, ktraer
        DO jl = 1, klon
          ztt1(jl, ja) = ztt(jl, ja)
        END DO
      END DO
    END DO
  END IF
  ! Mise a zero de securite pour PCTS en cas de LEVOIGT
  IF (levoigt) THEN
    DO jstra = 1, llm
      DO jl = 1, klon
        pcts(jl, jstra) = 0.
      END DO
    END DO
  END IF


  ! *         2.5     EXCHANGE WITH LOWER LIMIT
  ! -------------------------


  DO jl = 1, klon
    zbgnd(jl) = pbsui(jl)*pemis(jl) - (1.-pemis(jl))*pfluc(jl, 2, 1) - &
      pbint(jl, 1)
  END DO

  jk = 1
  in = (jk-1)*ng1p1 + 1

  DO jl = 1, klon
    zcnsol(jl) = pbsur(jl, 1) + pbsur(jl, 2) + pbsur(jl, 3) + pbsur(jl, 4) + &
      pbsur(jl, 5) + pbsur(jl, 6)
    zcnsol(jl) = zcnsol(jl)*zbgnd(jl)/pbsui(jl)
    zfu(jl) = zcnsol(jl) + pbint(jl, jk) - pdisu(jl, jk) - padju(jl, jk)
    pfluc(jl, 1, jk) = zfu(jl)
  END DO

  DO jk = 2, llm + 1
    in = (jk-1)*ng1p1 + 1


    DO ja = 1, kuaer
      DO jl = 1, klon
        zuu(jl, ja) = pabcu(jl, ja, 1) - pabcu(jl, ja, in)
      END DO
    END DO


    CALL lwtt(pgasur(1,1,1), pgbsur(1,1,1), zuu, ztt)

    DO jl = 1, klon
      zcnsol(jl) = pbsur(jl, 1)*ztt(jl, 1)*ztt(jl, 10) + &
        pbsur(jl, 2)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
        pbsur(jl, 3)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
        pbsur(jl, 4)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
        pbsur(jl, 5)*ztt(jl, 3)*ztt(jl, 14) + pbsur(jl, 6)*ztt(jl, 6)*ztt(jl, &
        15)
      zcnsol(jl) = zcnsol(jl)*zbgnd(jl)/pbsui(jl)
      zfu(jl) = zcnsol(jl) + pbint(jl, jk) - pdisu(jl, jk) - padju(jl, jk)
      pfluc(jl, 1, jk) = zfu(jl)
    END DO


  END DO



  ! *         2.7     CLEAR-SKY FLUXES
  ! ----------------


  IF (.NOT. levoigt) THEN
    DO jl = 1, klon
      zfn10(jl) = pfluc(jl, 1, jlim) + pfluc(jl, 2, jlim)
    END DO
    DO jk = jlim + 1, llm + 1
      DO jl = 1, klon
        zfn10(jl) = zfn10(jl) + pcts(jl, jk-1)
        pfluc(jl, 1, jk) = zfn10(jl)
        pfluc(jl, 2, jk) = 0.
      END DO
    END DO
  END IF

  ! ------------------------------------------------------------------

  RETURN
END SUBROUTINE lwvb
