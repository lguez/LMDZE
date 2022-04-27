module lwv_m

  IMPLICIT NONE

contains

  SUBROUTINE lwv(kuaer, ktraer, klim, pabcu, pb, pbint, pbsuin, pbsur, pbtop, &
       pdbsl, pemis, ppmb, pga, pgb, pgasur, pgbsur, pgatop, pgbtop, &
       pcntrb, pcts, pfluc)
    USE dimensions
    USE dimphy
    use lwvd_m, only: lwvd
    use lwvn_m, only: lwvn
    USE suphec_m
    USE raddimlw

    ! -----------------------------------------------------------------------
    ! PURPOSE.
    ! --------
    ! CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
    ! FLUXES OR RADIANCES

    ! METHOD.
    ! -------

    ! 1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
    ! CONTRIBUTIONS BY -  THE NEARBY LAYERS
    ! -  THE DISTANT LAYERS
    ! -  THE BOUNDARY TERMS
    ! 2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.

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
    INTEGER kuaer, ktraer, klim

    DOUBLE PRECISION pabcu(klon, nua, 3*llm+1) ! EFFECTIVE ABSORBER AMOUNTS
    DOUBLE PRECISION pb(klon, ninter, llm+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
    DOUBLE PRECISION pbint(klon, llm+1) ! HALF-LEVEL PLANCK FUNCTIONS
    DOUBLE PRECISION pbsur(klon, ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION pbsuin(klon) ! SURFACE PLANCK FUNCTION
    DOUBLE PRECISION pbtop(klon, ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION pdbsl(klon, ninter, llm*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
    DOUBLE PRECISION pemis(klon) ! SURFACE EMISSIVITY
    DOUBLE PRECISION ppmb(klon, llm+1) ! HALF-LEVEL PRESSURE (MB)
    DOUBLE PRECISION pga(klon, 8, 2, llm) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgb(klon, 8, 2, llm) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgasur(klon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgbsur(klon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgatop(klon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgbtop(klon, 8, 2) ! PADE APPROXIMANTS

    DOUBLE PRECISION pcntrb(klon, llm+1, llm+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
    DOUBLE PRECISION pcts(klon, llm) ! COOLING-TO-SPACE TERM
    DOUBLE PRECISION pfluc(klon, 2, llm+1) ! CLEAR-SKY RADIATIVE FLUXES
    ! -----------------------------------------------------------------------
    ! LOCAL VARIABLES:
    DOUBLE PRECISION zadjd(klon, llm+1)
    DOUBLE PRECISION zadju(klon, llm+1)
    DOUBLE PRECISION zdbdt(klon, ninter, llm)
    DOUBLE PRECISION zdisd(klon, llm+1)
    DOUBLE PRECISION zdisu(klon, llm+1)

    INTEGER jk, jl
    ! -----------------------------------------------------------------------

    DO jk = 1, llm + 1
       DO jl = 1, klon
          zadjd(jl, jk) = 0.
          zadju(jl, jk) = 0.
          zdisd(jl, jk) = 0.
          zdisu(jl, jk) = 0.
       END DO
    END DO

    DO jk = 1, llm
       DO jl = 1, klon
          pcts(jl, jk) = 0.
       END DO
    END DO

    ! * CONTRIBUTION FROM ADJACENT LAYERS

    CALL lwvn(kuaer, pabcu, pdbsl, pga, pgb, zadjd, zadju, pcntrb, zdbdt)
    ! * CONTRIBUTION FROM DISTANT LAYERS

    CALL lwvd(ktraer, pabcu, zdbdt, pga, pgb, pcntrb, zdisd, zdisu)

    ! * EXCHANGE WITH THE BOUNDARIES

    CALL lwvb(kuaer, ktraer, klim, pabcu, zadjd, zadju, pb, pbint, pbsuin, &
         pbsur, pbtop, zdisd, zdisu, pemis, ppmb, pga, pgb, pgasur, pgbsur, &
         pgatop, pgbtop, pcts, pfluc)


    RETURN
  END SUBROUTINE lwv

end module lwv_m
