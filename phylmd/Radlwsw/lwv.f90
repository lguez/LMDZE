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
    use conf_phys_m, only: kdlon
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

    DOUBLE PRECISION pabcu(kdlon, nua, 3*llm+1) ! EFFECTIVE ABSORBER AMOUNTS
    DOUBLE PRECISION pb(kdlon, ninter, llm+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
    DOUBLE PRECISION pbint(kdlon, llm+1) ! HALF-LEVEL PLANCK FUNCTIONS
    DOUBLE PRECISION pbsur(kdlon, ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION pbsuin(kdlon) ! SURFACE PLANCK FUNCTION
    DOUBLE PRECISION pbtop(kdlon, ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION pdbsl(kdlon, ninter, llm*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
    DOUBLE PRECISION pemis(kdlon) ! SURFACE EMISSIVITY
    DOUBLE PRECISION ppmb(kdlon, llm+1) ! HALF-LEVEL PRESSURE (MB)
    DOUBLE PRECISION pga(kdlon, 8, 2, llm) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgb(kdlon, 8, 2, llm) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgasur(kdlon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgbsur(kdlon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgatop(kdlon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgbtop(kdlon, 8, 2) ! PADE APPROXIMANTS

    DOUBLE PRECISION pcntrb(kdlon, llm+1, llm+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
    DOUBLE PRECISION pcts(kdlon, llm) ! COOLING-TO-SPACE TERM
    DOUBLE PRECISION pfluc(kdlon, 2, llm+1) ! CLEAR-SKY RADIATIVE FLUXES
    ! -----------------------------------------------------------------------
    ! LOCAL VARIABLES:
    DOUBLE PRECISION zadjd(kdlon, llm+1)
    DOUBLE PRECISION zadju(kdlon, llm+1)
    DOUBLE PRECISION zdbdt(kdlon, ninter, llm)
    DOUBLE PRECISION zdisd(kdlon, llm+1)
    DOUBLE PRECISION zdisu(kdlon, llm+1)

    INTEGER jk, jl
    ! -----------------------------------------------------------------------

    DO jk = 1, llm + 1
       DO jl = 1, kdlon
          zadjd(jl, jk) = 0.
          zadju(jl, jk) = 0.
          zdisd(jl, jk) = 0.
          zdisu(jl, jk) = 0.
       END DO
    END DO

    DO jk = 1, llm
       DO jl = 1, kdlon
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
