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
    USE raddim
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

    DOUBLE PRECISION pabcu(kdlon, nua, 3*kflev+1) ! EFFECTIVE ABSORBER AMOUNTS
    DOUBLE PRECISION pb(kdlon, ninter, kflev+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
    DOUBLE PRECISION pbint(kdlon, kflev+1) ! HALF-LEVEL PLANCK FUNCTIONS
    DOUBLE PRECISION pbsur(kdlon, ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION pbsuin(kdlon) ! SURFACE PLANCK FUNCTION
    DOUBLE PRECISION pbtop(kdlon, ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
    DOUBLE PRECISION pdbsl(kdlon, ninter, kflev*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
    DOUBLE PRECISION pemis(kdlon) ! SURFACE EMISSIVITY
    DOUBLE PRECISION ppmb(kdlon, kflev+1) ! HALF-LEVEL PRESSURE (MB)
    DOUBLE PRECISION pga(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgb(kdlon, 8, 2, kflev) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgasur(kdlon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgbsur(kdlon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgatop(kdlon, 8, 2) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgbtop(kdlon, 8, 2) ! PADE APPROXIMANTS

    DOUBLE PRECISION pcntrb(kdlon, kflev+1, kflev+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
    DOUBLE PRECISION pcts(kdlon, kflev) ! COOLING-TO-SPACE TERM
    DOUBLE PRECISION pfluc(kdlon, 2, kflev+1) ! CLEAR-SKY RADIATIVE FLUXES
    ! -----------------------------------------------------------------------
    ! LOCAL VARIABLES:
    DOUBLE PRECISION zadjd(kdlon, kflev+1)
    DOUBLE PRECISION zadju(kdlon, kflev+1)
    DOUBLE PRECISION zdbdt(kdlon, ninter, kflev)
    DOUBLE PRECISION zdisd(kdlon, kflev+1)
    DOUBLE PRECISION zdisu(kdlon, kflev+1)

    INTEGER jk, jl
    ! -----------------------------------------------------------------------

    DO jk = 1, kflev + 1
       DO jl = 1, kdlon
          zadjd(jl, jk) = 0.
          zadju(jl, jk) = 0.
          zdisd(jl, jk) = 0.
          zdisu(jl, jk) = 0.
       END DO
    END DO

    DO jk = 1, kflev
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
