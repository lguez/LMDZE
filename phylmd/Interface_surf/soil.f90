module soil_m

  IMPLICIT NONE

  private fz, compute_beta

contains

  SUBROUTINE soil(nisurf, snow, tsurf, tsoil, soilcap, soilflux)

    ! From LMDZ4/libf/phylmd/soil.F, version 1.1.1.1, 2004 May 19th

    ! Author: Frederic Hourdin, January 30th, 1992

    ! Object: computation of the soil temperature evolution, the heat
    ! capacity per unit surface and the surface conduction flux

    ! Method: implicit time integration

    ! Consecutive ground temperatures are related by:
    ! T(k + 1) = BETA(k) + ALPHA(k) * T(k) (equation 1)
    ! The coefficients BETA and ALPHA are computed at the t - dt time-step.
    ! Structure of the procedure:
    ! 1) BETA and ALPHA coefficients are computed from the old temperature
    ! 2) new temperatures are computed using equation 1
    ! 3) BETA and ALPHA coefficients are computed from the new temperature
    ! profile for the t + dt time-step
    ! 4) the coefficients A and B are computed where the diffusive
    ! fluxes at the t + dt time-step is given by
    ! Fdiff = A + B Ts(t + dt)
    ! or 
    ! Fdiff = F0 + Soilcap (Ts(t + dt) - Ts(t)) / dt
    ! with 
    ! F0 = A + B (Ts(t))
    ! Soilcap = B * dt

    use nr_util, only: pi

    use comconst, only: dtphys
    USE dimsoil, only: nsoilmx
    USE indicesol, only: is_lic, is_oce, is_sic, is_ter
    USE suphec_m, only: rtt
    use unit_nml_m, only: unit_nml

    INTEGER, intent(in):: nisurf ! surface type index
    REAL, intent(in):: snow(:) ! (knon)
    REAL, intent(in):: tsurf(:) ! (knon) surface temperature at time-step t (K)

    real, intent(inout):: tsoil(:, :) ! (knon, nsoilmx)
    ! temperature inside the ground (K), layer 1 nearest to the surface

    REAL, intent(out):: soilcap(:) ! (knon)
    ! specific heat per unit surface (in J m-2 K-1)

    REAL, intent(out):: soilflux(:) ! (knon) 
    ! surface diffusive flux from ground (W m-2)

    ! Local:

    INTEGER jk
    REAL min_period ! no dimension

    REAL fz1
    ! e-folding depth for a wave of period "min_period times 1 s",
    ! divided by e-folding depth for a wave of period one day

    real depth_ratio ! rapport entre les \'epaisseurs de 2 couches successives
    real, save:: delta(nsoilmx - 1)
    REAL therm_i(size(tsurf)) ! (knon) thermal inertia, in W m-2 K-1
    REAL tempor
    REAL, save:: d(nsoilmx - 1), dz2(nsoilmx), c(nsoilmx)
    REAL beta(size(tsurf), nsoilmx - 1) ! (knon, nsoilmx - 1)
    REAL, save:: alpha(nsoilmx - 1)
    REAL, save:: mu
    LOGICAL:: first_call = .TRUE.

    REAL:: inertie_sol = 2000., inertie_sno = 2000., inertie_sic = 2000.
    ! in W m-2 K-1

    namelist /soil_nml/ min_period, depth_ratio, inertie_sol, inertie_sno, &
         inertie_sic

    !-----------------------------------------------------------------------

    IF (first_call) THEN
       ! Compute ground levels:

       ! Default values:
       min_period = 1800.
       depth_ratio = 2.

       print *, "Enter namelist 'soil_nml'."
       read (unit = *, nml = soil_nml)
       write(unit_nml, nml = soil_nml)
       fz1 = sqrt(min_period / pi)

       forall (jk = 1:nsoilmx) dz2(jk) = fz(real(jk), depth_ratio, fz1) &
            - fz(jk - 1., depth_ratio, fz1)

       ! Hourdin 1992 k1078, equation A.12:
       forall (jk = 1:nsoilmx - 1) d(jk) = 1. &
            / (fz(jk + 0.5, depth_ratio, fz1) - fz(jk - 0.5, depth_ratio, fz1))

       ! Hourdin 1992 k1078, equation A.28:
       mu = fz(0.5, depth_ratio, fz1) * d(1)

       c = dz2 / dtphys
       ! Hourdin 1992 k1078, equation A.18:
       delta(nsoilmx - 1) = c(nsoilmx) + d(nsoilmx - 1)

       ! Hourdin 1992 k1078, equation A.16:
       alpha(nsoilmx - 1) = d(nsoilmx - 1) / delta(nsoilmx - 1)

       DO jk = nsoilmx - 1, 2, - 1
          ! Hourdin 1992 k1078, equation A.21:
          delta(jk - 1) = c(jk) + d(jk - 1) + d(jk) * (1. - alpha(jk))
          
          ! Hourdin 1992 k1078, equation A.19:
          alpha(jk - 1) = d(jk - 1) / delta(jk - 1)
       END DO

       first_call = .FALSE.
    END IF

    ! Calcul de l'inertie thermique. On initialise \`a inertie_sic m\^eme
    ! au-dessus d'un point de mer pour le cas o\`u le point de mer
    ! deviendrait point de glace au pas suivant. On corrige si on a un
    ! point de terre avec ou sans glace.

    select case (nisurf)
    case (is_sic)
       where (snow > 0.)
          therm_i = inertie_sno
       elsewhere
          therm_i = inertie_sic
       end where
    case (is_lic)
       where (snow > 0.)
          therm_i = inertie_sno
       elsewhere
          therm_i = inertie_sic
       end where
    case (is_ter)
       where (snow > 0.)
          therm_i = inertie_sno
       elsewhere
          therm_i = inertie_sol
       end where
    case (is_oce)
       therm_i = inertie_sic
    case default
       PRINT *, 'soil: unexpected subscript value:', nisurf
       STOP 1
    END select

    beta = compute_beta(c, d, delta, tsoil)

    ! Hourdin 1992 k1078, equation A.34:
    tsoil(:, 1) = (mu * beta(:, 1) + tsurf) / (mu * (1. - alpha(1)) + 1.)

    DO jk = 1, nsoilmx - 1
       ! Hourdin 1992 k1078, equation A.15:
       tsoil(:, jk + 1) = beta(:, jk) + alpha(jk) * tsoil(:, jk)
    END DO

    IF (nisurf == is_sic) tsoil(:, nsoilmx) = rtt - 1.8
    beta = compute_beta(c, d, delta, tsoil)

    ! Computation of the surface diffusive flux from ground and
    ! calorific capacity of the ground:

    tempor = mu * (1. - alpha(1)) + 1.

    ! Hourdin 1992 k1078, equation A.30:
    soilcap = therm_i * (dz2(1) + dtphys * (1. - alpha(1)) * d(1)) / tempor

    ! Hourdin 1992 k1078, equation A.31:
    soilflux = therm_i * d(1) * (beta(:, 1) + (alpha(1) - 1.) * tsoil(:, 1)) &
         + soilcap * (tsoil(:, 1) * tempor - mu * beta(:, 1) - tsurf) / dtphys

  END SUBROUTINE soil

  !****************************************************************

  pure real function fz(rk, depth_ratio, fz1)

    real, intent(in):: rk

    real, intent(in):: depth_ratio
    ! rapport entre les \'epaisseurs de 2 couches successives

    real, intent(in):: fz1
    ! e-folding depth for a wave of period "min_period times 1 s",
    ! divided by e-folding depth for a wave of period one day

    !-----------------------------------------

    fz = fz1 * (depth_ratio**rk - 1.) / (depth_ratio - 1.)
    ! Hourdin 1992 k1078, equation A.5

  end function fz

  !****************************************************************

  pure function compute_beta(c, d, delta, tsoil) result (beta)

    ! Computation of the coefficient Beta for the next step.

    USE dimsoil, only: nsoilmx
    
    REAL, intent(in):: c(:) ! (nsoilmx)
    REAL, intent(in):: d(:) ! (nsoilmx - 1)
    real, intent(in):: delta(:) ! (nsoilmx - 1)

    real, intent(in):: tsoil(:, :) ! (knon, nsoilmx)
    ! temperature inside the ground (K), layer 1 nearest to the surface

    REAL beta(size(tsoil, 1), size(d)) ! (knon, nsoilmx - 1)

    ! Local:
    integer jk

    !------------------------------------------------------------------

    ! Hourdin 1992 k1078, equation A.17:
    beta(:, nsoilmx - 1) = c(nsoilmx) * tsoil(:, nsoilmx) / delta(nsoilmx - 1)

    DO jk = nsoilmx - 1, 2, - 1
       ! Hourdin 1992 k1078, equation A.20:
       beta(:, jk - 1) = (tsoil(:, jk) * c(jk) + d(jk) * beta(:, jk)) &
            / delta(jk - 1)
    END DO

  end function compute_beta

end module soil_m
