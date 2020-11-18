module soil_m

  IMPLICIT NONE

  private compute_beta

contains

  SUBROUTINE soil(nisurf, snow, tsurf, tsoil, soilcap, soilflux)

    ! From LMDZ4/libf/phylmd/soil.F, version 1.1.1.1, 2004 May 19th

    ! Author: Frederic Hourdin, January 30th, 1992

    ! Object: computation of the soil temperature evolution, the heat
    ! capacity per unit surface and the surface conduction flux

    ! Method: implicit time integration

    ! Consecutive ground temperatures are related by:

    ! T(k + 1) = BETA(k) + ALPHA(k) * T(k) (equation 1)

    ! The coefficient BETA is computed using soil temperatures from
    ! the previous time-step.

    ! Structure of the procedure:
    
    ! 1) BETA is computed from input soil temperature which is
    ! supposed to be the soil temperature profile two time steps
    ! before. This is $\beta^{t - \delta t}$.

    ! 2) New soil temperature is computed using equation 1. This is
    ! supposed to be the soil temperature profile for one time step
    ! before, since it is associated to the surface temperature for
    ! one time step before.
    
    ! 3) BETA is computed again from the new soil temperature
    ! profile. This is $\beta^t$.
    
    ! 4) The coefficients soilflux and Soilcap are computed in the
    ! following relation between surface thermal conduction flux at
    ! time-step t, Fc0, and surface temperature at t:
    
    ! Fc0 + soilflux  = Soilcap (Ts(t) - Ts(t - delta t)) / (delta t)

    use nr_util, only: pi

    use comconst, only: dtphys
    USE dimsoil, only: nsoilmx
    USE indicesol, only: is_lic, is_oce, is_sic, is_ter
    USE suphec_m, only: rtt
    use unit_nml_m, only: unit_nml

    INTEGER, intent(in):: nisurf ! surface type index
    REAL, intent(in):: snow(:) ! (knon)

    REAL, intent(in):: tsurf(:) ! (knon)
    ! surface temperature at previous time-step (K)

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
    ! normalized e-folding depth for a wave of period "min_period times 1 s"

    real depth_ratio ! rapport entre les \'epaisseurs de 2 couches successives
    real, save:: delta(nsoilmx - 1)
    REAL therm_i(size(tsurf)) ! (knon) thermal inertia, in W m-2 K-1
    REAL, save:: tempor
    REAL, save:: d(nsoilmx - 1), c(nsoilmx)
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
       print *, "fz1 = ", fz1

       ! Hourdin 1992 k1078, equation A.11:
       
       c(1) = fz1 / dtphys

       do jk = 2, nsoilmx
          c(jk) = c(jk - 1) * depth_ratio
       end do

       ! Hourdin 1992 k1078, equation A.12:

       d(1) = 1. / (sqrt(depth_ratio) * fz1)

       do jk = 2, nsoilmx - 1
          d(jk) = d(jk - 1) / depth_ratio
       end do

       ! Hourdin 1992 k1078, equation A.28:
       mu = 1. / (depth_ratio + sqrt(depth_ratio))
       print *, "mu = ", mu

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

       tempor = mu * (1. - alpha(1)) + 1.
       print *, "tempor = ", tempor
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

    ! Hourdin 1992 k1078, equation A.30:
    soilcap = therm_i * dtphys * (c(1) + (1. - alpha(1)) * d(1)) / tempor

    ! Hourdin 1992 k1078, equation A.31:
    soilflux = therm_i * d(1) * (beta(:, 1) + (alpha(1) - 1.) * tsoil(:, 1)) &
         + soilcap * (tsoil(:, 1) * tempor - mu * beta(:, 1) - tsurf) / dtphys

  END SUBROUTINE soil

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
