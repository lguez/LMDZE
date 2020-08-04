module soil_m

  IMPLICIT NONE

  private fz, compute_c_d

contains

  SUBROUTINE soil(nisurf, snow, tsurf, tsoil, soilcap, soilflux)

    ! From LMDZ4/libf/phylmd/soil.F, version 1.1.1.1, 2004 May 19th

    ! Author: Frederic Hourdin, January 30th, 1992

    ! Object: computation of the soil temperature evolution, the heat
    ! capacity per unit surface and the surface conduction flux

    ! Method: implicit time integration

    ! Consecutive ground temperatures are related by:
    ! T(k + 1) = C(k) + D(k) * T(k) (equation 1)
    ! The coefficients C and D are computed at the t - dt time-step.
    ! Structure of the procedure:
!   ! 1) C and D coefficients are computed from the old temperature
    ! 2) new temperatures are computed using equation 1
    ! 3) C and D coefficients are computed from the new temperature
    ! profile for the t + dt time-step
    ! 4) the coefficients A and B are computed where the diffusive
    ! fluxes at the t + dt time-step is given by
    ! Fdiff = A + B Ts(t + dt)
    ! or 
    ! Fdiff = F0 + Soilcap (Ts(t + dt) - Ts(t)) / dt
    ! with 
    ! F0 = A + B (Ts(t))
    ! Soilcap = B * dt

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
    
    INTEGER knon, ig, jk
    REAL zdz2(nsoilmx), z1
    REAL min_period ! in s
    real depth_ratio ! rapport entre les \'epaisseurs de 2 couches successives
    REAL ztherm_i(size(tsurf)) ! (knon)
    REAL, save:: dz1(nsoilmx - 1), dz2(nsoilmx)
    REAL zc(size(tsurf), nsoilmx), zd(size(tsurf), nsoilmx) ! (knon, nsoilmx)
    REAL, save:: mu
    LOGICAL:: first_call = .TRUE.
    REAL, parameter:: isol = 2000., isno = 2000., iice = 2000.
    REAL fz1 ! depth

    namelist /soil_nml/ min_period, depth_ratio

    !-----------------------------------------------------------------------

    knon = size(tsurf)

    IF (first_call) THEN
       ! ground levels
       ! z / l where l is the skin depth of the diurnal cycle

       ! Default values:
       min_period = 1800.
       depth_ratio = 2.
       
       print *, "Enter namelist 'soil_nml'."
       read (unit = *, nml = soil_nml)
       write(unit_nml, nml = soil_nml)

       ! La premi\`ere couche repr\'esente un dixi\`eme de cycle diurne :
       fz1 = sqrt(min_period / 3.14)

       forall (jk = 1:nsoilmx) dz2(jk) = fz(real(jk), depth_ratio, fz1) &
            - fz(jk - 1., depth_ratio, fz1)
       forall (jk = 1:nsoilmx - 1) dz1(jk) = 1. &
            / (fz(jk + 0.5, depth_ratio, fz1) - fz(jk - 0.5, depth_ratio, fz1))
       mu = fz(0.5, depth_ratio, fz1) * dz1(1)
       first_call = .FALSE.
    END IF

    ! Calcul de l'inertie thermique. On initialise \`a iice m\^eme
    ! au-dessus d'un point de mer pour le cas o\`u le point de mer
    ! deviendrait point de glace au pas suivant. On corrige si on a un
    ! point de terre avec ou sans glace.

    select case (nisurf)
    case (is_sic)
       DO ig = 1, knon
          IF (snow(ig) > 0.) then
             ztherm_i(ig) = isno
          else
             ztherm_i(ig) = iice
          end IF
       END DO
    case (is_lic)
       DO ig = 1, knon
          IF (snow(ig) > 0.) then
             ztherm_i(ig) = isno
          else
             ztherm_i(ig) = iice
          end IF
       END DO
    case (is_ter)
       DO ig = 1, knon
          IF (snow(ig) > 0.) then
             ztherm_i(ig) = isno
          else
             ztherm_i(ig) = isol
          end IF
       END DO
    case (is_oce)
       DO ig = 1, knon
          ztherm_i(ig) = iice
       END DO
    case default
       PRINT *, 'soil: unexpected subscript value:', nisurf
       STOP 1
    END select

    DO jk = 1, nsoilmx
       zdz2(jk) = dz2(jk) / dtphys
    END DO

    call compute_c_d(zdz2, dz1, zc, zd, tsoil)

    ! Computation of the soil temperatures using the Zc and Zd
    ! coefficient computed above:

    ! Surface temperature (Hourdin 1992 k1078, equation A.27):
    tsoil(:, 1) = (mu * zc(:, 1) + tsurf(:)) / (mu * (1. - zd(:, 1)) + 1.)

    ! Other temperatures:
    DO jk = 1, nsoilmx - 1
       DO ig = 1, knon
          tsoil(ig, jk + 1) = zc(ig, jk) + zd(ig, jk) * tsoil(ig, jk)
       END DO
    END DO

    IF (nisurf == is_sic) THEN
       DO ig = 1, knon
          tsoil(ig, nsoilmx) = rtt - 1.8
       END DO
    END IF

    call compute_c_d(zdz2, dz1, zc, zd, tsoil)

    ! Computation of the surface diffusive flux from ground and
    ! calorific capacity of the ground:

    DO ig = 1, knon
       soilflux(ig) = ztherm_i(ig) * dz1(1) * (zc(ig, 1) &
            + (zd(ig, 1) - 1.) * tsoil(ig, 1))
       soilcap(ig) = ztherm_i(ig) * (dz2(1) &
            + dtphys * (1. - zd(ig, 1)) * dz1(1))
       z1 = mu * (1. - zd(ig, 1)) + 1.
       soilcap(ig) = soilcap(ig) / z1
       soilflux(ig) = soilflux(ig) + soilcap(ig) * (tsoil(ig, 1) * z1 - mu &
            * zc(ig, 1) - tsurf(ig)) / dtphys
    END DO

  END SUBROUTINE soil

  !****************************************************************

  pure real function fz(rk, depth_ratio, fz1)

    real, intent(in):: rk

    real, intent(in):: depth_ratio
    ! rapport entre les \'epaisseurs de 2 couches successives

    real, intent(in):: fz1 ! depth

    !-----------------------------------------

    fz = fz1 * (depth_ratio**rk - 1.) / (depth_ratio - 1.)
    ! Hourdin 1992 k1078, equation A.5

  end function fz

  !****************************************************************

  subroutine compute_c_d(zdz2, dz1, zc, zd, tsoil)

    ! Computation of the coefficients Zc and Zd for the next step.

    USE dimsoil, only: nsoilmx
    
    REAL, intent(in):: zdz2(:) ! (nsoilmx)
    REAL, intent(in):: dz1(:) ! (nsoilmx - 1)
    REAL, intent(inout):: zc(:, :), zd(:, :) ! (knon, nsoilmx)

    real, intent(in):: tsoil(:, :) ! (knon, nsoilmx)
    ! temperature inside the ground (K), layer 1 nearest to the surface

    ! Local:
    integer ig, knon, jk
    real z1

    !------------------------------------------------------------------

    knon  = size(tsoil, 1)
    z1 = zdz2(nsoilmx) + dz1(nsoilmx - 1)
    zc(:, nsoilmx - 1) = zdz2(nsoilmx) * tsoil(:, nsoilmx) / z1
    zd(:, nsoilmx - 1) = dz1(nsoilmx - 1) / z1

    DO jk = nsoilmx - 1, 2, - 1
       DO ig = 1, knon
          z1 = 1. / (zdz2(jk) + dz1(jk - 1) + dz1(jk) * (1. - zd(ig, jk)))
          zc(ig, jk - 1) = (tsoil(ig, jk) * zdz2(jk) &
               + dz1(jk) * zc(ig, jk)) * z1
          zd(ig, jk - 1) = dz1(jk - 1) * z1
       END DO
    END DO

  end subroutine compute_c_d

end module soil_m
