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
    ! 1) new temperatures are computed using equation 1
    ! 2) C and D coefficients are computed from the new temperature
    ! profile for the t + dt time-step
    ! 3) the coefficients A and B are computed where the diffusive
    ! fluxes at the t + dt time-step is given by
    ! Fdiff = A + B Ts(t + dt)
    ! or 
    ! Fdiff = F0 + Soilcap (Ts(t + dt) - Ts(t)) / dt
    ! with 
    ! F0 = A + B (Ts(t))
    ! Soilcap = B * dt

    ! Libraries:
    use jumble, only: new_unit

    use comconst, only: dtphys
    USE dimphy, only: klon
    USE dimsoil, only: nsoilmx
    USE indicesol, only: nbsrf, is_lic, is_oce, is_sic, is_ter
    USE suphec_m, only: rtt

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
    INTEGER knon, ig, jk, unit
    REAL zdz2(nsoilmx)
    real z1(size(tsurf), nbsrf) ! (knon, nbsrf)
    REAL min_period ! in s
    real dalph_soil ! rapport entre les \'epaisseurs de 2 couches successives
    REAL ztherm_i(size(tsurf)) ! (knon)
    REAL, save:: dz1(nsoilmx), dz2(nsoilmx)
    REAL, allocatable, save:: zc(:, :, :), zd(:, :, :) ! (klon, nsoilmx, nbsrf)
    REAL, save:: lambda
    LOGICAL:: firstsurf(nbsrf) = .TRUE.
    REAL, parameter:: isol = 2000., isno = 2000., iice = 2000.
    REAL rk, fz1, rk1, rk2 ! depths

    !-----------------------------------------------------------------------

    knon = size(tsurf)

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

    if (all(firstsurf)) allocate(zc(klon, nsoilmx, nbsrf), &
         zd(klon, nsoilmx, nbsrf))
    
    IF (firstsurf(nisurf)) THEN
       ! ground levels
       ! grnd=z / l where l is the skin depth of the diurnal cycle:

       min_period = 1800.
       dalph_soil = 2.
       call new_unit(unit)
       OPEN(unit, FILE = 'soil.def', STATUS = 'old', action = "read", &
            position = 'rewind', ERR = 9999)
       READ(unit, fmt = *) min_period
       READ(unit, fmt = *) dalph_soil
       PRINT *, 'Discretization for the soil model'
       PRINT *, 'First level e-folding depth', min_period, ' dalph', &
            dalph_soil
       CLOSE(unit)
9999   CONTINUE

       ! La premi\`ere couche repr\'esente un dixi\`eme de cycle diurne :
       fz1 = sqrt(min_period / 3.14)

       DO jk = 1, nsoilmx
          rk1 = jk
          rk2 = jk - 1
          dz2(jk) = fz(rk1, dalph_soil, fz1) - fz(rk2, dalph_soil, fz1)
       END DO

       DO jk = 1, nsoilmx - 1
          rk1 = jk + .5
          rk2 = jk - .5
          dz1(jk) = 1. / (fz(rk1, dalph_soil, fz1) - fz(rk2, dalph_soil, fz1))
       END DO

       lambda = fz(.5, dalph_soil, fz1) * dz1(1)
       PRINT *, 'full layers, intermediate layers (seconds)'

       DO jk = 1, nsoilmx
          rk = jk
          rk1 = jk + .5
          rk2 = jk - .5
          PRINT *, 'fz=', fz(rk1, dalph_soil, fz1) * fz(rk2, dalph_soil, fz1) &
               * 3.14, fz(rk, dalph_soil, fz1) * fz(rk, dalph_soil, fz1) * 3.14
       END DO

       firstsurf(nisurf) = .FALSE.
    ELSE
       ! Computation of the soil temperatures using the Zc and Zd
       ! coefficient computed at the previous time-step:

       ! Surface temperature:
       DO ig = 1, knon
          tsoil(ig, 1) = (lambda * zc(ig, 1, nisurf) + tsurf(ig)) &
               / (lambda * (1. - zd(ig, 1, nisurf)) + 1.)
       END DO

       ! Other temperatures:
       DO jk = 1, nsoilmx - 1
          DO ig = 1, knon
             tsoil(ig, jk + 1) = zc(ig, jk, nisurf) &
                  + zd(ig, jk, nisurf) * tsoil(ig, jk)
          END DO
       END DO
    END IF

    IF (nisurf==is_sic) THEN
       DO ig = 1, knon
          tsoil(ig, nsoilmx) = rtt - 1.8
       END DO
    END IF

    DO jk = 1, nsoilmx
       zdz2(jk) = dz2(jk) / dtphys
    END DO

    call compute_c_d(zdz2, dz1, zc(:knon, :, nisurf), zd(:knon, :, nisurf), &
         tsoil)

    ! Computation of the surface diffusive flux from ground and
    ! calorific capacity of the ground:

    DO ig = 1, knon
       soilflux(ig) = ztherm_i(ig) * dz1(1) * (zc(ig, 1, nisurf) &
            + (zd(ig, 1, nisurf) - 1.) * tsoil(ig, 1))
       soilcap(ig) = ztherm_i(ig) * (dz2(1) &
            + dtphys * (1. - zd(ig, 1, nisurf)) * dz1(1))
       z1(ig, nisurf) = lambda * (1. - zd(ig, 1, nisurf)) + 1.
       soilcap(ig) = soilcap(ig) / z1(ig, nisurf)
       soilflux(ig) = soilflux(ig) + soilcap(ig) * (tsoil(ig, 1) &
            * z1(ig, nisurf) - lambda * zc(ig, 1, nisurf) - tsurf(ig)) / dtphys
    END DO

  END SUBROUTINE soil

  !****************************************************************

  pure real function fz(rk, dalph_soil, fz1)

    real, intent(in):: rk

    real, intent(in):: dalph_soil
    ! rapport entre les \'epaisseurs de 2 couches successives

    real, intent(in):: fz1 ! depth

    !-----------------------------------------

    fz = fz1 * (dalph_soil**rk - 1.) / (dalph_soil - 1.)

  end function fz

  !****************************************************************

  subroutine compute_c_d(zdz2, dz1, zc, zd, tsoil)

    ! Computation of the coefficients Zc and Zd for the next step.

    USE dimsoil, only: nsoilmx
    
    REAL, intent(in):: zdz2(:) ! (nsoilmx)
    REAL, intent(in):: dz1(:) ! (nsoilmx)
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
