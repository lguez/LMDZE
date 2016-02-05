module soil_m

  IMPLICIT NONE

contains

  SUBROUTINE soil(dtime, nisurf, knon, snow, tsurf, tsoil, soilcap, soilflux)

    ! From LMDZ4/libf/phylmd/soil.F, version 1.1.1.1 2004/05/19

    ! Author: Frederic Hourdin 30/01/92

    ! Object: computation of the soil temperature evolution, the
    ! surfacic heat capacity "Soilcap" and the surface conduction flux

    ! Method: implicit time integration

    ! Consecutive ground temperatures are related by:
    ! T(k+1) = C(k) + D(k)*T(k) (1)
    ! the coefficients C and D are computed at the t-dt time-step.
    ! Routine structure:
    ! 1) new temperatures are computed using (1)
    ! 2) C and D coefficients are computed from the new temperature
    ! profile for the t+dt time-step
    ! 3) the coefficients A and B are computed where the diffusive
    ! fluxes at the t+dt time-step is given by
    ! Fdiff = A + B Ts(t+dt)
    ! or Fdiff = F0 + Soilcap (Ts(t+dt)-Ts(t))/dt
    ! with F0 = A + B (Ts(t))
    ! Soilcap = B*dt

    USE dimens_m, only:
    USE indicesol
    USE dimphy
    USE dimsoil
    USE suphec_m

    ! Interface:
    ! ----------

    ! Arguments:
    ! ----------
    ! dtime physical timestep (s)
    ! indice sub-surface index
    ! snow(klon, nbsrf) snow
    ! tsurf(knon) surface temperature at time-step t (K)
    ! tsoil(klon, nsoilmx) temperature inside the ground (K)
    ! soilcap(klon) surfacic specific heat (W*m-2*s*K-1)
    ! soilflux(klon) surface diffusive flux from ground (Wm-2)

    ! declarations:
    ! -------------

    ! -----------------------------------------------------------------------
    ! arguments
    ! ---------

    REAL dtime
    INTEGER nisurf, knon
    REAL tsurf(knon), tsoil(klon, nsoilmx), snow(klon)
    REAL soilcap(klon), soilflux(klon)

    ! -----------------------------------------------------------------------
    ! local arrays
    ! ------------

    INTEGER ig, jk
    ! $$$ REAL zdz2(nsoilmx), z1(klon)
    REAL zdz2(nsoilmx), z1(klon, nbsrf)
    REAL min_period, dalph_soil
    REAL ztherm_i(klon)

    ! local saved variables:
    ! ----------------------
    REAL dz1(nsoilmx), dz2(nsoilmx)
    ! $$$ REAL zc(klon, nsoilmx), zd(klon, nsoilmx)
    REAL zc(klon, nsoilmx, nbsrf), zd(klon, nsoilmx, nbsrf)
    REAL lambda
    SAVE dz1, dz2, zc, zd, lambda
    LOGICAL firstsurf(nbsrf)
    SAVE firstsurf
    REAL isol, isno, iice
    SAVE isol, isno, iice

    DATA firstsurf/.TRUE., .TRUE., .TRUE., .TRUE./

    DATA isol, isno, iice/2000., 2000., 2000./

    ! -----------------------------------------------------------------------
    ! Depthts:
    ! --------

    REAL rk, fz1, rk1, rk2

    soilflux(:) = 0.
    ! calcul de l'inertie thermique a partir de la variable rnat.
    ! on initialise a iice meme au-dessus d'un point de mer au cas
    ! ou le point de mer devienne point de glace au pas suivant
    ! on corrige si on a un point de terre avec ou sans glace

    IF (nisurf==is_sic) THEN
       DO ig = 1, knon
          ztherm_i(ig) = iice
          IF (snow(ig)>0.0) ztherm_i(ig) = isno
       END DO
    ELSE IF (nisurf==is_lic) THEN
       DO ig = 1, knon
          ztherm_i(ig) = iice
          IF (snow(ig)>0.0) ztherm_i(ig) = isno
       END DO
    ELSE IF (nisurf==is_ter) THEN
       DO ig = 1, knon
          ztherm_i(ig) = isol
          IF (snow(ig)>0.0) ztherm_i(ig) = isno
       END DO
    ELSE IF (nisurf==is_oce) THEN
       DO ig = 1, knon
          ztherm_i(ig) = iice
       END DO
    ELSE
       PRINT *, 'valeur d indice non prevue', nisurf
       STOP 1
    END IF

    IF (firstsurf(nisurf)) THEN

       ! -----------------------------------------------------------------------
       ! ground levels
       ! grnd=z/l where l is the skin depth of the diurnal cycle:
       ! --------------------------------------------------------

       min_period = 1800. ! en secondes
       dalph_soil = 2. ! rapport entre les epaisseurs de 2 couches succ.

       OPEN (99, FILE='soil.def', STATUS='old', FORM='formatted', ERR=9999)
       READ (99, *) min_period
       READ (99, *) dalph_soil
       PRINT *, 'Discretization for the soil model'
       PRINT *, 'First level e-folding depth', min_period, ' dalph', &
            dalph_soil
       CLOSE (99)
9999   CONTINUE

       ! la premiere couche represente un dixieme de cycle diurne
       fz1 = sqrt(min_period/3.14)

       DO jk = 1, nsoilmx
          rk1 = jk
          rk2 = jk - 1
          dz2(jk) = fz(rk1) - fz(rk2)
       END DO
       DO jk = 1, nsoilmx - 1
          rk1 = jk + .5
          rk2 = jk - .5
          dz1(jk) = 1./(fz(rk1)-fz(rk2))
       END DO
       lambda = fz(.5)*dz1(1)
       PRINT *, 'full layers, intermediate layers (seconds)'
       DO jk = 1, nsoilmx
          rk = jk
          rk1 = jk + .5
          rk2 = jk - .5
          PRINT *, 'fz=', fz(rk1)*fz(rk2)*3.14, fz(rk)*fz(rk)*3.14
       END DO
       ! PB
       firstsurf(nisurf) = .FALSE.

       ! Initialisations:
       ! ----------------

    ELSE
       ! -----------------------------------------------------------------------
       ! Computation of the soil temperatures using the Cgrd and Dgrd
       ! coefficient computed at the previous time-step:
       ! -----------------------------------------------

       ! surface temperature
       DO ig = 1, knon
          tsoil(ig, 1) = (lambda*zc(ig, 1, nisurf)+tsurf(ig))/(lambda*(1.-zd(ig, 1, &
               nisurf))+1.)
       END DO

       ! other temperatures
       DO jk = 1, nsoilmx - 1
          DO ig = 1, knon
             tsoil(ig, jk+1) = zc(ig, jk, nisurf) + zd(ig, jk, nisurf)*tsoil(ig, &
                  jk)
          END DO
       END DO

    END IF
    ! -----------------------------------------------------------------------
    ! Computation of the Cgrd and Dgrd coefficient for the next step:
    ! ---------------------------------------------------------------

    ! $$$ PB ajout pour cas glace de mer
    IF (nisurf==is_sic) THEN
       DO ig = 1, knon
          tsoil(ig, nsoilmx) = rtt - 1.8
       END DO
    END IF

    DO jk = 1, nsoilmx
       zdz2(jk) = dz2(jk)/dtime
    END DO

    DO ig = 1, knon
       z1(ig, nisurf) = zdz2(nsoilmx) + dz1(nsoilmx-1)
       zc(ig, nsoilmx-1, nisurf) = zdz2(nsoilmx)*tsoil(ig, nsoilmx)/ &
            z1(ig, nisurf)
       zd(ig, nsoilmx-1, nisurf) = dz1(nsoilmx-1)/z1(ig, nisurf)
    END DO

    DO jk = nsoilmx - 1, 2, -1
       DO ig = 1, knon
          z1(ig, nisurf) = 1./(zdz2(jk)+dz1(jk-1)+dz1(jk)*(1.-zd(ig, jk, nisurf)))
          zc(ig, jk-1, nisurf) = (tsoil(ig, jk)*zdz2(jk)+dz1(jk)*zc(ig, jk, nisurf) &
               )*z1(ig, nisurf)
          zd(ig, jk-1, nisurf) = dz1(jk-1)*z1(ig, nisurf)
       END DO
    END DO

    ! -----------------------------------------------------------------------
    ! computation of the surface diffusive flux from ground and
    ! calorific capacity of the ground:
    ! ---------------------------------

    DO ig = 1, knon
       soilflux(ig) = ztherm_i(ig)*dz1(1)*(zc(ig, 1, nisurf)+(zd(ig, 1, &
            nisurf)-1.)*tsoil(ig, 1))
       soilcap(ig) = ztherm_i(ig)*(dz2(1)+dtime*(1.-zd(ig, 1, nisurf))*dz1(1))
       z1(ig, nisurf) = lambda*(1.-zd(ig, 1, nisurf)) + 1.
       soilcap(ig) = soilcap(ig)/z1(ig, nisurf)
       soilflux(ig) = soilflux(ig) + soilcap(ig)*(tsoil(ig, 1)*z1(ig, nisurf)- &
            lambda*zc(ig, 1, nisurf)-tsurf(ig))/dtime
    END DO

  contains

    real function fz(rk)
      real rk
      fz = fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)
    end function fz

  END SUBROUTINE soil

end module soil_m
