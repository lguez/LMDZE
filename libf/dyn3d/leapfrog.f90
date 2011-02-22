module leapfrog_m

  IMPLICIT NONE

contains

  SUBROUTINE leapfrog(ucov, vcov, teta, ps, masse, phis, q, time_0)

    ! From dyn3d/leapfrog.F, version 1.6, 2005/04/13 08:58:34
    ! Authors: P. Le Van, L. Fairhead, F. Hourdin
    ! Matsuno-leapfrog scheme.

    use addfi_m, only: addfi
    use bilan_dyn_m, only: bilan_dyn
    use caladvtrac_m, only: caladvtrac
    USE calfis_m, ONLY: calfis
    USE com_io_dyn, ONLY: histaveid
    USE comconst, ONLY: daysec, dtphys, dtvr
    USE comgeom, ONLY: aire_2d, apoln, apols
    USE comvert, ONLY: ap, bp
    USE conf_gcm_m, ONLY: day_step, iconser, iperiod, iphysiq, nday, offline, &
         periodav
    USE dimens_m, ONLY: iim, jjm, llm, nqmx
    USE dynetat0_m, ONLY: day_ini
    use dynredem1_m, only: dynredem1
    USE exner_hyb_m, ONLY: exner_hyb
    use filtreg_m, only: filtreg
    USE guide_m, ONLY: guide
    use inidissip_m, only: idissip
    use integrd_m, only: integrd
    USE logic, ONLY: iflag_phys, ok_guide
    USE paramet_m, ONLY: ip1jmp1
    USE pressure_var, ONLY: p3d
    USE temps, ONLY: itau_dyn

    ! Variables dynamiques:
    REAL, intent(inout):: ucov(ip1jmp1, llm) ! vent covariant
    REAL, intent(inout):: vcov((iim + 1) * jjm, llm) ! vent covariant
    REAL, intent(inout):: teta(iim + 1, jjm + 1, llm) ! potential temperature
    REAL, intent(inout):: ps(iim + 1, jjm + 1) ! pression au sol, en Pa
    REAL masse(ip1jmp1, llm) ! masse d'air
    REAL phis(ip1jmp1) ! geopotentiel au sol

    REAL, intent(inout):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    ! mass fractions of advected fields

    REAL, intent(in):: time_0

    ! Variables local to the procedure:

    ! Variables dynamiques:

    REAL pks(ip1jmp1) ! exner au sol
    REAL pk(iim + 1, jjm + 1, llm) ! exner au milieu des couches
    REAL pkf(ip1jmp1, llm) ! exner filt.au milieu des couches
    REAL phi(ip1jmp1, llm) ! geopotential
    REAL w(ip1jmp1, llm) ! vitesse verticale

    ! variables dynamiques intermediaire pour le transport
    REAL pbaru(ip1jmp1, llm), pbarv((iim + 1) * jjm, llm) !flux de masse

    ! variables dynamiques au pas - 1
    REAL vcovm1((iim + 1) * jjm, llm), ucovm1(ip1jmp1, llm)
    REAL tetam1(iim + 1, jjm + 1, llm), psm1(iim + 1, jjm + 1)
    REAL massem1(ip1jmp1, llm)

    ! tendances dynamiques
    REAL dv((iim + 1) * jjm, llm), du(ip1jmp1, llm)
    REAL dteta(ip1jmp1, llm), dq(ip1jmp1, llm, nqmx), dp(ip1jmp1)

    ! tendances de la dissipation
    REAL dvdis((iim + 1) * jjm, llm), dudis(ip1jmp1, llm)
    REAL dtetadis(iim + 1, jjm + 1, llm)

    ! tendances physiques
    REAL dvfi((iim + 1) * jjm, llm), dufi(ip1jmp1, llm)
    REAL dtetafi(ip1jmp1, llm), dqfi(ip1jmp1, llm, nqmx), dpfi(ip1jmp1)

    ! variables pour le fichier histoire

    INTEGER itau ! index of the time step of the dynamics, starts at 0
    INTEGER itaufin
    REAL time ! time of day, as a fraction of day length
    real finvmaold(ip1jmp1, llm)
    INTEGER l
    REAL rdayvrai, rdaym_ini

    ! Variables test conservation energie
    REAL ecin(iim + 1, jjm + 1, llm), ecin0(iim + 1, jjm + 1, llm)
    ! Tendance de la temp. potentiel d (theta) / d t due a la 
    ! tansformation d'energie cinetique en energie thermique
    ! cree par la dissipation
    REAL dtetaecdt(iim + 1, jjm + 1, llm)
    REAL vcont((iim + 1) * jjm, llm), ucont(ip1jmp1, llm)
    logical leapf
    real dt

    !---------------------------------------------------

    print *, "Call sequence information: leapfrog"

    itaufin = nday * day_step
    ! "day_step" is a multiple of "iperiod", therefore "itaufin" is one too

    dq = 0.

    ! On initialise la pression et la fonction d'Exner :
    forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
    CALL exner_hyb(ps, p3d, pks, pk, pkf)

    time_integration: do itau = 0, itaufin - 1
       leapf = mod(itau, iperiod) /= 0
       if (leapf) then
          dt = 2 * dtvr
       else
          ! Matsuno
          dt = dtvr
          if (ok_guide .and. (itaufin - itau - 1) * dtvr > 21600.) &
               call guide(itau, ucov, vcov, teta, q, masse, ps)
          vcovm1 = vcov
          ucovm1 = ucov
          tetam1 = teta
          massem1 = masse
          psm1 = ps
          finvmaold = masse
          CALL filtreg(finvmaold, jjm + 1, llm, - 2, 2, .TRUE., 1)
       end if

       ! Calcul des tendances dynamiques:
       CALL geopot(ip1jmp1, teta, pk, pks, phis, phi)
       CALL caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
            MOD(itau, iconser) == 0, du, dv, dteta, dp, w, pbaru, pbarv, &
            time_0)

       ! Calcul des tendances advection des traceurs (dont l'humidité)
       CALL caladvtrac(q, pbaru, pbarv, p3d, masse, dq, teta, pk)

       ! Stokage du flux de masse pour traceurs offline:
       IF (offline) CALL fluxstokenc(pbaru, pbarv, masse, teta, phi, phis, &
            dtvr, itau)

       ! integrations dynamique et traceurs:
       CALL integrd(vcovm1, ucovm1, tetam1, psm1, massem1, dv, du, dteta, dp, &
            vcov, ucov, teta, q(:, :, :, :2), ps, masse, finvmaold, dt, leapf)

       if (.not. leapf) then
          ! Matsuno backward
          forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
          CALL exner_hyb(ps, p3d, pks, pk, pkf)

          ! Calcul des tendances dynamiques:
          CALL geopot(ip1jmp1, teta, pk, pks, phis, phi)
          CALL caldyn(itau + 1, ucov, vcov, teta, ps, masse, pk, pkf, phis, &
               phi, .false., du, dv, dteta, dp, w, pbaru, pbarv, time_0)

          ! integrations dynamique et traceurs:
          CALL integrd(vcovm1, ucovm1, tetam1, psm1, massem1, dv, du, dteta, &
               dp, vcov, ucov, teta, q(:, :, :, :2), ps, masse, finvmaold, &
               dtvr, leapf=.false.)
       end if

       IF (MOD(itau + 1, iphysiq) == 0 .AND. iflag_phys /= 0) THEN
          ! calcul des tendances physiques:

          forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
          CALL exner_hyb(ps, p3d, pks, pk, pkf)

          rdaym_ini = itau * dtvr / daysec
          rdayvrai = rdaym_ini + day_ini
          time = REAL(mod(itau, day_step)) / day_step + time_0
          IF (time > 1.) time = time - 1.

          CALL calfis(rdayvrai, time, ucov, vcov, teta, q, masse, ps, pk, &
               phis, phi, du, dv, dteta, dq, w, dufi, dvfi, dtetafi, dqfi, &
               dpfi, lafin=itau+1==itaufin)

          ! ajout des tendances physiques:
          CALL addfi(nqmx, dtphys, ucov, vcov, teta, q, ps, dufi, dvfi, &
               dtetafi, dqfi, dpfi)
       ENDIF

       forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
       CALL exner_hyb(ps, p3d, pks, pk, pkf)

       IF (MOD(itau + 1, idissip) == 0) THEN
          ! dissipation horizontale et verticale des petites echelles:

          ! calcul de l'energie cinetique avant dissipation
          call covcont(llm, ucov, vcov, ucont, vcont)
          call enercin(vcov, ucov, vcont, ucont, ecin0)

          ! dissipation
          CALL dissip(vcov, ucov, teta, p3d, dvdis, dudis, dtetadis)
          ucov=ucov + dudis
          vcov=vcov + dvdis

          ! On rajoute la tendance due à la transformation Ec -> E
          ! thermique créée lors de la dissipation
          call covcont(llm, ucov, vcov, ucont, vcont)
          call enercin(vcov, ucov, vcont, ucont, ecin)
          dtetaecdt= (ecin0 - ecin) / pk
          dtetadis=dtetadis + dtetaecdt
          teta=teta + dtetadis

          ! Calcul de la valeur moyenne aux pôles :
          forall (l = 1: llm)
             teta(:, 1, l) = SUM(aire_2d(:iim, 1) * teta(:iim, 1, l)) &
                  / apoln
             teta(:, jjm + 1, l) = SUM(aire_2d(:iim, jjm+1) &
                  * teta(:iim, jjm + 1, l)) / apols
          END forall

          ps(:, 1) = SUM(aire_2d(:iim, 1) * ps(:iim, 1)) / apoln
          ps(:, jjm + 1) = SUM(aire_2d(:iim, jjm+1) * ps(:iim, jjm + 1)) &
               / apols
       END IF

       IF (MOD(itau + 1, iperiod) == 0) THEN
          ! Écriture du fichier histoire moyenne:
          CALL writedynav(histaveid, nqmx, itau + 1, vcov, ucov, teta, pk, &
               phi, q, masse, ps, phis)
          call bilan_dyn(ps, masse, pk, pbaru, pbarv, teta, phi, ucov, vcov, &
               q(:, :, :, 1), dt_app = dtvr * iperiod, &
               dt_cum = dtvr * day_step * periodav)
       ENDIF
    end do time_integration

    CALL dynredem1("restart.nc", vcov, ucov, teta, q, masse, ps, &
         itau=itau_dyn+itaufin)

    ! Calcul des tendances dynamiques:
    CALL geopot(ip1jmp1, teta, pk, pks, phis, phi)
    CALL caldyn(itaufin, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
         MOD(itaufin, iconser) == 0, du, dv, dteta, dp, w, pbaru, pbarv, &
         time_0)

  END SUBROUTINE leapfrog

end module leapfrog_m
