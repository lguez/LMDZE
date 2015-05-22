module leapfrog_m

  IMPLICIT NONE

contains

  SUBROUTINE leapfrog(ucov, vcov, teta, ps, masse, phis, q)

    ! From dyn3d/leapfrog.F, version 1.6, 2005/04/13 08:58:34 revision 616
    ! Authors: P. Le Van, L. Fairhead, F. Hourdin

    ! Intégration temporelle du modèle : Matsuno-leapfrog scheme.

    use addfi_m, only: addfi
    use bilan_dyn_m, only: bilan_dyn
    use caladvtrac_m, only: caladvtrac
    use caldyn_m, only: caldyn
    USE calfis_m, ONLY: calfis
    USE comconst, ONLY: daysec, dtvr
    USE comgeom, ONLY: aire_2d, apoln, apols
    USE disvert_m, ONLY: ap, bp
    USE conf_gcm_m, ONLY: day_step, iconser, iperiod, iphysiq, nday, offline, &
         iflag_phys, iecri
    USE conf_guide_m, ONLY: ok_guide
    USE dimens_m, ONLY: iim, jjm, llm, nqmx
    use dissip_m, only: dissip
    USE dynetat0_m, ONLY: day_ini
    use dynredem1_m, only: dynredem1
    USE exner_hyb_m, ONLY: exner_hyb
    use filtreg_scal_m, only: filtreg_scal
    use fluxstokenc_m, only: fluxstokenc
    use geopot_m, only: geopot
    USE guide_m, ONLY: guide
    use inidissip_m, only: idissip
    use integrd_m, only: integrd
    use nr_util, only: assert
    USE pressure_var, ONLY: p3d
    USE temps, ONLY: itau_dyn
    use writedynav_m, only: writedynav
    use writehist_m, only: writehist

    ! Variables dynamiques:
    REAL, intent(inout):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, intent(inout):: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant

    REAL, intent(inout):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! potential temperature

    REAL, intent(inout):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol, en Pa
    REAL, intent(inout):: masse(:, :, :) ! (iim + 1, jjm + 1, llm) masse d'air
    REAL, intent(in):: phis(:, :) ! (iim + 1, jjm + 1) surface geopotential

    REAL, intent(inout):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    ! mass fractions of advected fields

    ! Local:

    ! Variables dynamiques:

    REAL pks(iim + 1, jjm + 1) ! exner au sol
    REAL pk(iim + 1, jjm + 1, llm) ! exner au milieu des couches
    REAL pkf(iim + 1, jjm + 1, llm) ! exner filtr\'e au milieu des couches
    REAL phi(iim + 1, jjm + 1, llm) ! geopotential
    REAL w(iim + 1, jjm + 1, llm) ! vitesse verticale

    ! Variables dynamiques intermediaire pour le transport
    ! Flux de masse :
    REAL pbaru(iim + 1, jjm + 1, llm), pbarv(iim + 1, jjm, llm)

    ! Variables dynamiques au pas - 1
    REAL vcovm1(iim + 1, jjm, llm), ucovm1(iim + 1, jjm + 1, llm)
    REAL tetam1(iim + 1, jjm + 1, llm), psm1(iim + 1, jjm + 1)
    REAL massem1(iim + 1, jjm + 1, llm)

    ! Tendances dynamiques
    REAL dv((iim + 1) * jjm, llm), dudyn(iim + 1, jjm + 1, llm)
    REAL dteta(iim + 1, jjm + 1, llm)
    real dp((iim + 1) * (jjm + 1))

    ! Tendances de la dissipation :
    REAL dvdis(iim + 1, jjm, llm), dudis(iim + 1, jjm + 1, llm)
    REAL dtetadis(iim + 1, jjm + 1, llm)

    ! Tendances physiques
    REAL dvfi(iim + 1, jjm, llm), dufi(iim + 1, jjm + 1, llm)
    REAL dtetafi(iim + 1, jjm + 1, llm), dqfi(iim + 1, jjm + 1, llm, nqmx)

    ! Variables pour le fichier histoire
    INTEGER itau ! index of the time step of the dynamics, starts at 0
    INTEGER itaufin
    real finvmaold(iim + 1, jjm + 1, llm)
    INTEGER l

    ! Variables test conservation \'energie
    REAL ecin(iim + 1, jjm + 1, llm), ecin0(iim + 1, jjm + 1, llm)

    REAL vcont((iim + 1) * jjm, llm), ucont((iim + 1) * (jjm + 1), llm)
    logical leapf
    real dt ! time step, in s

    !---------------------------------------------------

    print *, "Call sequence information: leapfrog"
    call assert(shape(ucov) == (/iim + 1, jjm + 1, llm/), "leapfrog")

    itaufin = nday * day_step
    ! "day_step" is a multiple of "iperiod", therefore so is "itaufin".

    ! On initialise la pression et la fonction d'Exner :
    forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
    CALL exner_hyb(ps, p3d, pks, pk)
    pkf = pk
    CALL filtreg_scal(pkf, direct = .true., intensive = .true.)

    time_integration: do itau = 0, itaufin - 1
       leapf = mod(itau, iperiod) /= 0
       if (leapf) then
          dt = 2 * dtvr
       else
          ! Matsuno
          dt = dtvr
          if (ok_guide) call guide(itau, ucov, vcov, teta, q(:, :, :, 1), ps)
          vcovm1 = vcov
          ucovm1 = ucov
          tetam1 = teta
          massem1 = masse
          psm1 = ps
          finvmaold = masse
          CALL filtreg_scal(finvmaold, direct = .false., intensive = .false.)
       end if

       ! Calcul des tendances dynamiques:
       CALL geopot(teta, pk, pks, phis, phi)
       CALL caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
            dudyn, dv, dteta, dp, w, pbaru, pbarv, &
            conser = MOD(itau, iconser) == 0)

       CALL caladvtrac(q, pbaru, pbarv, p3d, masse, teta, pk)

       ! Stokage du flux de masse pour traceurs offline:
       IF (offline) CALL fluxstokenc(pbaru, pbarv, masse, teta, phi, phis, &
            dtvr, itau)

       ! Int\'egrations dynamique et traceurs:
       CALL integrd(vcovm1, ucovm1, tetam1, psm1, massem1, dv, dudyn, dteta, &
            dp, vcov, ucov, teta, q(:, :, :, :2), ps, masse, finvmaold, dt, &
            leapf)

       forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
       CALL exner_hyb(ps, p3d, pks, pk)
       pkf = pk
       CALL filtreg_scal(pkf, direct = .true., intensive = .true.)

       if (.not. leapf) then
          ! Matsuno backward
          ! Calcul des tendances dynamiques:
          CALL geopot(teta, pk, pks, phis, phi)
          CALL caldyn(itau + 1, ucov, vcov, teta, ps, masse, pk, pkf, phis, &
               phi, dudyn, dv, dteta, dp, w, pbaru, pbarv, conser = .false.)

          ! integrations dynamique et traceurs:
          CALL integrd(vcovm1, ucovm1, tetam1, psm1, massem1, dv, dudyn, &
               dteta, dp, vcov, ucov, teta, q(:, :, :, :2), ps, masse, &
               finvmaold, dtvr, leapf=.false.)

          forall (l = 1: llm + 1) p3d(:, :, l) = ap(l) + bp(l) * ps
          CALL exner_hyb(ps, p3d, pks, pk)
          pkf = pk
          CALL filtreg_scal(pkf, direct = .true., intensive = .true.)
       end if

       IF (MOD(itau + 1, iphysiq) == 0 .AND. iflag_phys /= 0) THEN
          CALL calfis(itau / day_step + day_ini, &
               REAL(mod(itau, day_step)) / day_step, ucov, vcov, teta, q, pk, &
               phis, phi, w, dufi, dvfi, dtetafi, dqfi, &
               lafin = itau + 1 == itaufin)

          CALL addfi(ucov, vcov, teta, q, dufi, dvfi, dtetafi, dqfi)
       ENDIF

       IF (MOD(itau + 1, idissip) == 0) THEN
          ! Dissipation horizontale et verticale des petites \'echelles

          ! calcul de l'\'energie cin\'etique avant dissipation
          call covcont(llm, ucov, vcov, ucont, vcont)
          call enercin(vcov, ucov, vcont, ucont, ecin0)

          ! dissipation
          CALL dissip(vcov, ucov, teta, p3d, dvdis, dudis, dtetadis)
          ucov = ucov + dudis
          vcov = vcov + dvdis

          ! On ajoute la tendance due \`a la transformation \'energie
          ! cin\'etique en \'energie thermique par la dissipation
          call covcont(llm, ucov, vcov, ucont, vcont)
          call enercin(vcov, ucov, vcont, ucont, ecin)
          dtetadis = dtetadis + (ecin0 - ecin) / pk
          teta = teta + dtetadis

          ! Calcul de la valeur moyenne aux p\^oles :
          forall (l = 1: llm)
             teta(:, 1, l) = SUM(aire_2d(:iim, 1) * teta(:iim, 1, l)) &
                  / apoln
             teta(:, jjm + 1, l) = SUM(aire_2d(:iim, jjm+1) &
                  * teta(:iim, jjm + 1, l)) / apols
          END forall
       END IF

       IF (MOD(itau + 1, iperiod) == 0) THEN
          ! \'Ecriture du fichier histoire moyenne:
          CALL writedynav(vcov, ucov, teta, pk, phi, q, masse, ps, phis, &
               time = itau + 1)
          call bilan_dyn(ps, masse, pk, pbaru, pbarv, teta, phi, ucov, vcov, &
               q(:, :, :, 1))
       ENDIF

       IF (MOD(itau + 1, iecri * day_step) == 0) THEN
          CALL geopot(teta, pk, pks, phis, phi)
          CALL writehist(itau, vcov, ucov, teta, phi, masse, ps)
       END IF
    end do time_integration

    CALL dynredem1("restart.nc", vcov, ucov, teta, q, masse, ps, &
         itau = itau_dyn + itaufin)

    ! Calcul des tendances dynamiques:
    CALL geopot(teta, pk, pks, phis, phi)
    CALL caldyn(itaufin, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
         dudyn, dv, dteta, dp, w, pbaru, pbarv, &
         conser = MOD(itaufin, iconser) == 0)

  END SUBROUTINE leapfrog

end module leapfrog_m
