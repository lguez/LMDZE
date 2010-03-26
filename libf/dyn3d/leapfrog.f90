module leapfrog_m

  IMPLICIT NONE

contains

  SUBROUTINE leapfrog(ucov, vcov, teta, ps, masse, phis, q, time_0)

    ! From dyn3d/leapfrog.F, version 1.6, 2005/04/13 08:58:34
    ! Authors: P. Le Van, L. Fairhead, F. Hourdin

    USE calfis_m, ONLY: calfis
    USE com_io_dyn, ONLY: histaveid
    USE comconst, ONLY: daysec, dtphys, dtvr
    USE comgeom, ONLY: aire, apoln, apols
    USE comvert, ONLY: ap, bp
    USE conf_gcm_m, ONLY: day_step, iconser, iperiod, iphysiq, nday, offline, &
         periodav
    USE dimens_m, ONLY: iim, llm, nqmx
    USE dynetat0_m, ONLY: day_ini
    use dynredem1_m, only: dynredem1
    USE exner_hyb_m, ONLY: exner_hyb
    use filtreg_m, only: filtreg
    USE guide_m, ONLY: guide
    use inidissip_m, only: idissip
    USE logic, ONLY: iflag_phys, ok_guide
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1, jjp1
    USE pression_m, ONLY: pression
    USE pressure_var, ONLY: p3d
    USE temps, ONLY: itau_dyn

    ! Variables dynamiques:
    REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm) ! vents covariants
    REAL teta(ip1jmp1, llm) ! temperature potentielle 
    REAL ps(ip1jmp1) ! pression au sol, en Pa

    REAL masse(ip1jmp1, llm) ! masse d'air
    REAL phis(ip1jmp1) ! geopotentiel au sol
    REAL q(ip1jmp1, llm, nqmx) ! mass fractions of advected fields
    REAL, intent(in):: time_0

    ! Variables local to the procedure:

    ! Variables dynamiques:

    REAL pks(ip1jmp1) ! exner au sol
    REAL pk(ip1jmp1, llm) ! exner au milieu des couches
    REAL pkf(ip1jmp1, llm) ! exner filt.au milieu des couches
    REAL phi(ip1jmp1, llm) ! geopotential
    REAL w(ip1jmp1, llm) ! vitesse verticale

    ! variables dynamiques intermediaire pour le transport
    REAL pbaru(ip1jmp1, llm), pbarv(ip1jm, llm) !flux de masse

    ! variables dynamiques au pas - 1
    REAL vcovm1(ip1jm, llm), ucovm1(ip1jmp1, llm)
    REAL tetam1(ip1jmp1, llm), psm1(ip1jmp1)
    REAL massem1(ip1jmp1, llm)

    ! tendances dynamiques
    REAL dv(ip1jm, llm), du(ip1jmp1, llm)
    REAL dteta(ip1jmp1, llm), dq(ip1jmp1, llm, nqmx), dp(ip1jmp1)

    ! tendances de la dissipation
    REAL dvdis(ip1jm, llm), dudis(ip1jmp1, llm)
    REAL dtetadis(ip1jmp1, llm)

    ! tendances physiques
    REAL dvfi(ip1jm, llm), dufi(ip1jmp1, llm)
    REAL dtetafi(ip1jmp1, llm), dqfi(ip1jmp1, llm, nqmx), dpfi(ip1jmp1)

    ! variables pour le fichier histoire

    REAL tppn(iim), tpps(iim), tpn, tps

    INTEGER itau ! index of the time step of the dynamics, starts at 0
    INTEGER itaufin
    INTEGER iday ! jour julien
    REAL time ! time of day, as a fraction of day length
    real finvmaold(ip1jmp1, llm)
    LOGICAL:: lafin=.false.
    INTEGER ij, l

    REAL rdayvrai, rdaym_ini

    ! Variables test conservation energie
    REAL ecin(ip1jmp1, llm), ecin0(ip1jmp1, llm)
    ! Tendance de la temp. potentiel d (theta) / d t due a la 
    ! tansformation d'energie cinetique en energie thermique
    ! cree par la dissipation
    REAL dtetaecdt(ip1jmp1, llm)
    REAL vcont(ip1jm, llm), ucont(ip1jmp1, llm)
    logical forward, leapf
    REAL dt

    !---------------------------------------------------

    print *, "Call sequence information: leapfrog"

    itaufin = nday * day_step
    itau = 0
    iday = day_ini
    time = time_0
    dq = 0.
    ! On initialise la pression et la fonction d'Exner :
    CALL pression(ip1jmp1, ap, bp, ps, p3d)
    CALL exner_hyb(ps, p3d, pks, pk, pkf)

    ! Début de l'integration temporelle :
    outer_loop:do
       if (ok_guide .and. (itaufin - itau - 1) * dtvr > 21600.) &
            call guide(itau, ucov, vcov, teta, q, masse, ps)
       vcovm1 = vcov
       ucovm1 = ucov
       tetam1 = teta
       massem1 = masse
       psm1 = ps
       forward = .TRUE.
       leapf = .FALSE.
       dt = dtvr
       finvmaold = masse
       CALL filtreg(finvmaold, jjp1, llm, - 2, 2, .TRUE., 1)

       do
          ! Calcul des tendances dynamiques:
          CALL geopot(ip1jmp1, teta, pk, pks, phis, phi)
          CALL caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
               MOD(itau, iconser) == 0, du, dv, dteta, dp, w, pbaru, pbarv, &
               time + iday - day_ini)

          IF (forward .OR. leapf) THEN
             ! Calcul des tendances advection des traceurs (dont l'humidité)
             CALL caladvtrac(q, pbaru, pbarv, p3d, masse, dq, teta, pk)
             IF (offline) THEN
                ! Stokage du flux de masse pour traceurs off-line
                CALL fluxstokenc(pbaru, pbarv, masse, teta, phi, phis, dtvr, &
                     itau)
             ENDIF
          ENDIF

          ! integrations dynamique et traceurs:
          CALL integrd(2, vcovm1, ucovm1, tetam1, psm1, massem1, dv, du, &
               dteta, dq, dp, vcov, ucov, teta, q, ps, masse, phis, &
               finvmaold, leapf, dt)

          IF (MOD(itau + 1, iphysiq) == 0 .AND. iflag_phys /= 0) THEN
             ! calcul des tendances physiques:
             IF (itau + 1 == itaufin) lafin = .TRUE.

             CALL pression(ip1jmp1, ap, bp, ps, p3d)
             CALL exner_hyb(ps, p3d, pks, pk, pkf)

             rdaym_ini = itau * dtvr / daysec
             rdayvrai = rdaym_ini + day_ini

             CALL calfis(nqmx, lafin, rdayvrai, time, ucov, vcov, teta, q, &
                  masse, ps, pk, phis, phi, du, dv, dteta, dq, w, &
                  dufi, dvfi, dtetafi, dqfi, dpfi)

             ! ajout des tendances physiques:
             CALL addfi(nqmx, dtphys, &
                  ucov, vcov, teta, q, ps, &
                  dufi, dvfi, dtetafi, dqfi, dpfi)
          ENDIF

          CALL pression(ip1jmp1, ap, bp, ps, p3d)
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

             ! Calcul de la valeur moyenne unique de h aux pôles
             DO l = 1, llm
                DO ij = 1, iim
                   tppn(ij) = aire(ij) * teta(ij, l)
                   tpps(ij) = aire(ij + ip1jm) * teta(ij + ip1jm, l)
                ENDDO
                tpn = SUM(tppn) / apoln
                tps = SUM(tpps) / apols

                DO ij = 1, iip1
                   teta(ij, l) = tpn
                   teta(ij + ip1jm, l) = tps
                ENDDO
             ENDDO

             DO ij = 1, iim
                tppn(ij) = aire(ij) * ps(ij)
                tpps(ij) = aire(ij + ip1jm) * ps(ij + ip1jm)
             ENDDO
             tpn = SUM(tppn) / apoln
             tps = SUM(tpps) / apols

             DO ij = 1, iip1
                ps(ij) = tpn
                ps(ij + ip1jm) = tps
             ENDDO
          END IF

          ! fin de l'intégration dynamique et physique pour le pas "itau"
          ! préparation du pas d'intégration suivant

          ! schema matsuno + leapfrog
          IF (forward .OR. leapf) THEN
             itau = itau + 1
             iday = day_ini + itau / day_step
             time = REAL(itau - (iday - day_ini) * day_step) / day_step &
                  + time_0
             IF (time > 1.) THEN
                time = time - 1.
                iday = iday + 1
             ENDIF
          ENDIF

          IF (itau == itaufin + 1) exit outer_loop

          IF (MOD(itau, iperiod) == 0 .OR. itau == itaufin) THEN
             ! ecriture du fichier histoire moyenne:
             CALL writedynav(histaveid, nqmx, itau, vcov, &
                  ucov, teta, pk, phi, q, masse, ps, phis)
             call bilan_dyn(2, dtvr * iperiod, dtvr * day_step * periodav, &
                  ps, masse, pk, pbaru, pbarv, teta, phi, ucov, vcov, q)
          ENDIF

          IF (itau == itaufin) THEN
             CALL dynredem1("restart.nc", vcov, ucov, teta, q, masse, ps, &
                  itau=itau_dyn+itaufin)
          ENDIF

          ! gestion de l'integration temporelle:
          IF (MOD(itau, iperiod) == 0) exit
          IF (MOD(itau - 1, iperiod) == 0) THEN
             IF (forward) THEN
                ! fin du pas forward et debut du pas backward
                forward = .FALSE.
                leapf = .FALSE.
             ELSE
                ! fin du pas backward et debut du premier pas leapfrog
                leapf = .TRUE.
                dt = 2. * dtvr
             END IF
          ELSE
             ! pas leapfrog
             leapf = .TRUE.
             dt = 2. * dtvr
          END IF
       end do
    end do outer_loop

  END SUBROUTINE leapfrog

end module leapfrog_m
