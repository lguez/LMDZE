module leapfrog_m

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT NONE

contains

  SUBROUTINE leapfrog(ucov, vcov, teta, ps, masse, phis, nq, q, clesphy0, &
       time_0)

    ! From dyn3d/leapfrog.F, version 1.6 2005/04/13 08:58:34

    ! Version du 10/01/98, avec coordonnees verticales hybrides, avec
    ! nouveaux operat. dissipation * (gradiv2, divgrad2, nxgraro2)

    ! Auteur: P. Le Van /L. Fairhead/F.Hourdin
    ! Objet:
    ! GCM LMD nouvelle grille

    ! ... Dans inigeom, nouveaux calculs pour les elongations cu, cv
    ! et possibilite d'appeler une fonction f(y) a derivee tangente
    ! hyperbolique a la place de la fonction a derivee sinusoidale.

    ! ... Possibilite de choisir le shema pour l'advection de
    ! q, en modifiant iadv dans "traceur.def" (10/02) .

    ! Pour Van-Leer + Vapeur d'eau saturee, iadv(1)=4. (F.Codron, 10/99)
    ! Pour Van-Leer iadv=10 

    use dimens_m, only: iim, llm, nqmx
    use paramet_m, only: ip1jmp1, ip1jm, llmp1, ijmllm, ijp1llm, jjp1, iip1, &
         iip2
    use comconst, only: dtvr, daysec, dtphys
    use comvert, only: ap, bp
    use conf_gcm_m, only: day_step, iconser, idissip, iphysiq, iperiod, nday, &
         offline, periodav
    use logic, only: ok_guide, apdiss, apphys, conser, forward, iflag_phys, &
         leapf, statcl
    use comgeom
    use serre
    use temps, only: itaufin, day_ini, dt
    use iniprint, only: prt_level
    use com_io_dyn
    use abort_gcm_m, only: abort_gcm
    use ener
    use calfis_m, only: calfis
    use exner_hyb_m, only: exner_hyb
    use guide_m, only: guide
    use pression_m, only: pression

    integer nq

    INTEGER longcles
    PARAMETER (longcles = 20)
    REAL clesphy0(longcles)

    ! variables dynamiques
    REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm) ! vents covariants
    REAL teta(ip1jmp1, llm) ! temperature potentielle 
    REAL q(ip1jmp1, llm, nqmx) ! mass fractions of advected fields
    REAL ps(ip1jmp1) ! pression au sol
    REAL p(ip1jmp1, llmp1) ! pression aux interfac.des couches
    REAL pks(ip1jmp1) ! exner au sol
    REAL pk(ip1jmp1, llm) ! exner au milieu des couches
    REAL pkf(ip1jmp1, llm) ! exner filt.au milieu des couches
    REAL masse(ip1jmp1, llm) ! masse d'air
    REAL phis(ip1jmp1) ! geopotentiel au sol
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

    INTEGER itau, itaufinp1
    INTEGER iday ! jour julien
    REAL time ! Heure de la journee en fraction d'1 jour

    REAL SSUM
    REAL time_0, finvmaold(ip1jmp1, llm)

    LOGICAL :: lafin=.false.
    INTEGER ij, l

    REAL rdayvrai, rdaym_ini
    LOGICAL callinigrads

    data callinigrads/.true./

    !+jld variables test conservation energie
    REAL ecin(ip1jmp1, llm), ecin0(ip1jmp1, llm)
    ! Tendance de la temp. potentiel d (theta) / d t due a la 
    ! tansformation d'energie cinetique en energie thermique
    ! cree par la dissipation
    REAL dtetaecdt(ip1jmp1, llm)
    REAL vcont(ip1jm, llm), ucont(ip1jmp1, llm)
    CHARACTER*15 ztit
    INTEGER ip_ebil_dyn ! PRINT level for energy conserv. diag.
    SAVE ip_ebil_dyn
    DATA ip_ebil_dyn /0/

    character(len=*), parameter:: modname = "leapfrog"
    character*80 abort_message

    logical dissip_conservative
    save dissip_conservative
    data dissip_conservative /.true./

    LOGICAL prem
    save prem
    DATA prem /.true./

    !---------------------------------------------------

    print *, "Call sequence information: leapfrog"

    itaufin = nday * day_step
    itaufinp1 = itaufin + 1

    itau = 0
    iday = day_ini
    time = time_0
    IF (time > 1.) THEN
       time = time - 1.
       iday = iday + 1
    ENDIF

    ! On initialise la pression et la fonction d'Exner :
    dq=0.
    CALL pression(ip1jmp1, ap, bp, ps, p)
    CALL exner_hyb(ps, p, pks, pk, pkf)

    ! Debut de l'integration temporelle:
    do
       if (ok_guide.and.(itaufin - itau - 1) * dtvr > 21600) then
          call guide(itau, ucov, vcov, teta, q, masse, ps)
       else
          IF (prt_level > 9) print *, &
               'Attention : on ne guide pas les 6 dernieres heures.'
       endif

       CALL SCOPY(ijmllm, vcov, 1, vcovm1, 1)
       CALL SCOPY(ijp1llm, ucov, 1, ucovm1, 1)
       CALL SCOPY(ijp1llm, teta, 1, tetam1, 1)
       CALL SCOPY(ijp1llm, masse, 1, massem1, 1)
       CALL SCOPY(ip1jmp1, ps, 1, psm1, 1)

       forward = .TRUE.
       leapf = .FALSE.
       dt = dtvr

       CALL SCOPY(ijp1llm, masse, 1, finvmaold, 1)
       CALL filtreg(finvmaold, jjp1, llm, - 2, 2, .TRUE., 1)

       do
          ! gestion des appels de la physique et des dissipations:

          apphys = .FALSE.
          statcl = .FALSE.
          conser = .FALSE.
          apdiss = .FALSE.

          IF (MOD(itau, iconser) == 0) conser = .TRUE.
          IF (MOD(itau + 1, idissip) == 0) apdiss = .TRUE.
          IF (MOD(itau + 1, iphysiq) == 0 .AND. iflag_phys /= 0) apphys=.TRUE.

          ! calcul des tendances dynamiques:

          CALL geopot(ip1jmp1, teta, pk, pks, phis, phi)

          CALL caldyn(itau, ucov, vcov, teta, ps, masse, pk, pkf, phis, phi, &
               conser, du, dv, dteta, dp, w, pbaru, pbarv, &
               time + iday - day_ini)

          ! calcul des tendances advection des traceurs (dont l'humidite)

          IF (forward .OR. leapf) THEN
             CALL caladvtrac(q, pbaru, pbarv, p, masse, dq, teta, pk)
             IF (offline) THEN
                !maf stokage du flux de masse pour traceurs OFF-LINE
                CALL fluxstokenc(pbaru, pbarv, masse, teta, phi, phis, dtvr, &
                     itau)
             ENDIF
          ENDIF

          ! integrations dynamique et traceurs:
          CALL integrd(2, vcovm1, ucovm1, tetam1, psm1, massem1, dv, du, &
               dteta, dq, dp, vcov, ucov, teta, q, ps, masse, phis, finvmaold)

          ! calcul des tendances physiques:

          IF (apphys) THEN
             IF (itau + 1 == itaufin) lafin = .TRUE.

             CALL pression(ip1jmp1, ap, bp, ps, p)
             CALL exner_hyb(ps, p, pks, pk, pkf)

             rdaym_ini = itau * dtvr / daysec
             rdayvrai = rdaym_ini + day_ini

             ! Interface avec les routines de phylmd (phymars ...)

             ! Diagnostique de conservation de l'énergie : initialisation
             IF (ip_ebil_dyn >= 1) THEN 
                ztit='bil dyn'
                CALL diagedyn(ztit, 2, 1, 1, dtphys &
                     , ucov, vcov, ps, p, pk, teta, q(:, :, 1), q(:, :, 2))
             ENDIF

             CALL calfis(nq, lafin, rdayvrai, time, ucov, vcov, teta, q, &
                  masse, ps, p, pk, phis, phi, du, dv, dteta, dq, w, &
                  clesphy0, dufi, dvfi, dtetafi, dqfi, dpfi)

             ! ajout des tendances physiques:
             CALL addfi(nqmx, dtphys, &
                  ucov, vcov, teta, q, ps, &
                  dufi, dvfi, dtetafi, dqfi, dpfi)

             ! Diagnostique de conservation de l'énergie : difference
             IF (ip_ebil_dyn >= 1) THEN 
                ztit = 'bil phys'
                CALL diagedyn(ztit, 2, 1, 1, dtphys, ucov, vcov, ps, p, pk, &
                     teta, q(:, :, 1), q(:, :, 2))
             ENDIF
          ENDIF

          CALL pression(ip1jmp1, ap, bp, ps, p)
          CALL exner_hyb(ps, p, pks, pk, pkf)

          ! dissipation horizontale et verticale des petites echelles:

          IF (apdiss) THEN
             ! calcul de l'energie cinetique avant dissipation
             call covcont(llm, ucov, vcov, ucont, vcont)
             call enercin(vcov, ucov, vcont, ucont, ecin0)

             ! dissipation
             CALL dissip(vcov, ucov, teta, p, dvdis, dudis, dtetadis)
             ucov=ucov + dudis
             vcov=vcov + dvdis

             if (dissip_conservative) then
                ! On rajoute la tendance due a la transform. Ec -> E
                ! therm. cree lors de la dissipation
                call covcont(llm, ucov, vcov, ucont, vcont)
                call enercin(vcov, ucov, vcont, ucont, ecin)
                dtetaecdt= (ecin0 - ecin) / pk
                dtetadis=dtetadis + dtetaecdt
             endif
             teta=teta + dtetadis

             ! Calcul de la valeur moyenne, unique de h aux poles .....

             DO l = 1, llm
                DO ij = 1, iim
                   tppn(ij) = aire(ij) * teta(ij, l)
                   tpps(ij) = aire(ij + ip1jm) * teta(ij + ip1jm, l)
                ENDDO
                tpn = SSUM(iim, tppn, 1) / apoln
                tps = SSUM(iim, tpps, 1) / apols

                DO ij = 1, iip1
                   teta(ij, l) = tpn
                   teta(ij + ip1jm, l) = tps
                ENDDO
             ENDDO

             DO ij = 1, iim
                tppn(ij) = aire(ij) * ps(ij)
                tpps(ij) = aire(ij + ip1jm) * ps(ij + ip1jm)
             ENDDO
             tpn = SSUM(iim, tppn, 1) / apoln
             tps = SSUM(iim, tpps, 1) / apols

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

          IF (itau == itaufinp1) then 
             abort_message = 'Simulation finished'
             call abort_gcm(modname, abort_message, 0)
          ENDIF

          ! ecriture du fichier histoire moyenne:

          ! Comment out the following calls when you do not want the output
          ! files "dyn_hist_ave.nc" and "dynzon.nc"
          IF (MOD(itau, iperiod) == 0 .OR. itau == itaufin) THEN
             CALL writedynav(histaveid, nqmx, itau, vcov, &
                  ucov, teta, pk, phi, q, masse, ps, phis)
             call bilan_dyn(2, dtvr * iperiod, dtvr * day_step * periodav, &
                  ps, masse, pk, pbaru, pbarv, teta, phi, ucov, vcov, q)
          ENDIF

          IF (itau == itaufin) THEN
             CALL dynredem1("restart.nc", 0., vcov, ucov, teta, q, masse, ps)
             CLOSE(99)
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
             ! ...... pas leapfrog .....
             leapf = .TRUE.
             dt = 2. * dtvr
          END IF
       end do
    end do

  END SUBROUTINE leapfrog

end module leapfrog_m
