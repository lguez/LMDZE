module interfoce_slab_m

  implicit none

contains

  SUBROUTINE interfoce_slab(klon, debut, itap, dtime, ijour,  &
       radsol, fluxo, fluxg, pctsrf,  &
       tslab, seaice, pctsrf_slab)

    ! Cette routine calcule la temperature d'un slab ocean, la glace de mer
    ! et les pourcentages de la maille couverte par l'ocean libre et/ou
    ! la glace de mer pour un "slab" ocean de 50m

    ! I. Musat 04.02.2005

    ! input:
    ! klon nombre total de points de grille
    ! debut logical: 1er appel a la physique
    ! itap numero du pas de temps
    ! dtime pas de temps de la physique (en s)
    ! ijour jour dans l'annee en cours
    ! radsol rayonnement net au sol (LW + SW)
    ! fluxo flux turbulent (sensible + latent) sur les mailles oceaniques
    ! fluxg flux de conduction entre la surface de la glace de mer et l'ocean
    ! pctsrf tableau des pourcentages de surface de chaque maille
    ! output:
    ! tslab temperature de l'ocean libre
    ! seaice glace de mer (kg/m2)
    ! pctsrf_slab "pourcentages" (valeurs entre 0. et 1.) surfaces issus du slab

    use indicesol
    use clesphys
    use abort_gcm_m, only: abort_gcm
    use SUPHEC_M

    ! Parametres d'entree
    integer, intent(IN) :: klon
    logical, intent(IN) :: debut
    INTEGER, intent(IN) :: itap
    REAL, intent(IN) :: dtime
    INTEGER, intent(IN) :: ijour
    REAL, dimension(klon), intent(IN) :: radsol
    REAL, dimension(klon), intent(IN) :: fluxo
    REAL, dimension(klon), intent(IN) :: fluxg
    real, dimension(klon, nbsrf), intent(IN) :: pctsrf
    ! Parametres de sortie
    real, dimension(klon), intent(INOUT) :: tslab
    real, dimension(klon), intent(INOUT) :: seaice ! glace de mer (kg/m2)
    real, dimension(klon, nbsrf), intent(OUT) :: pctsrf_slab

    ! Variables locales :
    INTEGER, save :: lmt_pas, julien, idayvrai
    REAL, parameter :: unjour=86400.
    real, allocatable, dimension(:), save :: tmp_tslab, tmp_seaice
    REAL, allocatable, dimension(:), save :: slab_bils
    REAL, allocatable, dimension(:), save :: lmt_bils
    logical, save :: check = .false.

    REAL, parameter :: cyang=50.0 * 4.228e+06 ! capacite calorifique volumetrique de l'eau J/(m2 K)
    REAL, parameter :: cbing=0.334e+05 ! J/kg
    real, dimension(klon) :: siceh !hauteur de la glace de mer (m)
    INTEGER :: i
    integer :: sum_error, error
    REAL :: zz, za, zb

    character (len = 80) :: abort_message
    character (len = 20) :: modname = 'interfoce_slab'

    julien = MOD(ijour, 360)
    sum_error = 0
    IF (debut) THEN
       allocate(slab_bils(klon), stat = error)
       sum_error = sum_error + error
       allocate(lmt_bils(klon), stat = error)
       sum_error = sum_error + error
       allocate(tmp_tslab(klon), stat = error)
       sum_error = sum_error + error
       allocate(tmp_seaice(klon), stat = error)
       sum_error = sum_error + error
       if (sum_error /= 0) then
          abort_message='Pb allocation var. slab_bils, lmt_bils, tmp_tslab, tmp_seaice'
          call abort_gcm(modname, abort_message, 1)
       endif
       tmp_tslab=tslab
       tmp_seaice=seaice
       lmt_pas = nint(86400./dtime * 1.0) ! pour une lecture une fois par jour

       IF (check) THEN
          PRINT*, 'interfoce_slab klon, debut, itap, dtime, ijour, &
               & lmt_pas ', klon, debut, itap, dtime, ijour, &
               lmt_pas
       ENDIF !check

       PRINT*, '************************'
       PRINT*, 'SLAB OCEAN est actif, prenez precautions !'
       PRINT*, '************************'

       ! a mettre un slab_bils aussi en force !!!

       DO i = 1, klon
          slab_bils(i) = 0.0
       ENDDO

    ENDIF !debut
    pctsrf_slab(1:klon, 1:nbsrf) = pctsrf(1:klon, 1:nbsrf)

    ! lecture du bilan au sol lmt_bils issu d'une simulation forcee en debut de journee

    IF (MOD(itap, lmt_pas) .EQ. 1) THEN !1er pas de temps de la journee
       idayvrai = ijour
       CALL condsurf(julien, idayvrai, lmt_bils)
    ENDIF !(MOD(itap-1, lmt_pas) .EQ. 0) THEN

    DO i = 1, klon
       IF((pctsrf_slab(i, is_oce).GT.epsfra).OR. &
            (pctsrf_slab(i, is_sic).GT.epsfra)) THEN

          ! fabriquer de la glace si congelation atteinte:

          IF (tmp_tslab(i).LT.(RTT-1.8)) THEN
             zz = (RTT-1.8)-tmp_tslab(i)
             tmp_seaice(i) = tmp_seaice(i) + cyang/cbing * zz
             seaice(i) = tmp_seaice(i)
             tmp_tslab(i) = RTT-1.8
          ENDIF

          ! faire fondre de la glace si temperature est superieure a 0:

          IF ((tmp_tslab(i).GT.RTT) .AND. (tmp_seaice(i).GT.0.0)) THEN
             zz = cyang/cbing * (tmp_tslab(i)-RTT)
             zz = MIN(zz, tmp_seaice(i))
             tmp_seaice(i) = tmp_seaice(i) - zz
             seaice(i) = tmp_seaice(i)
             tmp_tslab(i) = tmp_tslab(i) - zz*cbing/cyang
          ENDIF

          ! limiter la glace de mer a 10 metres (10000 kg/m2)

          IF(tmp_seaice(i).GT.45.) THEN
             tmp_seaice(i) = MIN(tmp_seaice(i), 10000.0)
          ELSE
             tmp_seaice(i) = 0.
          ENDIF
          seaice(i) = tmp_seaice(i)
          siceh(i)=tmp_seaice(i)/1000. !en metres

          ! determiner la nature du sol (glace de mer ou ocean libre):

          ! on fait dependre la fraction de seaice "pctsrf(i, is_sic)"
          ! de l'epaisseur de seaice :
          ! pctsrf(i, is_sic)=1. si l'epaisseur de la glace de mer est >= a 20cm
          ! et pctsrf(i, is_sic) croit lineairement avec seaice de 0. a 20cm d'epaisseur

          pctsrf_slab(i, is_sic)=MIN(siceh(i)/0.20, &
               1.-(pctsrf_slab(i, is_ter)+pctsrf_slab(i, is_lic)))
          pctsrf_slab(i, is_oce)=1.0 - &
               (pctsrf_slab(i, is_ter)+pctsrf_slab(i, is_lic)+pctsrf_slab(i, is_sic))
       ENDIF !pctsrf
    ENDDO

    ! Calculer le bilan du flux de chaleur au sol :

    DO i = 1, klon
       za = radsol(i) + fluxo(i)
       zb = fluxg(i)
       IF((pctsrf_slab(i, is_oce).GT.epsfra).OR. &
            (pctsrf_slab(i, is_sic).GT.epsfra)) THEN
          slab_bils(i)=slab_bils(i)+(za*pctsrf_slab(i, is_oce) &
               +zb*pctsrf_slab(i, is_sic))/ FLOAT(lmt_pas)
       ENDIF
    ENDDO !klon

    ! calcul tslab

    IF (MOD(itap, lmt_pas).EQ.0) THEN !fin de journee
       DO i = 1, klon
          IF ((pctsrf_slab(i, is_oce).GT.epsfra).OR. &
               (pctsrf_slab(i, is_sic).GT.epsfra)) THEN
             tmp_tslab(i) = tmp_tslab(i) + &
                  (slab_bils(i)-lmt_bils(i)) &
                  /cyang*unjour
             ! on remet l'accumulation a 0
             slab_bils(i) = 0.
          ENDIF !pctsrf
       ENDDO !klon
    ENDIF !(MOD(itap, lmt_pas).EQ.0) THEN

    tslab = tmp_tslab
    seaice = tmp_seaice
  END SUBROUTINE interfoce_slab

end module interfoce_slab_m
