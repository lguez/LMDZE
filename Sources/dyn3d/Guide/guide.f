MODULE guide_m

  ! From dyn3d/guide.F, version 1.3 2005/05/25 13:10:09
  ! and dyn3d/guide.h, version 1.1.1.1 2004/05/19 12:53:06

  IMPLICIT NONE

CONTAINS

  SUBROUTINE guide(itau, ucov, vcov, teta, q, ps)

    ! Author: F.Hourdin

    USE comconst, ONLY: cpp, kappa
    USE conf_gcm_m, ONLY: day_step
    use conf_guide_m, only: guide_u, guide_v, guide_t, guide_q, ncep, &
         ini_anal, tau_min_u, tau_max_u, tau_min_v, tau_max_v, tau_min_t, &
         tau_max_t, tau_min_q, tau_max_q, online, factt
    USE dimens_m, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: ap, bp, preff
    use dynetat0_m, only: grossismx, grossismy, rlatu, rlatv
    USE exner_hyb_m, ONLY: exner_hyb
    use init_tau2alpha_m, only: init_tau2alpha
    use netcdf, only: nf90_nowrite
    use netcdf95, only: nf95_close, nf95_inq_dimid, nf95_inquire_dimension, &
         nf95_open
    use nr_util, only: pi
    USE paramet_m, ONLY: iip1, ip1jmp1, jjp1, llmp1
    USE q_sat_m, ONLY: q_sat
    use read_reanalyse_m, only: read_reanalyse
    use tau2alpha_m, only: tau2alpha
    use writefield_m, only: writefield

    INTEGER, INTENT(IN):: itau
    REAL, intent(inout):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, intent(inout):: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant

    REAL, intent(inout):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! température potentielle 

    REAL, intent(inout):: q(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol

    ! Local:

    ! variables dynamiques pour les réanalyses

    REAL, save:: ucovrea1(iim + 1, jjm + 1, llm), vcovrea1(iim + 1, jjm, llm)
    ! vents covariants reanalyses

    REAL, save:: tetarea1(iim + 1, jjm + 1, llm) ! temp pot reales
    REAL, save:: qrea1(iim + 1, jjm + 1, llm) ! temp pot reales

    REAL, save:: ucovrea2(iim + 1, jjm + 1, llm), vcovrea2(iim + 1, jjm, llm)
    ! vents covariants reanalyses

    REAL, save:: tetarea2(iim + 1, jjm + 1, llm) ! temp pot reales
    REAL, save:: qrea2(iim + 1, jjm + 1, llm) ! temp pot reales
    REAL, save:: masserea2(ip1jmp1, llm) ! masse

    ! alpha détermine la part des injections de données à chaque étape
    ! alpha=0 signifie pas d'injection
    ! alpha=1 signifie injection totale
    REAL, save:: alpha_q(iim + 1, jjm + 1)
    REAL, save:: alpha_t(iim + 1, jjm + 1)
    REAL, save:: alpha_u(iim + 1, jjm + 1), alpha_v(iim + 1, jjm)

    INTEGER, save:: step_rea, count_no_rea

    INTEGER l
    INTEGER ncid, dimid
    REAL tau
    INTEGER, SAVE:: nlev

    ! TEST SUR QSAT
    REAL p(iim + 1, jjm + 1, llmp1)
    real pk(iim + 1, jjm + 1, llm), pks(iim + 1, jjm + 1)
    REAL qsat(iim + 1, jjm + 1, llm)

    REAL dxdys(iip1, jjp1), dxdyu(iip1, jjp1), dxdyv(iip1, jjm)

    !-----------------------------------------------------------------------

    !!PRINT *, 'Call sequence information: guide'

    first_call: IF (itau == 0) THEN
       IF (online) THEN
          IF (abs(grossismx - 1.) < 0.1 .OR. abs(grossismy - 1.) < 0.1) THEN
             ! grille regulière
             if (guide_u) alpha_u = factt / tau_max_u
             if (guide_v) alpha_v = factt / tau_max_v
             if (guide_t) alpha_t = factt / tau_max_t
             if (guide_q) alpha_q = factt / tau_max_q
          else
             call init_tau2alpha(dxdys, dxdyu, dxdyv)

             if (guide_u) then
                CALL tau2alpha(dxdyu, rlatu, tau_min_u, tau_max_u, alpha_u)
                CALL writefield("alpha_u", alpha_u)
             end if

             if (guide_v) then
                CALL tau2alpha(dxdyv, rlatv, tau_min_v, tau_max_v, alpha_v)
                CALL writefield("alpha_v", alpha_v)
             end if

             if (guide_t) then
                CALL tau2alpha(dxdys, rlatu, tau_min_t, tau_max_t, alpha_t)
                CALL writefield("alpha_t", alpha_t)
             end if

             if (guide_q)  then
                CALL tau2alpha(dxdys, rlatu, tau_min_q, tau_max_q, alpha_q)
                CALL writefield("alpha_q", alpha_q)
             end if
          end IF
       ELSE
          ! Cas où on force exactement par les variables analysées
          if (guide_u) alpha_u = 1.
          if (guide_v) alpha_v = 1.
          if (guide_t) alpha_t = 1.
          if (guide_q) alpha_q = 1.
       END IF

       step_rea = 1
       count_no_rea = 0

       ! lecture d'un fichier netcdf pour determiner le nombre de niveaux :

       if (guide_u) then
          call nf95_open('u.nc',Nf90_NOWRITe,ncid)
       else if (guide_v) then
          call nf95_open('v.nc',nf90_nowrite,ncid)
       else if (guide_T) then
          call nf95_open('T.nc',nf90_nowrite,ncid)
       else
          call nf95_open('hur.nc',nf90_nowrite, ncid)
       end if

       IF (ncep) THEN
          call nf95_inq_dimid(ncid, 'LEVEL', dimid)
       ELSE
          call nf95_inq_dimid(ncid, 'PRESSURE', dimid)
       END IF
       call nf95_inquire_dimension(ncid, dimid, nclen=nlev)
       PRINT *, 'nlev = ', nlev
       call nf95_close(ncid)

       ! Lecture du premier état des réanalyses :
       CALL read_reanalyse(1, ps, ucovrea2, vcovrea2, tetarea2, qrea2, &
            masserea2, nlev)
       qrea2 = max(qrea2, 0.1)
    END IF first_call

    ! IMPORTATION DES VENTS, PRESSION ET TEMPERATURE REELS:

    ! Nudging fields are given 4 times per day:
    IF (mod(itau, day_step / 4) == 0) THEN
       vcovrea1 = vcovrea2
       ucovrea1 = ucovrea2
       tetarea1 = tetarea2
       qrea1 = qrea2

       PRINT *, 'Lecture fichiers guidage, pas ', step_rea, 'apres ', &
            count_no_rea, ' non lectures'
       step_rea = step_rea + 1
       CALL read_reanalyse(step_rea, ps, ucovrea2, vcovrea2, tetarea2, qrea2, &
            masserea2, nlev)
       qrea2 = max(qrea2, 0.1)

       if (guide_u) then
          CALL writefield("ucov", ucov)
          CALL writefield("ucovrea2", ucovrea2)
       end if

       if (guide_t) then
          CALL writefield("teta", teta)
          CALL writefield("tetarea2", tetarea2)
       end if

       if (guide_q) then
          CALL writefield("qrea2", qrea2)
          CALL writefield("q", q)
       end if
    ELSE
       count_no_rea = count_no_rea + 1
    END IF

    ! Guidage

    tau = mod(real(itau) / real(day_step / 4), 1.)

    ! x_gcm = a * x_gcm + (1 - a) * x_reanalyses

    IF (guide_u) THEN
       IF (itau == 0 .AND. ini_anal) then
          ucov = ucovrea1
       else
          forall (l = 1: llm) ucov(:, :, l) = (1. - alpha_u) * ucov(:, :, l) &
               + alpha_u * ((1. - tau) * ucovrea1(:, :, l) &
               + tau * ucovrea2(:, :, l))
       end IF
    END IF

    IF (guide_t) THEN
       IF (itau == 0 .AND. ini_anal) then
          teta = tetarea1
       else
          forall (l = 1: llm) teta(:, :, l) = (1. - alpha_t) * teta(:, :, l) &
               + alpha_t * ((1. - tau) * tetarea1(:, :, l) &
               + tau * tetarea2(:, :, l))
       end IF
    END IF

    IF (guide_q) THEN
       ! Calcul de l'humidité saturante :
       forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * ps
       CALL exner_hyb(ps, p, pks, pk)
       qsat = q_sat(pk * teta / cpp, preff * (pk / cpp)**(1. / kappa))

       ! humidité relative en % -> humidité spécifique
       IF (itau == 0 .AND. ini_anal) then
          q = qsat * qrea1 * 0.01
       else
          forall (l = 1: llm) q(:, :, l) = (1. - alpha_q) * q(:, :, l) &
               + alpha_q * (qsat(:, :, l) * ((1. - tau) * qrea1(:, :, l) &
               + tau * qrea2(:, :, l)) * 0.01)
       end IF
    END IF

    IF (guide_v) THEN
       IF (itau == 0 .AND. ini_anal) then
          vcov = vcovrea1
       else
          forall (l = 1: llm) vcov(:, :, l) = (1. - alpha_v) * vcov(:, :, l) &
               + alpha_v * ((1. - tau) * vcovrea1(:, :, l) &
               + tau * vcovrea2(:, :, l))
       end IF
    END IF

  END SUBROUTINE guide

END MODULE guide_m