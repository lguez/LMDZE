MODULE guide_m

  ! From dyn3d/guide.F, version 1.3, 2005/05/25 13:10:09
  ! and dyn3d/guide.h, version 1.1.1.1, 2004/05/19 12:53:06

  IMPLICIT NONE

CONTAINS

  SUBROUTINE guide(itau, ucov, vcov, teta, q, ps)

    ! Author: F. Hourdin

    USE comconst, ONLY: cpp, kappa
    USE conf_gcm_m, ONLY: day_step
    use conf_guide_m, only: guide_u, guide_v, guide_t, guide_q, ini_anal, &
         alpha_u, alpha_v, alpha_t, alpha_q
    USE dimens_m, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: ap, bp, preff
    USE exner_hyb_m, ONLY: exner_hyb
    USE q_sat_m, ONLY: q_sat
    use read_reanalyse_m, only: read_reanalyse
    use writefield_m, only: writefield

    INTEGER, INTENT(IN):: itau
    REAL, intent(inout):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, intent(inout):: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant

    REAL, intent(inout):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! température potentielle 

    REAL, intent(inout):: q(:, :, :) ! (iim + 1, jjm + 1, llm)
    REAL, intent(in):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol

    ! Local:

    ! Variables dynamiques pour les réanalyses

    REAL, save:: ucovrea1(iim + 1, jjm + 1, llm), vcovrea1(iim + 1, jjm, llm)
    ! vents covariants r\'eanalyses

    REAL, save:: tetarea1(iim + 1, jjm + 1, llm)
    ! potential temperture from reanalysis

    REAL, save:: qrea1(iim + 1, jjm + 1, llm)

    REAL, save:: ucovrea2(iim + 1, jjm + 1, llm), vcovrea2(iim + 1, jjm, llm)
    ! vents covariants reanalyses

    REAL, save:: tetarea2(iim + 1, jjm + 1, llm)
    ! potential temperture from reanalysis

    REAL, save:: qrea2(iim + 1, jjm + 1, llm)

    INTEGER l
    REAL tau

    ! TEST SUR QSAT
    REAL p(iim + 1, jjm + 1, llm + 1)
    real pk(iim + 1, jjm + 1, llm), pks(iim + 1, jjm + 1)
    REAL qsat(iim + 1, jjm + 1, llm)

    !-----------------------------------------------------------------------

    IF (itau == 0) THEN
       ! Lecture du premier état des réanalyses :
       CALL read_reanalyse(ps, ucovrea2, vcovrea2, tetarea2, qrea2)
       qrea2 = max(qrea2, 0.1)

       if (ini_anal) then
          IF (guide_u) ucov = ucovrea2
          IF (guide_v) vcov = vcovrea2
          IF (guide_t) teta = tetarea2

          IF (guide_q) then
             ! Calcul de l'humidité saturante :
             forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * ps
             CALL exner_hyb(ps, p, pks, pk)
             q = q_sat(pk * teta / cpp, preff * (pk / cpp)**(1. / kappa)) &
                  * qrea2 * 0.01
          end IF
       end if
    END IF

    ! Importation des vents, pression et temp\'erature r\'eels :

    ! Nudging fields are given 4 times per day:
    IF (mod(itau, day_step / 4) == 0) THEN
       vcovrea1 = vcovrea2
       ucovrea1 = ucovrea2
       tetarea1 = tetarea2
       qrea1 = qrea2

       CALL read_reanalyse(ps, ucovrea2, vcovrea2, tetarea2, qrea2)
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
    END IF

    ! Guidage

    tau = mod(real(itau) / real(day_step / 4), 1.)

    ! x_gcm = a * x_gcm + (1 - a) * x_reanalyses

    IF (guide_u) forall (l = 1: llm) ucov(:, :, l) = (1. - alpha_u) &
         * ucov(:, :, l) + alpha_u * ((1. - tau) * ucovrea1(:, :, l) + tau &
         * ucovrea2(:, :, l))

    IF (guide_v) forall (l = 1: llm) vcov(:, :, l) = (1. - alpha_v) &
         * vcov(:, :, l) + alpha_v * ((1. - tau) * vcovrea1(:, :, l) + tau &
         * vcovrea2(:, :, l))

    IF (guide_t) forall (l = 1: llm) teta(:, :, l) = (1. - alpha_t) &
         * teta(:, :, l) + alpha_t * ((1. - tau) * tetarea1(:, :, l) + tau &
         * tetarea2(:, :, l))

    IF (guide_q) THEN
       ! Calcul de l'humidité saturante :
       forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * ps
       CALL exner_hyb(ps, p, pks, pk)
       qsat = q_sat(pk * teta / cpp, preff * (pk / cpp)**(1. / kappa))

       ! humidité relative en % -> humidité spécifique
       forall (l = 1: llm) q(:, :, l) = (1. - alpha_q) * q(:, :, l) &
            + alpha_q * (qsat(:, :, l) * ((1. - tau) * qrea1(:, :, l) &
            + tau * qrea2(:, :, l)) * 0.01)
    END IF

  END SUBROUTINE guide

END MODULE guide_m
