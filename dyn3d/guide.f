MODULE guide_m

  ! From dyn3d/guide.F, version 1.3 2005/05/25 13:10:09
  ! and dyn3d/guide.h, version 1.1.1.1 2004/05/19 12:53:06

  IMPLICIT NONE

  REAL aire_min, aire_max

CONTAINS

  SUBROUTINE guide(itau, ucov, vcov, teta, q, ps)

    ! Author: F.Hourdin

    USE comconst, ONLY: cpp, daysec, dtvr, kappa
    USE comgeom, ONLY: aire, rlatu, rlonv
    USE conf_gcm_m, ONLY: day_step, iperiod
    use conf_guide_m, only: conf_guide, guide_u, guide_v, guide_t, guide_q, &
         ncep, ini_anal, tau_min_u, tau_max_u, tau_min_v, tau_max_v, &
         tau_min_t, tau_max_t, tau_min_q, tau_max_q, tau_min_p, tau_max_p, &
         online
    USE dimens_m, ONLY: iim, jjm, llm
    USE disvert_m, ONLY: ap, bp, preff, presnivs
    use dump2d_m, only: dump2d
    USE exner_hyb_m, ONLY: exner_hyb
    USE inigrads_m, ONLY: inigrads
    use massdair_m, only: massdair
    use netcdf, only: nf90_nowrite, nf90_close, nf90_inq_dimid
    use netcdf95, only: nf95_inquire_dimension, nf95_open
    use nr_util, only: pi
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1, jjp1, llmp1
    USE q_sat_m, ONLY: q_sat
    use read_reanalyse_m, only: read_reanalyse
    USE serre, ONLY: clat, clon
    use tau2alpha_m, only: tau2alpha, dxdys

    INTEGER, INTENT(IN):: itau

    ! variables dynamiques

    REAL, intent(inout):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) vent covariant
    REAL, intent(inout):: vcov(:, :, :) ! (iim + 1, jjm, llm) ! vent covariant

    REAL, intent(inout):: teta(iim + 1, jjm + 1, llm) ! température potentielle 
    REAL, intent(inout):: q(iim + 1, jjm + 1, llm)
    REAL, intent(in):: ps(:, :) ! (iim + 1, jjm + 1) pression au sol

    ! Local:

    ! variables dynamiques pour les reanalyses.

    REAL, save:: ucovrea1(iim + 1, jjm + 1, llm), vcovrea1(iim + 1, jjm, llm)
    ! vents covariants reanalyses

    REAL, save:: tetarea1(iim + 1, jjm + 1, llm) ! temp pot reales
    REAL, save:: qrea1(iim + 1, jjm + 1, llm) ! temp pot reales

    REAL, save:: ucovrea2(iim + 1, jjm + 1, llm), vcovrea2(iim + 1, jjm, llm)
    ! vents covariants reanalyses

    REAL, save:: tetarea2(iim + 1, jjm + 1, llm) ! temp pot reales
    REAL, save:: qrea2(iim + 1, jjm + 1, llm) ! temp pot reales
    REAL, save:: masserea2(ip1jmp1, llm) ! masse

    ! alpha determine la part des injections de donnees a chaque etape
    ! alpha=1 signifie pas d'injection
    ! alpha=0 signifie injection totale
    REAL, save:: alpha_q(iim + 1, jjm + 1)
    REAL, save:: alpha_t(iim + 1, jjm + 1), alpha_p(ip1jmp1)
    REAL, save:: alpha_u(iim + 1, jjm + 1), alpha_v(iim + 1, jjm)

    INTEGER, save:: step_rea, count_no_rea

    INTEGER ilon, ilat
    REAL factt ! pas de temps entre deux appels au guidage, en fraction de jour
    real ztau(iim + 1, jjm + 1)

    INTEGER ij, l
    INTEGER ncidpl, status
    INTEGER rcod, rid
    REAL tau
    INTEGER, SAVE:: nlev

    ! TEST SUR QSAT
    REAL p(iim + 1, jjm + 1, llmp1)
    real pk(iim + 1, jjm + 1, llm), pks(iim + 1, jjm + 1)

    REAL qsat(iim + 1, jjm + 1, llm)

    INTEGER, parameter:: igrads = 2
    REAL:: dtgrads = 100.

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: guide'

    first_call: IF (itau == 0) THEN
       CALL conf_guide
       CALL inigrads(igrads, rlonv, 180. / pi, -180., 180., rlatu, -90., &
            90., 180. / pi, presnivs, 1., dtgrads, 'guide', 'dyn_zon ')

       IF (online) THEN
          ! Constantes de temps de rappel en jour

          ! coordonnees du centre du zoom
          CALL coordij(clon, clat, ilon, ilat)
          ! aire de la maille au centre du zoom
          aire_min = aire(ilon+(ilat - 1) * iip1)
          ! aire maximale de la maille
          aire_max = 0.
          DO ij = 1, ip1jmp1
             aire_max = max(aire_max, aire(ij))
          END DO

          factt = dtvr * iperiod / daysec

          CALL tau2alpha(3, iip1, jjm, factt, tau_min_v, tau_max_v, alpha_v)
          CALL tau2alpha(2, iip1, jjp1, factt, tau_min_u, tau_max_u, alpha_u)
          CALL tau2alpha(1, iip1, jjp1, factt, tau_min_t, tau_max_t, alpha_t)
          CALL tau2alpha(1, iip1, jjp1, factt, tau_min_p, tau_max_p, alpha_p)
          CALL tau2alpha(1, iip1, jjp1, factt, tau_min_q, tau_max_q, alpha_q)

          CALL dump2d(iip1, jjp1, aire, 'AIRE MAILLe ')
          CALL dump2d(iip1, jjp1, alpha_u, 'COEFF U ')
          CALL dump2d(iip1, jjp1, alpha_t, 'COEFF T ')
       ELSE
          ! Cas ou on force exactement par les variables analysees
          alpha_t = 0.
          alpha_u = 0.
          alpha_v = 0.
          alpha_p = 0.
       END IF

       step_rea = 1
       count_no_rea = 0
       ncidpl = -99

       ! lecture d'un fichier netcdf pour determiner le nombre de niveaux
       if (guide_u) call nf95_open('u.nc',Nf90_NOWRITe,ncidpl)
       if (guide_v) call nf95_open('v.nc',nf90_nowrite,ncidpl)
       if (guide_T) call nf95_open('T.nc',nf90_nowrite,ncidpl)
       if (guide_Q) call nf95_open('hur.nc',nf90_nowrite, ncidpl)

       IF (ncep) THEN
          status = nf90_inq_dimid(ncidpl, 'LEVEL', rid)
       ELSE
          status = nf90_inq_dimid(ncidpl, 'PRESSURE', rid)
       END IF
       call nf95_inquire_dimension(ncidpl, rid, nclen=nlev)
       PRINT *, 'nlev', nlev
       rcod = nf90_close(ncidpl)
       ! Lecture du premier etat des reanalyses.
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

       PRINT *, 'LECTURE REANALYSES, pas ', step_rea, 'apres ', &
            count_no_rea, ' non lectures'
       step_rea = step_rea + 1
       CALL read_reanalyse(step_rea, ps, ucovrea2, vcovrea2, tetarea2, qrea2, &
            masserea2, nlev)
       qrea2 = max(qrea2, 0.1)
       factt = dtvr * iperiod / daysec
       ztau = factt / max(alpha_t, 1E-10)
       CALL wrgrads(igrads, 1, aire, 'aire ', 'aire ')
       CALL wrgrads(igrads, 1, dxdys, 'dxdy ', 'dxdy ')
       CALL wrgrads(igrads, 1, alpha_u, 'au ', 'au ')
       CALL wrgrads(igrads, 1, alpha_t, 'at ', 'at ')
       CALL wrgrads(igrads, 1, ztau, 'taut ', 'taut ')
       CALL wrgrads(igrads, llm, ucov, 'u ', 'u ')
       CALL wrgrads(igrads, llm, ucovrea2, 'ua ', 'ua ')
       CALL wrgrads(igrads, llm, teta, 'T ', 'T ')
       CALL wrgrads(igrads, llm, tetarea2, 'Ta ', 'Ta ')
       CALL wrgrads(igrads, llm, qrea2, 'Qa ', 'Qa ')
       CALL wrgrads(igrads, llm, q, 'Q ', 'Q ')
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
