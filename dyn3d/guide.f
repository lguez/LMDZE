MODULE guide_m

  ! From dyn3d/guide.F, version 1.3 2005/05/25 13:10:09
  ! and dyn3d/guide.h, version 1.1.1.1 2004/05/19 12:53:06

  IMPLICIT NONE

  REAL aire_min, aire_max

CONTAINS

  SUBROUTINE guide(itau, ucov, vcov, teta, q, masse, ps)

    ! Author: F.Hourdin

    USE comconst, ONLY: cpp, daysec, dtvr, kappa
    USE comgeom, ONLY: aire, rlatu, rlonv
    USE conf_gcm_m, ONLY: day_step, iperiod
    use conf_guide_m, only: conf_guide, guide_u, guide_v, guide_t, guide_q, &
         ncep, ini_anal, tau_min_u, tau_max_u, tau_min_v, tau_max_v, &
         tau_min_t, tau_max_t, tau_min_q, tau_max_q, tau_min_p, tau_max_p, &
         online
    USE dimens_m, ONLY: jjm, llm
    USE disvert_m, ONLY: ap, bp, preff, presnivs
    USE exner_hyb_m, ONLY: exner_hyb
    USE inigrads_m, ONLY: inigrads
    use massdair_m, only: massdair
    use netcdf, only: nf90_nowrite, nf90_open, nf90_close, nf90_inq_dimid, &
         nf90_inquire_dimension
    use nr_util, only: pi
    USE paramet_m, ONLY: iip1, ip1jm, ip1jmp1, jjp1, llmp1
    USE q_sat_m, ONLY: q_sat
    USE serre, ONLY: clat, clon
    use tau2alpha_m, only: tau2alpha, dxdys

    INTEGER, INTENT(IN):: itau

    ! variables dynamiques
    REAL ucov(ip1jmp1, llm), vcov(ip1jm, llm) ! vents covariants
    REAL, intent(inout):: teta(ip1jmp1, llm) ! temperature potentielle 
    REAL q(ip1jmp1, llm) ! temperature potentielle 
    REAL, intent(out):: masse(ip1jmp1, llm) ! masse d'air
    REAL, intent(in):: ps(ip1jmp1) ! pression au sol

    ! Local:

    ! variables dynamiques pour les reanalyses.
    REAL, save:: ucovrea1(ip1jmp1, llm), vcovrea1(ip1jm, llm) !vts cov reas
    REAL, save:: tetarea1(ip1jmp1, llm) ! temp pot reales
    REAL, save:: qrea1(ip1jmp1, llm) ! temp pot reales
    REAL, save:: ucovrea2(ip1jmp1, llm), vcovrea2(ip1jm, llm) !vts cov reas
    REAL, save:: tetarea2(ip1jmp1, llm) ! temp pot reales
    REAL, save:: qrea2(ip1jmp1, llm) ! temp pot reales
    REAL, save:: masserea2(ip1jmp1, llm) ! masse

    REAL, save:: alpha_q(ip1jmp1)
    REAL, save:: alpha_t(ip1jmp1), alpha_p(ip1jmp1)
    REAL, save:: alpha_u(ip1jmp1), alpha_v(ip1jm)
    REAL dday_step, toto, reste
    real, save:: itau_test
    INTEGER, save:: step_rea, count_no_rea

    INTEGER ilon, ilat
    REAL factt, ztau(ip1jmp1)

    INTEGER ij, l
    INTEGER ncidpl, varidpl, status
    INTEGER rcod, rid
    REAL ditau, tau, a
    INTEGER, SAVE:: nlev

    ! TEST SUR QSAT
    REAL p(ip1jmp1, llmp1), pk(ip1jmp1, llm), pks(ip1jmp1)
    REAL pkf(ip1jmp1, llm)
    REAL pres(ip1jmp1, llm)

    REAL qsat(ip1jmp1, llm)
    REAL unskap
    REAL tnat(ip1jmp1, llm)

    LOGICAL:: first = .TRUE.
    CHARACTER(len=10) file
    INTEGER:: igrads = 2
    REAL:: dtgrads = 100.

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: guide'

    ! calcul de l'humidite saturante

    forall (l = 1: llm + 1) p(:, l) = ap(l) + bp(l) * ps
    CALL massdair(p, masse)
    CALL exner_hyb(ps, p, pks, pk, pkf)
    tnat(:, :) = pk(:, :)*teta(:, :)/cpp
    unskap = 1./kappa
    pres(:, :) = preff*(pk(:, :)/cpp)**unskap
    qsat = q_sat(tnat, pres)

    ! initialisations pour la lecture des reanalyses.
    ! alpha determine la part des injections de donnees a chaque etape
    ! alpha=1 signifie pas d'injection
    ! alpha=0 signifie injection totale

    IF (online==-1) THEN
       RETURN
    END IF

    IF (first) THEN
       CALL conf_guide
       file = 'guide'
       CALL inigrads(igrads, rlonv, 180./pi, -180., 180., rlatu, -90., 90., &
            180./pi, presnivs, 1., dtgrads, file, 'dyn_zon ')
       PRINT *, '1: en-ligne, 0: hors-ligne (x=x_rea), -1: climat (x=x_gcm)'
       IF (online==-1) RETURN

       IF (online==1) THEN
          ! Constantes de temps de rappel en jour
          ! 0.1 c'est en gros 2h30. 
          ! 1e10 est une constante infinie donc en gros pas de guidage

          ! coordonnees du centre du zoom
          CALL coordij(clon, clat, ilon, ilat)
          ! aire de la maille au centre du zoom
          aire_min = aire(ilon+(ilat-1)*iip1)
          ! aire maximale de la maille
          aire_max = 0.
          DO ij = 1, ip1jmp1
             aire_max = max(aire_max, aire(ij))
          END DO
          ! factt = pas de temps en fraction de jour
          factt = dtvr*iperiod/daysec

          CALL tau2alpha(3, iip1, jjm, factt, tau_min_v, tau_max_v, alpha_v)
          CALL tau2alpha(2, iip1, jjp1, factt, tau_min_u, tau_max_u, alpha_u)
          CALL tau2alpha(1, iip1, jjp1, factt, tau_min_t, tau_max_t, alpha_t)
          CALL tau2alpha(1, iip1, jjp1, factt, tau_min_p, tau_max_p, alpha_p)
          CALL tau2alpha(1, iip1, jjp1, factt, tau_min_q, tau_max_q, alpha_q)

          CALL dump2d(iip1, jjp1, aire, 'AIRE MAILLe ')
          CALL dump2d(iip1, jjp1, alpha_u, 'COEFF U ')
          CALL dump2d(iip1, jjp1, alpha_t, 'COEFF T ')

          ! Cas ou on force exactement par les variables analysees
       ELSE
          alpha_t = 0.
          alpha_u = 0.
          alpha_v = 0.
          alpha_p = 0.
          ! physic=.false.
       END IF

       itau_test = 1001
       step_rea = 1
       count_no_rea = 0
       ncidpl = -99

       ! itau_test montre si l'importation a deja ete faite au rang itau
       ! lecture d'un fichier netcdf pour determiner le nombre de niveaux
       if (guide_u) then
          if (ncidpl.eq.-99) rcod=nf90_open('u.nc',Nf90_NOWRITe,ncidpl)
       endif

       if (guide_v) then
          if (ncidpl.eq.-99) rcod=nf90_open('v.nc',nf90_nowrite,ncidpl)
       endif

       if (guide_T) then
          if (ncidpl.eq.-99) rcod=nf90_open('T.nc',nf90_nowrite,ncidpl)
       endif

       if (guide_Q) then
          if (ncidpl.eq.-99) rcod=nf90_open('hur.nc',nf90_nowrite, ncidpl)
       endif

       IF (ncep) THEN
          status = nf90_inq_dimid(ncidpl, 'LEVEL', rid)
       ELSE
          status = nf90_inq_dimid(ncidpl, 'PRESSURE', rid)
       END IF
       status = nf90_inquire_dimension(ncidpl, rid, len=nlev)
       PRINT *, 'nlev', nlev
       rcod = nf90_close(ncidpl)
       ! Lecture du premier etat des reanalyses.
       CALL read_reanalyse(1, ps, ucovrea2, vcovrea2, tetarea2, qrea2, &
            masserea2, 1, nlev)
       qrea2(:, :) = max(qrea2(:, :), 0.1)

       ! Debut de l'integration temporelle:
    END IF ! first

    ! IMPORTATION DES VENTS, PRESSION ET TEMPERATURE REELS:

    ditau = real(itau)
    dday_step = real(day_step)
    WRITE (*, *) 'ditau, dday_step'
    WRITE (*, *) ditau, dday_step
    toto = 4*ditau/dday_step
    reste = toto - aint(toto)

    IF (reste==0.) THEN
       IF (itau_test==itau) THEN
          WRITE (*, *) 'deuxieme passage de advreel a itau=', itau
          STOP
       ELSE
          vcovrea1(:, :) = vcovrea2(:, :)
          ucovrea1(:, :) = ucovrea2(:, :)
          tetarea1(:, :) = tetarea2(:, :)
          qrea1(:, :) = qrea2(:, :)

          PRINT *, 'LECTURE REANALYSES, pas ', step_rea, 'apres ', &
               count_no_rea, ' non lectures'
          step_rea = step_rea + 1
          itau_test = itau
          CALL read_reanalyse(step_rea, ps, ucovrea2, vcovrea2, tetarea2, &
               qrea2, masserea2, 1, nlev)
          qrea2(:, :) = max(qrea2(:, :), 0.1)
          factt = dtvr*iperiod/daysec
          ztau(:) = factt/max(alpha_t(:), 1.E-10)
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

          CALL wrgrads(igrads, llm, qsat, 'QSAT ', 'QSAT ')

       END IF
    ELSE
       count_no_rea = count_no_rea + 1
    END IF

    ! Guidage
    ! x_gcm = a * x_gcm + (1-a) * x_reanalyses

    IF (ini_anal) PRINT *, 'ATTENTION !!! ON PART DU GUIDAGE'

    ditau = real(itau)
    dday_step = real(day_step)

    tau = 4*ditau/dday_step
    tau = tau - aint(tau)

    ! ucov
    IF (guide_u) THEN
       DO l = 1, llm
          DO ij = 1, ip1jmp1
             a = (1.-tau)*ucovrea1(ij, l) + tau*ucovrea2(ij, l)
             ucov(ij, l) = (1.-alpha_u(ij))*ucov(ij, l) + alpha_u(ij)*a
             IF (first .AND. ini_anal) ucov(ij, l) = a
          END DO
       END DO
    END IF

    IF (guide_t) THEN
       DO l = 1, llm
          DO ij = 1, ip1jmp1
             a = (1.-tau)*tetarea1(ij, l) + tau*tetarea2(ij, l)
             teta(ij, l) = (1.-alpha_t(ij))*teta(ij, l) + alpha_t(ij)*a
             IF (first .AND. ini_anal) teta(ij, l) = a
          END DO
       END DO
    END IF

    IF (guide_q) THEN
       DO l = 1, llm
          DO ij = 1, ip1jmp1
             a = (1.-tau)*qrea1(ij, l) + tau*qrea2(ij, l)
             ! hum relative en % -> hum specif
             a = qsat(ij, l)*a*0.01
             q(ij, l) = (1.-alpha_q(ij))*q(ij, l) + alpha_q(ij)*a
             IF (first .AND. ini_anal) q(ij, l) = a
          END DO
       END DO
    END IF

    ! vcov
    IF (guide_v) THEN
       DO l = 1, llm
          DO ij = 1, ip1jm
             a = (1.-tau)*vcovrea1(ij, l) + tau*vcovrea2(ij, l)
             vcov(ij, l) = (1.-alpha_v(ij))*vcov(ij, l) + alpha_v(ij)*a
             IF (first .AND. ini_anal) vcov(ij, l) = a
          END DO
          IF (first .AND. ini_anal) vcov(ij, l) = a
       END DO
    END IF

    first = .FALSE.

  END SUBROUTINE guide

END MODULE guide_m
