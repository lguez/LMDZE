MODULE guide_m

  ! From dyn3d/guide.F, v 1.3 2005/05/25 13:10:09
  ! and dyn3d/guide.h, v 1.1.1.1 2004/05/19 12:53:06

  REAL :: tau_min_u, tau_max_u
  REAL :: tau_min_v, tau_max_v
  REAL :: tau_min_t, tau_max_t
  REAL :: tau_min_q, tau_max_q
  REAL :: tau_min_p, tau_max_p
  REAL :: aire_min, aire_max


  LOGICAL :: guide_u, guide_v, guide_t, guide_q, guide_p
  REAL :: lat_min_guide, lat_max_guide

  LOGICAL :: ncep, ini_anal
  INTEGER :: online

CONTAINS

SUBROUTINE guide(itau, ucov, vcov, teta, q, masse, ps)

  USE dimens_m, ONLY : jjm, llm
  USE paramet_m, ONLY : iip1, ip1jm, ip1jmp1, jjp1, llmp1
  USE comconst, ONLY : cpp, daysec, dtvr, kappa, pi
  USE comvert, ONLY : ap, bp, preff, presnivs
  USE conf_gcm_m, ONLY : day_step, iperiod
  USE comgeom, ONLY : aire, rlatu, rlonv
  USE serre, ONLY : clat, clon
  USE q_sat_m, ONLY : q_sat
  USE exner_hyb_m, ONLY : exner_hyb
  USE pression_m, ONLY : pression
  USE inigrads_m, ONLY : inigrads
  use netcdf, only: nf90_nowrite, nf90_open, nf90_close

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  !      ......   Version  du 10/01/98    ..........

  !             avec  coordonnees  verticales hybrides 
  !   avec nouveaux operat. dissipation * ( gradiv2, divgrad2, nxgraro2 )

  !=======================================================================

  !   Auteur:  F.Hourdin
  !   -------

  !   Objet:
  !   ------

  !   GCM LMD nouvelle grille

  !=======================================================================

  !   ...  Dans inigeom , nouveaux calculs pour les elongations  cu , cv 
  !        et possibilite d'appeler une fonction f(y)  a derivee tangente 
  !        hyperbolique a la  place de la fonction a derivee sinusoidale.         

  !   ...  Possibilite de choisir le shema de Van-leer pour l'advection de
  !         q  , en faisant iadv = 10  dans   traceur  (29/04/97) .

  !-----------------------------------------------------------------------
  !   Declarations:
  !   -------------


  !   variables dynamiques
  REAL :: vcov(ip1jm, llm), ucov(ip1jmp1, llm) ! vents covariants
  REAL :: teta(ip1jmp1, llm) ! temperature potentielle 
  REAL :: q(ip1jmp1, llm) ! temperature potentielle 
  REAL :: ps(ip1jmp1) ! pression  au sol
  REAL :: masse(ip1jmp1, llm) ! masse d'air

  !   common passe pour des sorties
  REAL :: dxdys(iip1, jjp1), dxdyu(iip1, jjp1), dxdyv(iip1, jjm)
  COMMON /comdxdy/dxdys, dxdyu, dxdyv

  !   variables dynamiques pour les reanalyses.
  REAL :: ucovrea1(ip1jmp1, llm), vcovrea1(ip1jm, llm) !vts cov reas
  REAL :: tetarea1(ip1jmp1, llm) ! temp pot  reales
  REAL :: qrea1(ip1jmp1, llm) ! temp pot  reales
  REAL :: psrea1(ip1jmp1) ! ps
  REAL :: ucovrea2(ip1jmp1, llm), vcovrea2(ip1jm, llm) !vts cov reas
  REAL :: tetarea2(ip1jmp1, llm) ! temp pot  reales
  REAL :: qrea2(ip1jmp1, llm) ! temp pot  reales
  REAL :: masserea2(ip1jmp1, llm) ! masse
  REAL :: psrea2(ip1jmp1) ! ps

  REAL :: alpha_q(ip1jmp1)
  REAL :: alpha_t(ip1jmp1), alpha_p(ip1jmp1)
  REAL :: alpha_u(ip1jmp1), alpha_v(ip1jm)
  REAL :: dday_step, toto, reste, itau_test
  INTEGER :: step_rea, count_no_rea

  !IM 180305   real aire_min, aire_max
  INTEGER :: ilon, ilat
  REAL :: factt, ztau(ip1jmp1)

  INTEGER, INTENT (IN) :: itau
  INTEGER :: ij, l
  INTEGER :: ncidpl, varidpl, nlev, status
  INTEGER :: rcod, rid
  REAL :: ditau, tau, a
  SAVE nlev

  !  TEST SUR QSAT
  REAL :: p(ip1jmp1, llmp1), pk(ip1jmp1, llm), pks(ip1jmp1)
  REAL :: pkf(ip1jmp1, llm)
  REAL :: pres(ip1jmp1, llm)

  REAL :: qsat(ip1jmp1, llm)
  REAL :: unskap
  REAL :: tnat(ip1jmp1, llm)
  !cccccccccccccccc


  LOGICAL :: first
  SAVE first
  DATA first/ .TRUE./

  SAVE ucovrea1, vcovrea1, tetarea1, psrea1, qrea1
  SAVE ucovrea2, vcovrea2, tetarea2, masserea2, psrea2, qrea2

  SAVE alpha_t, alpha_q, alpha_u, alpha_v, alpha_p, itau_test
  SAVE step_rea, count_no_rea

  CHARACTER (10) :: file
  INTEGER :: igrads
  REAL :: dtgrads
  SAVE igrads, dtgrads
  DATA igrads, dtgrads/2, 100./

  PRINT *, 'Call sequence information: guide'

  !-----------------------------------------------------------------------
  ! calcul de l'humidite saturante
  !-----------------------------------------------------------------------
  CALL pression(ip1jmp1, ap, bp, ps, p)
  CALL massdair(p, masse)
  PRINT *, 'OK1'
  CALL exner_hyb(ps, p, pks, pk, pkf)
  PRINT *, 'OK2'
  tnat(:, :) = pk(:, :)*teta(:, :)/cpp
  PRINT *, 'OK3'
  unskap = 1./kappa
  pres(:, :) = preff*(pk(:, :)/cpp)**unskap
  PRINT *, 'OK4'
  qsat = q_sat(tnat, pres)

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !   initialisations pour la lecture des reanalyses.
  !    alpha determine la part des injections de donnees a chaque etape
  !    alpha=1 signifie pas d'injection
  !    alpha=0 signifie injection totale
  !-----------------------------------------------------------------------

  PRINT *, 'ONLINE=', online
  IF (online==-1) THEN
     RETURN
  END IF

  IF (first) THEN

     PRINT *, 'initialisation du guide '
     CALL conf_guide
     PRINT *, 'apres conf_guide'

     file = 'guide'
     CALL inigrads(igrads, rlonv, 180./pi, -180., 180., rlatu, -90., 90., &
          180./pi, presnivs, 1., dtgrads, file, 'dyn_zon ')

     PRINT *, &
          '1: en-ligne, 0: hors-ligne (x=x_rea), -1: climat (x=x_gcm)'

     IF (online==-1) RETURN
     IF (online==1) THEN

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !  Constantes de temps de rappel en jour
        !  0.1 c'est en gros 2h30. 
        !  1e10  est une constante infinie donc en gros pas de guidage
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !   coordonnees du centre du zoom
        CALL coordij(clon, clat, ilon, ilat)
        !   aire de la maille au centre du zoom
        aire_min = aire(ilon+(ilat-1)*iip1)
        !   aire maximale de la maille
        aire_max = 0.
        DO ij = 1, ip1jmp1
           aire_max = max(aire_max, aire(ij))
        END DO
        !  factt = pas de temps en fraction de jour
        factt = dtvr*iperiod/daysec

        !     subroutine tau2alpha(type, im, jm, factt, taumin, taumax, alpha)
        CALL tau2alpha(3, iip1, jjm, factt, tau_min_v, tau_max_v, alpha_v)
        CALL tau2alpha(2, iip1, jjp1, factt, tau_min_u, tau_max_u, alpha_u)
        CALL tau2alpha(1, iip1, jjp1, factt, tau_min_t, tau_max_t, alpha_t)
        CALL tau2alpha(1, iip1, jjp1, factt, tau_min_p, tau_max_p, alpha_p)
        CALL tau2alpha(1, iip1, jjp1, factt, tau_min_q, tau_max_q, alpha_q)

        CALL dump2d(iip1, jjp1, aire, 'AIRE MAILLe ')
        CALL dump2d(iip1, jjp1, alpha_u, 'COEFF U   ')
        CALL dump2d(iip1, jjp1, alpha_t, 'COEFF T   ')

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !   Cas ou on force exactement par les variables analysees
     ELSE
        alpha_t = 0.
        alpha_u = 0.
        alpha_v = 0.
        alpha_p = 0.
        !           physic=.false.
     END IF

     itau_test = 1001
     step_rea = 1
     count_no_rea = 0
     ncidpl = -99

     !    itau_test    montre si l'importation a deja ete faite au rang itau
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
        status = nf_inq_dimid(ncidpl, 'LEVEL', rid)
     ELSE
        status = nf_inq_dimid(ncidpl, 'PRESSURE', rid)
     END IF
     status = nf_inq_dimlen(ncidpl, rid, nlev)
     PRINT *, 'nlev', nlev
     rcod = nf90_close(ncidpl)
     !   Lecture du premier etat des reanalyses.
     CALL read_reanalyse(1, ps, ucovrea2, vcovrea2, tetarea2, qrea2, masserea2, &
          psrea2, 1, nlev)
     qrea2(:, :) = max(qrea2(:, :), 0.1)


     !-----------------------------------------------------------------------
     !   Debut de l'integration temporelle:
     !   ----------------------------------

  END IF ! first

  !-----------------------------------------------------------------------
  !----- IMPORTATION DES VENTS, PRESSION ET TEMPERATURE REELS:
  !-----------------------------------------------------------------------

  ditau = real(itau)
  dday_step = real(day_step)
  WRITE (*, *) 'ditau, dday_step'
  WRITE (*, *) ditau, dday_step
  toto = 4*ditau/dday_step
  reste = toto - aint(toto)
  !     write(*, *)'toto, reste', toto, reste

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
        CALL read_reanalyse(step_rea, ps, ucovrea2, vcovrea2, tetarea2, qrea2, &
             masserea2, psrea2, 1, nlev)
        qrea2(:, :) = max(qrea2(:, :), 0.1)
        factt = dtvr*iperiod/daysec
        ztau(:) = factt/max(alpha_t(:), 1.E-10)
        CALL wrgrads(igrads, 1, aire, 'aire      ', 'aire      ')
        CALL wrgrads(igrads, 1, dxdys, 'dxdy      ', 'dxdy      ')
        CALL wrgrads(igrads, 1, alpha_u, 'au        ', 'au        ')
        CALL wrgrads(igrads, 1, alpha_t, 'at        ', 'at        ')
        CALL wrgrads(igrads, 1, ztau, 'taut      ', 'taut      ')
        CALL wrgrads(igrads, llm, ucov, 'u         ', 'u         ')
        CALL wrgrads(igrads, llm, ucovrea2, 'ua        ', 'ua        ')
        CALL wrgrads(igrads, llm, teta, 'T         ', 'T         ')
        CALL wrgrads(igrads, llm, tetarea2, 'Ta        ', 'Ta        ')
        CALL wrgrads(igrads, llm, qrea2, 'Qa        ', 'Qa        ')
        CALL wrgrads(igrads, llm, q, 'Q         ', 'Q         ')

        CALL wrgrads(igrads, llm, qsat, 'QSAT      ', 'QSAT      ')

     END IF
  ELSE
     count_no_rea = count_no_rea + 1
  END IF

  !-----------------------------------------------------------------------
  !   Guidage
  !    x_gcm = a * x_gcm + (1-a) * x_reanalyses
  !-----------------------------------------------------------------------

  IF (ini_anal) PRINT *, 'ATTENTION !!! ON PART DU GUIDAGE'

  ditau = real(itau)
  dday_step = real(day_step)


  tau = 4*ditau/dday_step
  tau = tau - aint(tau)

  !  ucov
  IF (guide_u) THEN
     DO l = 1, llm
        DO ij = 1, ip1jmp1
           a = (1.-tau)*ucovrea1(ij, l) + tau*ucovrea2(ij, l)
           ucov(ij, l) = (1.-alpha_u(ij))*ucov(ij, l) + alpha_u(ij)*a
           IF (first .AND. ini_anal) ucov(ij, l) = a
        END DO
     END DO
  END IF

  !  teta
  IF (guide_t) THEN
     DO l = 1, llm
        DO ij = 1, ip1jmp1
           a = (1.-tau)*tetarea1(ij, l) + tau*tetarea2(ij, l)
           teta(ij, l) = (1.-alpha_t(ij))*teta(ij, l) + alpha_t(ij)*a
           IF (first .AND. ini_anal) teta(ij, l) = a
        END DO
     END DO
  END IF

  !  P
  IF (guide_p) THEN
     DO ij = 1, ip1jmp1
        a = (1.-tau)*psrea1(ij) + tau*psrea2(ij)
        ps(ij) = (1.-alpha_p(ij))*ps(ij) + alpha_p(ij)*a
        IF (first .AND. ini_anal) ps(ij) = a
     END DO
     CALL pression(ip1jmp1, ap, bp, ps, p)
     CALL massdair(p, masse)
  END IF


  !  q
  IF (guide_q) THEN
     DO l = 1, llm
        DO ij = 1, ip1jmp1
           a = (1.-tau)*qrea1(ij, l) + tau*qrea2(ij, l)
           !   hum relative en % -> hum specif
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

  !     call dump2d(iip1, jjp1, tetarea1, 'TETA REA 1     ')
  !     call dump2d(iip1, jjp1, tetarea2, 'TETA REA 2     ')
  !     call dump2d(iip1, jjp1, teta, 'TETA           ')

  first = .FALSE.

  RETURN
END SUBROUTINE guide

  !=======================================================================
  SUBROUTINE tau2alpha(type, pim, pjm, factt, taumin, taumax, alpha)
    !=======================================================================

    USE dimens_m, ONLY : iim, jjm
    USE paramet_m, ONLY : iip1, jjp1
    USE comconst, ONLY : pi
    USE comgeom, ONLY : cu_2d, cv_2d, rlatu, rlatv
    USE serre, ONLY : clat, clon, grossismx, grossismy
    IMPLICIT NONE

    !   arguments :
    INTEGER :: type
    INTEGER :: pim, pjm
    REAL :: factt, taumin, taumax
    REAL :: dxdy_, alpha(pim, pjm)
    REAL :: dxdy_min, dxdy_max

    !  local :
    REAL :: alphamin, alphamax, gamma, xi
    SAVE gamma
    INTEGER :: i, j, ilon, ilat

    LOGICAL :: first
    SAVE first
    DATA first/ .TRUE./

    REAL :: zdx(iip1, jjp1), zdy(iip1, jjp1)

    REAL :: zlat
    REAL :: dxdys(iip1, jjp1), dxdyu(iip1, jjp1), dxdyv(iip1, jjm)
    COMMON /comdxdy/dxdys, dxdyu, dxdyv

    IF (first) THEN
       DO j = 2, jjm
          DO i = 2, iip1
             zdx(i, j) = 0.5*(cu_2d(i-1, j)+cu_2d(i, j))/cos(rlatu(j))
          END DO
          zdx(1, j) = zdx(iip1, j)
       END DO
       DO j = 2, jjm
          DO i = 1, iip1
             zdy(i, j) = 0.5*(cv_2d(i, j-1)+cv_2d(i, j))
          END DO
       END DO
       DO i = 1, iip1
          zdx(i, 1) = zdx(i, 2)
          zdx(i, jjp1) = zdx(i, jjm)
          zdy(i, 1) = zdy(i, 2)
          zdy(i, jjp1) = zdy(i, jjm)
       END DO
       DO j = 1, jjp1
          DO i = 1, iip1
             dxdys(i, j) = sqrt(zdx(i, j)*zdx(i, j)+zdy(i, j)*zdy(i, j))
          END DO
       END DO
       DO j = 1, jjp1
          DO i = 1, iim
             dxdyu(i, j) = 0.5*(dxdys(i, j)+dxdys(i+1, j))
          END DO
          dxdyu(iip1, j) = dxdyu(1, j)
       END DO
       DO j = 1, jjm
          DO i = 1, iip1
             dxdyv(i, j) = 0.5*(dxdys(i, j)+dxdys(i+1, j))
          END DO
       END DO

       CALL dump2d(iip1, jjp1, dxdys, 'DX2DY2 SCAL  ')
       CALL dump2d(iip1, jjp1, dxdyu, 'DX2DY2 U     ')
       CALL dump2d(iip1, jjp1, dxdyv, 'DX2DY2 v     ')

       !   coordonnees du centre du zoom
       CALL coordij(clon, clat, ilon, ilat)
       !   aire de la maille au centre du zoom
       dxdy_min = dxdys(ilon, ilat)
       !   dxdy maximale de la maille
       dxdy_max = 0.
       DO j = 1, jjp1
          DO i = 1, iip1
             dxdy_max = max(dxdy_max, dxdys(i, j))
          END DO
       END DO

       IF (abs(grossismx-1.)<0.1 .OR. abs(grossismy-1.)<0.1) THEN
          PRINT *, 'ATTENTION modele peu zoome'
          PRINT *, 'ATTENTION on prend une constante de guidage cste'
          gamma = 0.
       ELSE
          gamma = (dxdy_max-2.*dxdy_min)/(dxdy_max-dxdy_min)
          PRINT *, 'gamma=', gamma
          IF (gamma<1.E-5) THEN
             PRINT *, 'gamma =', gamma, '<1e-5'
             STOP
          END IF
          PRINT *, 'gamma=', gamma
          gamma = log(0.5)/log(gamma)
       END IF
    END IF

    alphamin = factt/taumax
    alphamax = factt/taumin

    DO j = 1, pjm
       DO i = 1, pim
          IF (type==1) THEN
             dxdy_ = dxdys(i, j)
             zlat = rlatu(j)*180./pi
          ELSE IF (type==2) THEN
             dxdy_ = dxdyu(i, j)
             zlat = rlatu(j)*180./pi
          ELSE IF (type==3) THEN
             dxdy_ = dxdyv(i, j)
             zlat = rlatv(j)*180./pi
          END IF
          IF (abs(grossismx-1.)<0.1 .OR. abs(grossismy-1.)<0.1) THEN
             !  pour une grille reguliere, xi=xxx**0=1 -> alpha=alphamin
             alpha(i, j) = alphamin
          ELSE
             xi = ((dxdy_max-dxdy_)/(dxdy_max-dxdy_min))**gamma
             xi = min(xi, 1.)
             IF (lat_min_guide<=zlat .AND. zlat<=lat_max_guide) THEN
                alpha(i, j) = xi*alphamin + (1.-xi)*alphamax
             ELSE
                alpha(i, j) = 0.
             END IF
          END IF
       END DO
    END DO


    RETURN
  END SUBROUTINE tau2alpha

END MODULE guide_m
