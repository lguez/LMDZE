
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advx.F,v 1.2 2005/05/25 13:10:09
! fairhead Exp $

SUBROUTINE advx(limit, dtx, pbaru, sm, s0, sx, sy, sz, lati, latf)
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  IMPLICIT NONE

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! first-order moments (FOM) advection of tracer in X direction  C
  ! C
  ! Source : Pascal Simon (Meteo,CNRM)                            C
  ! Adaptation : A.Armengaud (LGGE) juin 94                       C
  ! C
  ! limit,dtx,pbaru,pbarv,sm,s0,sx,sy,sz                       C
  ! sont des arguments d'entree pour le s-pg...                   C
  ! C
  ! sm,s0,sx,sy,sz                                                C
  ! sont les arguments de sortie pour le s-pg                     C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! parametres principaux du modele


  ! Arguments :
  ! -----------
  ! dtx : frequence fictive d'appel du transport
  ! pbaru, pbarv : flux de masse en x et y en Pa.m2.s-1

  INTEGER ntra
  PARAMETER (ntra=1)

  ! ATTENTION partout ou on trouve ntra, insertion de boucle
  ! possible dans l'avenir.

  REAL dtx
  REAL, INTENT (IN) :: pbaru(iip1, jjp1, llm)

  ! moments: SM  total mass in each grid box
  ! S0  mass of tracer in each grid box
  ! Si  1rst order moment in i direction

  REAL sm(iip1, jjp1, llm), s0(iip1, jjp1, llm, ntra)
  REAL sx(iip1, jjp1, llm, ntra), sy(iip1, jjp1, llm, ntra)
  REAL sz(iip1, jjp1, llm, ntra)

  ! Local :
  ! -------

  ! mass fluxes across the boundaries (UGRI,VGRI,WGRI)
  ! mass fluxes in kg
  ! declaration :

  REAL ugri(iip1, jjp1, llm)

  ! Rem : VGRI et WGRI ne sont pas utilises dans
  ! cette subroutine ( advection en x uniquement )

  ! Ti are the moments for the current latitude and level

  REAL tm(iim)
  REAL t0(iim, ntra), tx(iim, ntra)
  REAL ty(iim, ntra), tz(iim, ntra)
  REAL temptm ! just a temporary variable

  ! the moments F are similarly defined and used as temporary
  ! storage for portions of the grid boxes in transit

  REAL fm(iim)
  REAL f0(iim, ntra), fx(iim, ntra)
  REAL fy(iim, ntra), fz(iim, ntra)

  ! work arrays

  REAL alf(iim), alf1(iim), alfq(iim), alf1q(iim)

  REAL smnew(iim), uext(iim)

  REAL sqi, sqf

  LOGICAL limit
  INTEGER num(jjp1), lonk, numk
  INTEGER lon, lati, latf, niv
  INTEGER i, i2, i3, j, jv, l, k, itrac

  lon = iim
  niv = llm

  ! *** Test de passage d'arguments ******


  ! -------------------------------------
  DO j = 1, jjp1
    num(j) = 1
  END DO
  sqi = 0.
  sqf = 0.

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        ! IM 240305            sqi = sqi + S0(i,j,l,9)
        sqi = sqi + s0(i, j, l, ntra)
      END DO
    END DO
  END DO
  PRINT *, '-------- DIAG DANS ADVX - ENTREE ---------'
  PRINT *, 'sqi=', sqi


  ! Interface : adaptation nouveau modele
  ! -------------------------------------

  ! ---------------------------------------------------------
  ! Conversion des flux de masses en kg/s
  ! pbaru est en N/s d'ou :
  ! ugri est en kg/s

  DO l = 1, llm
    DO j = 1, jjm + 1
      DO i = 1, iip1
        ! ugri (i,j,llm+1-l) = pbaru (i,j,l) * ( dsig(l) / g )
        ugri(i, j, llm+1-l) = pbaru(i, j, l)
      END DO
    END DO
  END DO


  ! ---------------------------------------------------------
  ! ---------------------------------------------------------
  ! ---------------------------------------------------------

  ! start here

  ! boucle principale sur les niveaux et les latitudes

  DO l = 1, niv
    DO k = lati, latf

      ! initialisation

      ! program assumes periodic boundaries in X

      DO i = 2, lon
        smnew(i) = sm(i, k, l) + (ugri(i-1,k,l)-ugri(i,k,l))*dtx
      END DO
      smnew(1) = sm(1, k, l) + (ugri(lon,k,l)-ugri(1,k,l))*dtx

      ! modifications for extended polar zones

      numk = num(k)
      lonk = lon/numk

      IF (numk>1) THEN

        DO i = 1, lon
          tm(i) = 0.
        END DO
        DO jv = 1, ntra
          DO i = 1, lon
            t0(i, jv) = 0.
            tx(i, jv) = 0.
            ty(i, jv) = 0.
            tz(i, jv) = 0.
          END DO
        END DO

        DO i2 = 1, numk

          DO i = 1, lonk
            i3 = (i-1)*numk + i2
            tm(i) = tm(i) + sm(i3, k, l)
            alf(i) = sm(i3, k, l)/tm(i)
            alf1(i) = 1. - alf(i)
          END DO

          DO jv = 1, ntra
            DO i = 1, lonk
              i3 = (i-1)*numk + i2
              temptm = -alf(i)*t0(i, jv) + alf1(i)*s0(i3, k, l, jv)
              t0(i, jv) = t0(i, jv) + s0(i3, k, l, jv)
              tx(i, jv) = alf(i)*sx(i3, k, l, jv) + alf1(i)*tx(i, jv) + &
                3.*temptm
              ty(i, jv) = ty(i, jv) + sy(i3, k, l, jv)
              tz(i, jv) = tz(i, jv) + sz(i3, k, l, jv)
            END DO
          END DO

        END DO

      ELSE

        DO i = 1, lon
          tm(i) = sm(i, k, l)
        END DO
        DO jv = 1, ntra
          DO i = 1, lon
            t0(i, jv) = s0(i, k, l, jv)
            tx(i, jv) = sx(i, k, l, jv)
            ty(i, jv) = sy(i, k, l, jv)
            tz(i, jv) = sz(i, k, l, jv)
          END DO
        END DO

      END IF

      DO i = 1, lonk
        uext(i) = ugri(i*numk, k, l)
      END DO

      ! place limits on appropriate moments before transport
      ! (if flux-limiting is to be applied)

      IF (.NOT. limit) GO TO 13

      DO jv = 1, ntra
        DO i = 1, lonk
          tx(i, jv) = sign(amin1(amax1(t0(i,jv),0.),abs(tx(i,jv))), tx(i,jv))
        END DO
      END DO

13    CONTINUE

      ! calculate flux and moments between adjacent boxes
      ! 1- create temporary moments/masses for partial boxes in transit
      ! 2- reajusts moments remaining in the box

      ! flux from IP to I if U(I).lt.0

      DO i = 1, lonk - 1
        IF (uext(i)<0.) THEN
          fm(i) = -uext(i)*dtx
          alf(i) = fm(i)/tm(i+1)
          tm(i+1) = tm(i+1) - fm(i)
        END IF
      END DO

      i = lonk
      IF (uext(i)<0.) THEN
        fm(i) = -uext(i)*dtx
        alf(i) = fm(i)/tm(1)
        tm(1) = tm(1) - fm(i)
      END IF

      ! flux from I to IP if U(I).gt.0

      DO i = 1, lonk
        IF (uext(i)>=0.) THEN
          fm(i) = uext(i)*dtx
          alf(i) = fm(i)/tm(i)
          tm(i) = tm(i) - fm(i)
        END IF
      END DO

      DO i = 1, lonk
        alfq(i) = alf(i)*alf(i)
        alf1(i) = 1. - alf(i)
        alf1q(i) = alf1(i)*alf1(i)
      END DO

      DO jv = 1, ntra
        DO i = 1, lonk - 1

          IF (uext(i)<0.) THEN

            f0(i, jv) = alf(i)*(t0(i+1,jv)-alf1(i)*tx(i+1,jv))
            fx(i, jv) = alfq(i)*tx(i+1, jv)
            fy(i, jv) = alf(i)*ty(i+1, jv)
            fz(i, jv) = alf(i)*tz(i+1, jv)

            t0(i+1, jv) = t0(i+1, jv) - f0(i, jv)
            tx(i+1, jv) = alf1q(i)*tx(i+1, jv)
            ty(i+1, jv) = ty(i+1, jv) - fy(i, jv)
            tz(i+1, jv) = tz(i+1, jv) - fz(i, jv)

          END IF

        END DO
      END DO

      i = lonk
      IF (uext(i)<0.) THEN

        DO jv = 1, ntra

          f0(i, jv) = alf(i)*(t0(1,jv)-alf1(i)*tx(1,jv))
          fx(i, jv) = alfq(i)*tx(1, jv)
          fy(i, jv) = alf(i)*ty(1, jv)
          fz(i, jv) = alf(i)*tz(1, jv)

          t0(1, jv) = t0(1, jv) - f0(i, jv)
          tx(1, jv) = alf1q(i)*tx(1, jv)
          ty(1, jv) = ty(1, jv) - fy(i, jv)
          tz(1, jv) = tz(1, jv) - fz(i, jv)

        END DO

      END IF

      DO jv = 1, ntra
        DO i = 1, lonk

          IF (uext(i)>=0.) THEN

            f0(i, jv) = alf(i)*(t0(i,jv)+alf1(i)*tx(i,jv))
            fx(i, jv) = alfq(i)*tx(i, jv)
            fy(i, jv) = alf(i)*ty(i, jv)
            fz(i, jv) = alf(i)*tz(i, jv)

            t0(i, jv) = t0(i, jv) - f0(i, jv)
            tx(i, jv) = alf1q(i)*tx(i, jv)
            ty(i, jv) = ty(i, jv) - fy(i, jv)
            tz(i, jv) = tz(i, jv) - fz(i, jv)

          END IF

        END DO
      END DO

      ! puts the temporary moments Fi into appropriate neighboring boxes

      DO i = 1, lonk
        IF (uext(i)<0.) THEN
          tm(i) = tm(i) + fm(i)
          alf(i) = fm(i)/tm(i)
        END IF
      END DO

      DO i = 1, lonk - 1
        IF (uext(i)>=0.) THEN
          tm(i+1) = tm(i+1) + fm(i)
          alf(i) = fm(i)/tm(i+1)
        END IF
      END DO

      i = lonk
      IF (uext(i)>=0.) THEN
        tm(1) = tm(1) + fm(i)
        alf(i) = fm(i)/tm(1)
      END IF

      DO i = 1, lonk
        alf1(i) = 1. - alf(i)
      END DO

      DO jv = 1, ntra
        DO i = 1, lonk

          IF (uext(i)<0.) THEN

            temptm = -alf(i)*t0(i, jv) + alf1(i)*f0(i, jv)
            t0(i, jv) = t0(i, jv) + f0(i, jv)
            tx(i, jv) = alf(i)*fx(i, jv) + alf1(i)*tx(i, jv) + 3.*temptm
            ty(i, jv) = ty(i, jv) + fy(i, jv)
            tz(i, jv) = tz(i, jv) + fz(i, jv)

          END IF

        END DO
      END DO

      DO jv = 1, ntra
        DO i = 1, lonk - 1

          IF (uext(i)>=0.) THEN

            temptm = alf(i)*t0(i+1, jv) - alf1(i)*f0(i, jv)
            t0(i+1, jv) = t0(i+1, jv) + f0(i, jv)
            tx(i+1, jv) = alf(i)*fx(i, jv) + alf1(i)*tx(i+1, jv) + 3.*temptm
            ty(i+1, jv) = ty(i+1, jv) + fy(i, jv)
            tz(i+1, jv) = tz(i+1, jv) + fz(i, jv)

          END IF

        END DO
      END DO

      i = lonk
      IF (uext(i)>=0.) THEN
        DO jv = 1, ntra
          temptm = alf(i)*t0(1, jv) - alf1(i)*f0(i, jv)
          t0(1, jv) = t0(1, jv) + f0(i, jv)
          tx(1, jv) = alf(i)*fx(i, jv) + alf1(i)*tx(1, jv) + 3.*temptm
          ty(1, jv) = ty(1, jv) + fy(i, jv)
          tz(1, jv) = tz(1, jv) + fz(i, jv)
        END DO
      END IF

      ! retour aux mailles d'origine (passage des Tij aux Sij)

      IF (numk>1) THEN

        DO i2 = 1, numk

          DO i = 1, lonk

            i3 = i2 + (i-1)*numk
            sm(i3, k, l) = smnew(i3)
            alf(i) = smnew(i3)/tm(i)
            tm(i) = tm(i) - smnew(i3)

            alfq(i) = alf(i)*alf(i)
            alf1(i) = 1. - alf(i)
            alf1q(i) = alf1(i)*alf1(i)

          END DO
        END DO

        DO jv = 1, ntra
          DO i = 1, lonk

            i3 = i2 + (i-1)*numk
            s0(i3, k, l, jv) = alf(i)*(t0(i,jv)-alf1(i)*tx(i,jv))
            sx(i3, k, l, jv) = alfq(i)*tx(i, jv)
            sy(i3, k, l, jv) = alf(i)*ty(i, jv)
            sz(i3, k, l, jv) = alf(i)*tz(i, jv)

            ! reajusts moments remaining in the box

            t0(i, jv) = t0(i, jv) - s0(i3, k, l, jv)
            tx(i, jv) = alf1q(i)*tx(i, jv)
            ty(i, jv) = ty(i, jv) - sy(i3, k, l, jv)
            tz(i, jv) = tz(i, jv) - sz(i3, k, l, jv)
          END DO
        END DO


      ELSE

        DO i = 1, lon
          sm(i, k, l) = tm(i)
        END DO
        DO jv = 1, ntra
          DO i = 1, lon
            s0(i, k, l, jv) = t0(i, jv)
            sx(i, k, l, jv) = tx(i, jv)
            sy(i, k, l, jv) = ty(i, jv)
            sz(i, k, l, jv) = tz(i, jv)
          END DO
        END DO

      END IF

    END DO
  END DO

  ! ---------- bouclage cyclique
  DO itrac = 1, ntra
    DO l = 1, llm
      DO j = lati, latf
        sm(iip1, j, l) = sm(1, j, l)
        s0(iip1, j, l, itrac) = s0(1, j, l, itrac)
        sx(iip1, j, l, itrac) = sx(1, j, l, itrac)
        sy(iip1, j, l, itrac) = sy(1, j, l, itrac)
        sz(iip1, j, l, itrac) = sz(1, j, l, itrac)
      END DO
    END DO
  END DO

  ! ----------- qqtite totale de traceur dans tte l'atmosphere
  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        ! IM 240405          sqf = sqf + S0(i,j,l,9)
        sqf = sqf + s0(i, j, l, ntra)
      END DO
    END DO
  END DO

  PRINT *, '------ DIAG DANS ADVX - SORTIE -----'
  PRINT *, 'sqf=', sqf
  ! -------------

  RETURN
END SUBROUTINE advx
! _________________________________________________________________
! _________________________________________________________________
