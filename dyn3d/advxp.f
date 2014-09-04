
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advxp.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE advxp(limit, dtx, pbaru, sm, s0, ssx, sy, sz, ssxx, ssxy, ssxz, &
    syy, syz, szz, ntra)
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  IMPLICIT NONE
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! second-order moments (SOM) advection of tracer in X direction  C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! parametres principaux du modele


  INTEGER ntra
  ! PARAMETER (ntra = 1)

  ! definition de la grille du modele

  REAL dtx
  REAL, INTENT (IN) :: pbaru(iip1, jjp1, llm)

  ! moments: SM  total mass in each grid box
  ! S0  mass of tracer in each grid box
  ! Si  1rst order moment in i direction
  ! Sij 2nd  order moment in i and j directions

  REAL sm(iip1, jjp1, llm), s0(iip1, jjp1, llm, ntra)
  REAL ssx(iip1, jjp1, llm, ntra), sy(iip1, jjp1, llm, ntra), &
    sz(iip1, jjp1, llm, ntra)
  REAL ssxx(iip1, jjp1, llm, ntra), ssxy(iip1, jjp1, llm, ntra), &
    ssxz(iip1, jjp1, llm, ntra), syy(iip1, jjp1, llm, ntra), &
    syz(iip1, jjp1, llm, ntra), szz(iip1, jjp1, llm, ntra)

  ! Local :
  ! -------

  ! mass fluxes across the boundaries (UGRI,VGRI,WGRI)
  ! mass fluxes in kg
  ! declaration :

  REAL ugri(iip1, jjp1, llm)

  ! Rem : VGRI et WGRI ne sont pas utilises dans
  ! cette subroutine ( advection en x uniquement )


  ! Tij are the moments for the current latitude and level

  REAL tm(iim)
  REAL t0(iim, ntra), tx(iim, ntra)
  REAL ty(iim, ntra), tz(iim, ntra)
  REAL txx(iim, ntra), txy(iim, ntra)
  REAL txz(iim, ntra), tyy(iim, ntra)
  REAL tyz(iim, ntra), tzz(iim, ntra)

  ! the moments F are similarly defined and used as temporary
  ! storage for portions of the grid boxes in transit

  REAL fm(iim)
  REAL f0(iim, ntra), fx(iim, ntra)
  REAL fy(iim, ntra), fz(iim, ntra)
  REAL fxx(iim, ntra), fxy(iim, ntra)
  REAL fxz(iim, ntra), fyy(iim, ntra)
  REAL fyz(iim, ntra), fzz(iim, ntra)

  ! work arrays

  REAL alf(iim), alf1(iim), alfq(iim), alf1q(iim)
  REAL alf2(iim), alf3(iim), alf4(iim)

  REAL smnew(iim), uext(iim)
  REAL sqi, sqf
  REAL temptm
  REAL slpmax
  REAL s1max, s1new, s2new

  LOGICAL limit
  INTEGER num(jjp1), lonk, numk
  INTEGER lon, lati, latf, niv
  INTEGER i, i2, i3, j, jv, l, k

  lon = iim
  lati = 2
  latf = jjm
  niv = llm

  ! *** Test : diagnostique de la qtite totale de traceur
  ! dans l'atmosphere avant l'advection

  sqi = 0.
  sqf = 0.

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        sqi = sqi + s0(i, j, l, ntra)
      END DO
    END DO
  END DO
  PRINT *, '------ DIAG DANS ADVX2 - ENTREE -----'
  PRINT *, 'sqi=', sqi
  ! test
  ! -------------------------------------
  DO j = 1, jjp1
    num(j) = 1
  END DO
  ! DO l=1,llm
  ! NUM(2,l)=6
  ! NUM(3,l)=6
  ! NUM(jjm-1,l)=6
  ! NUM(jjm,l)=6
  ! ENDDO
  ! DO j=2,6
  ! NUM(j)=12
  ! ENDDO
  ! DO j=jjm-5,jjm-1
  ! NUM(j)=12
  ! ENDDO

  ! Interface : adaptation nouveau modele
  ! -------------------------------------

  ! ---------------------------------------------------------
  ! Conversion des flux de masses en kg/s
  ! pbaru est en N/s d'ou :
  ! ugri est en kg/s

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        ugri(i, j, llm+1-l) = pbaru(i, j, l)
      END DO
    END DO
  END DO

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
            txx(i, jv) = 0.
            txy(i, jv) = 0.
            txz(i, jv) = 0.
            tyy(i, jv) = 0.
            tyz(i, jv) = 0.
            tzz(i, jv) = 0.
          END DO
        END DO

        DO i2 = 1, numk

          DO i = 1, lonk
            i3 = (i-1)*numk + i2
            tm(i) = tm(i) + sm(i3, k, l)
            alf(i) = sm(i3, k, l)/tm(i)
            alf1(i) = 1. - alf(i)
            alfq(i) = alf(i)*alf(i)
            alf1q(i) = alf1(i)*alf1(i)
            alf2(i) = alf1(i) - alf(i)
            alf3(i) = alf(i)*alf1(i)
          END DO

          DO jv = 1, ntra
            DO i = 1, lonk
              i3 = (i-1)*numk + i2
              temptm = -alf(i)*t0(i, jv) + alf1(i)*s0(i3, k, l, jv)
              t0(i, jv) = t0(i, jv) + s0(i3, k, l, jv)
              txx(i, jv) = alfq(i)*ssxx(i3, k, l, jv) + alf1q(i)*txx(i, jv) + &
                5.*(alf3(i)*(ssx(i3,k,l,jv)-tx(i,jv))+alf2(i)*temptm)
              tx(i, jv) = alf(i)*ssx(i3, k, l, jv) + alf1(i)*tx(i, jv) + &
                3.*temptm
              txy(i, jv) = alf(i)*ssxy(i3, k, l, jv) + alf1(i)*txy(i, jv) + &
                3.*(alf1(i)*sy(i3,k,l,jv)-alf(i)*ty(i,jv))
              txz(i, jv) = alf(i)*ssxz(i3, k, l, jv) + alf1(i)*txz(i, jv) + &
                3.*(alf1(i)*sz(i3,k,l,jv)-alf(i)*tz(i,jv))
              ty(i, jv) = ty(i, jv) + sy(i3, k, l, jv)
              tz(i, jv) = tz(i, jv) + sz(i3, k, l, jv)
              tyy(i, jv) = tyy(i, jv) + syy(i3, k, l, jv)
              tyz(i, jv) = tyz(i, jv) + syz(i3, k, l, jv)
              tzz(i, jv) = tzz(i, jv) + szz(i3, k, l, jv)
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
            tx(i, jv) = ssx(i, k, l, jv)
            ty(i, jv) = sy(i, k, l, jv)
            tz(i, jv) = sz(i, k, l, jv)
            txx(i, jv) = ssxx(i, k, l, jv)
            txy(i, jv) = ssxy(i, k, l, jv)
            txz(i, jv) = ssxz(i, k, l, jv)
            tyy(i, jv) = syy(i, k, l, jv)
            tyz(i, jv) = syz(i, k, l, jv)
            tzz(i, jv) = szz(i, k, l, jv)
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
          IF (t0(i,jv)>0.) THEN
            slpmax = t0(i, jv)
            s1max = 1.5*slpmax
            s1new = amin1(s1max, amax1(-s1max,tx(i,jv)))
            s2new = amin1(2.*slpmax-abs(s1new)/3., amax1(abs( &
              s1new)-slpmax,txx(i,jv)))
            tx(i, jv) = s1new
            txx(i, jv) = s2new
            txy(i, jv) = amin1(slpmax, amax1(-slpmax,txy(i,jv)))
            txz(i, jv) = amin1(slpmax, amax1(-slpmax,txz(i,jv)))
          ELSE
            tx(i, jv) = 0.
            txx(i, jv) = 0.
            txy(i, jv) = 0.
            txz(i, jv) = 0.
          END IF
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
        alf2(i) = alf1(i) - alf(i)
        alf3(i) = alf(i)*alfq(i)
        alf4(i) = alf1(i)*alf1q(i)
      END DO

      DO jv = 1, ntra
        DO i = 1, lonk - 1

          IF (uext(i)<0.) THEN

            f0(i, jv) = alf(i)*(t0(i+1,jv)-alf1(i)*(tx(i+1, &
              jv)-alf2(i)*txx(i+1,jv)))
            fx(i, jv) = alfq(i)*(tx(i+1,jv)-3.*alf1(i)*txx(i+1,jv))
            fxx(i, jv) = alf3(i)*txx(i+1, jv)
            fy(i, jv) = alf(i)*(ty(i+1,jv)-alf1(i)*txy(i+1,jv))
            fz(i, jv) = alf(i)*(tz(i+1,jv)-alf1(i)*txz(i+1,jv))
            fxy(i, jv) = alfq(i)*txy(i+1, jv)
            fxz(i, jv) = alfq(i)*txz(i+1, jv)
            fyy(i, jv) = alf(i)*tyy(i+1, jv)
            fyz(i, jv) = alf(i)*tyz(i+1, jv)
            fzz(i, jv) = alf(i)*tzz(i+1, jv)

            t0(i+1, jv) = t0(i+1, jv) - f0(i, jv)
            tx(i+1, jv) = alf1q(i)*(tx(i+1,jv)+3.*alf(i)*txx(i+1,jv))
            txx(i+1, jv) = alf4(i)*txx(i+1, jv)
            ty(i+1, jv) = ty(i+1, jv) - fy(i, jv)
            tz(i+1, jv) = tz(i+1, jv) - fz(i, jv)
            tyy(i+1, jv) = tyy(i+1, jv) - fyy(i, jv)
            tyz(i+1, jv) = tyz(i+1, jv) - fyz(i, jv)
            tzz(i+1, jv) = tzz(i+1, jv) - fzz(i, jv)
            txy(i+1, jv) = alf1q(i)*txy(i+1, jv)
            txz(i+1, jv) = alf1q(i)*txz(i+1, jv)

          END IF

        END DO
      END DO

      i = lonk
      IF (uext(i)<0.) THEN

        DO jv = 1, ntra

          f0(i, jv) = alf(i)*(t0(1,jv)-alf1(i)*(tx(1,jv)-alf2(i)*txx(1,jv)))
          fx(i, jv) = alfq(i)*(tx(1,jv)-3.*alf1(i)*txx(1,jv))
          fxx(i, jv) = alf3(i)*txx(1, jv)
          fy(i, jv) = alf(i)*(ty(1,jv)-alf1(i)*txy(1,jv))
          fz(i, jv) = alf(i)*(tz(1,jv)-alf1(i)*txz(1,jv))
          fxy(i, jv) = alfq(i)*txy(1, jv)
          fxz(i, jv) = alfq(i)*txz(1, jv)
          fyy(i, jv) = alf(i)*tyy(1, jv)
          fyz(i, jv) = alf(i)*tyz(1, jv)
          fzz(i, jv) = alf(i)*tzz(1, jv)

          t0(1, jv) = t0(1, jv) - f0(i, jv)
          tx(1, jv) = alf1q(i)*(tx(1,jv)+3.*alf(i)*txx(1,jv))
          txx(1, jv) = alf4(i)*txx(1, jv)
          ty(1, jv) = ty(1, jv) - fy(i, jv)
          tz(1, jv) = tz(1, jv) - fz(i, jv)
          tyy(1, jv) = tyy(1, jv) - fyy(i, jv)
          tyz(1, jv) = tyz(1, jv) - fyz(i, jv)
          tzz(1, jv) = tzz(1, jv) - fzz(i, jv)
          txy(1, jv) = alf1q(i)*txy(1, jv)
          txz(1, jv) = alf1q(i)*txz(1, jv)

        END DO

      END IF

      DO jv = 1, ntra
        DO i = 1, lonk

          IF (uext(i)>=0.) THEN

            f0(i, jv) = alf(i)*(t0(i,jv)+alf1(i)*(tx(i,jv)+alf2(i)*txx(i, &
              jv)))
            fx(i, jv) = alfq(i)*(tx(i,jv)+3.*alf1(i)*txx(i,jv))
            fxx(i, jv) = alf3(i)*txx(i, jv)
            fy(i, jv) = alf(i)*(ty(i,jv)+alf1(i)*txy(i,jv))
            fz(i, jv) = alf(i)*(tz(i,jv)+alf1(i)*txz(i,jv))
            fxy(i, jv) = alfq(i)*txy(i, jv)
            fxz(i, jv) = alfq(i)*txz(i, jv)
            fyy(i, jv) = alf(i)*tyy(i, jv)
            fyz(i, jv) = alf(i)*tyz(i, jv)
            fzz(i, jv) = alf(i)*tzz(i, jv)

            t0(i, jv) = t0(i, jv) - f0(i, jv)
            tx(i, jv) = alf1q(i)*(tx(i,jv)-3.*alf(i)*txx(i,jv))
            txx(i, jv) = alf4(i)*txx(i, jv)
            ty(i, jv) = ty(i, jv) - fy(i, jv)
            tz(i, jv) = tz(i, jv) - fz(i, jv)
            tyy(i, jv) = tyy(i, jv) - fyy(i, jv)
            tyz(i, jv) = tyz(i, jv) - fyz(i, jv)
            tzz(i, jv) = tzz(i, jv) - fzz(i, jv)
            txy(i, jv) = alf1q(i)*txy(i, jv)
            txz(i, jv) = alf1q(i)*txz(i, jv)

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
        alfq(i) = alf(i)*alf(i)
        alf1q(i) = alf1(i)*alf1(i)
        alf2(i) = alf1(i) - alf(i)
        alf3(i) = alf(i)*alf1(i)
      END DO

      DO jv = 1, ntra
        DO i = 1, lonk

          IF (uext(i)<0.) THEN

            temptm = -alf(i)*t0(i, jv) + alf1(i)*f0(i, jv)
            t0(i, jv) = t0(i, jv) + f0(i, jv)
            txx(i, jv) = alfq(i)*fxx(i, jv) + alf1q(i)*txx(i, jv) + &
              5.*(alf3(i)*(fx(i,jv)-tx(i,jv))+alf2(i)*temptm)
            tx(i, jv) = alf(i)*fx(i, jv) + alf1(i)*tx(i, jv) + 3.*temptm
            txy(i, jv) = alf(i)*fxy(i, jv) + alf1(i)*txy(i, jv) + &
              3.*(alf1(i)*fy(i,jv)-alf(i)*ty(i,jv))
            txz(i, jv) = alf(i)*fxz(i, jv) + alf1(i)*txz(i, jv) + &
              3.*(alf1(i)*fz(i,jv)-alf(i)*tz(i,jv))
            ty(i, jv) = ty(i, jv) + fy(i, jv)
            tz(i, jv) = tz(i, jv) + fz(i, jv)
            tyy(i, jv) = tyy(i, jv) + fyy(i, jv)
            tyz(i, jv) = tyz(i, jv) + fyz(i, jv)
            tzz(i, jv) = tzz(i, jv) + fzz(i, jv)

          END IF

        END DO
      END DO

      DO jv = 1, ntra
        DO i = 1, lonk - 1

          IF (uext(i)>=0.) THEN

            temptm = alf(i)*t0(i+1, jv) - alf1(i)*f0(i, jv)
            t0(i+1, jv) = t0(i+1, jv) + f0(i, jv)
            txx(i+1, jv) = alfq(i)*fxx(i, jv) + alf1q(i)*txx(i+1, jv) + &
              5.*(alf3(i)*(tx(i+1,jv)-fx(i,jv))-alf2(i)*temptm)
            tx(i+1, jv) = alf(i)*fx(i, jv) + alf1(i)*tx(i+1, jv) + 3.*temptm
            txy(i+1, jv) = alf(i)*fxy(i, jv) + alf1(i)*txy(i+1, jv) + &
              3.*(alf(i)*ty(i+1,jv)-alf1(i)*fy(i,jv))
            txz(i+1, jv) = alf(i)*fxz(i, jv) + alf1(i)*txz(i+1, jv) + &
              3.*(alf(i)*tz(i+1,jv)-alf1(i)*fz(i,jv))
            ty(i+1, jv) = ty(i+1, jv) + fy(i, jv)
            tz(i+1, jv) = tz(i+1, jv) + fz(i, jv)
            tyy(i+1, jv) = tyy(i+1, jv) + fyy(i, jv)
            tyz(i+1, jv) = tyz(i+1, jv) + fyz(i, jv)
            tzz(i+1, jv) = tzz(i+1, jv) + fzz(i, jv)

          END IF

        END DO
      END DO

      i = lonk
      IF (uext(i)>=0.) THEN
        DO jv = 1, ntra
          temptm = alf(i)*t0(1, jv) - alf1(i)*f0(i, jv)
          t0(1, jv) = t0(1, jv) + f0(i, jv)
          txx(1, jv) = alfq(i)*fxx(i, jv) + alf1q(i)*txx(1, jv) + &
            5.*(alf3(i)*(tx(1,jv)-fx(i,jv))-alf2(i)*temptm)
          tx(1, jv) = alf(i)*fx(i, jv) + alf1(i)*tx(1, jv) + 3.*temptm
          txy(1, jv) = alf(i)*fxy(i, jv) + alf1(i)*txy(1, jv) + &
            3.*(alf(i)*ty(1,jv)-alf1(i)*fy(i,jv))
          txz(1, jv) = alf(i)*fxz(i, jv) + alf1(i)*txz(1, jv) + &
            3.*(alf(i)*tz(1,jv)-alf1(i)*fz(i,jv))
          ty(1, jv) = ty(1, jv) + fy(i, jv)
          tz(1, jv) = tz(1, jv) + fz(i, jv)
          tyy(1, jv) = tyy(1, jv) + fyy(i, jv)
          tyz(1, jv) = tyz(1, jv) + fyz(i, jv)
          tzz(1, jv) = tzz(1, jv) + fzz(i, jv)
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
            alf2(i) = alf1(i) - alf(i)
            alf3(i) = alf(i)*alfq(i)
            alf4(i) = alf1(i)*alf1q(i)

          END DO

          DO jv = 1, ntra
            DO i = 1, lonk

              i3 = i2 + (i-1)*numk
              s0(i3, k, l, jv) = alf(i)*(t0(i,jv)-alf1(i)*(tx(i, &
                jv)-alf2(i)*txx(i,jv)))
              ssx(i3, k, l, jv) = alfq(i)*(tx(i,jv)-3.*alf1(i)*txx(i,jv))
              ssxx(i3, k, l, jv) = alf3(i)*txx(i, jv)
              sy(i3, k, l, jv) = alf(i)*(ty(i,jv)-alf1(i)*txy(i,jv))
              sz(i3, k, l, jv) = alf(i)*(tz(i,jv)-alf1(i)*txz(i,jv))
              ssxy(i3, k, l, jv) = alfq(i)*txy(i, jv)
              ssxz(i3, k, l, jv) = alfq(i)*txz(i, jv)
              syy(i3, k, l, jv) = alf(i)*tyy(i, jv)
              syz(i3, k, l, jv) = alf(i)*tyz(i, jv)
              szz(i3, k, l, jv) = alf(i)*tzz(i, jv)

              ! reajusts moments remaining in the box

              t0(i, jv) = t0(i, jv) - s0(i3, k, l, jv)
              tx(i, jv) = alf1q(i)*(tx(i,jv)+3.*alf(i)*txx(i,jv))
              txx(i, jv) = alf4(i)*txx(i, jv)
              ty(i, jv) = ty(i, jv) - sy(i3, k, l, jv)
              tz(i, jv) = tz(i, jv) - sz(i3, k, l, jv)
              tyy(i, jv) = tyy(i, jv) - syy(i3, k, l, jv)
              tyz(i, jv) = tyz(i, jv) - syz(i3, k, l, jv)
              tzz(i, jv) = tzz(i, jv) - szz(i3, k, l, jv)
              txy(i, jv) = alf1q(i)*txy(i, jv)
              txz(i, jv) = alf1q(i)*txz(i, jv)

            END DO
          END DO

        END DO

      ELSE

        DO i = 1, lon
          sm(i, k, l) = tm(i)
        END DO
        DO jv = 1, ntra
          DO i = 1, lon
            s0(i, k, l, jv) = t0(i, jv)
            ssx(i, k, l, jv) = tx(i, jv)
            sy(i, k, l, jv) = ty(i, jv)
            sz(i, k, l, jv) = tz(i, jv)
            ssxx(i, k, l, jv) = txx(i, jv)
            ssxy(i, k, l, jv) = txy(i, jv)
            ssxz(i, k, l, jv) = txz(i, jv)
            syy(i, k, l, jv) = tyy(i, jv)
            syz(i, k, l, jv) = tyz(i, jv)
            szz(i, k, l, jv) = tzz(i, jv)
          END DO
        END DO

      END IF

    END DO
  END DO

  ! ---------- bouclage cyclique

  DO l = 1, llm
    DO j = 1, jjp1
      sm(iip1, j, l) = sm(1, j, l)
      s0(iip1, j, l, ntra) = s0(1, j, l, ntra)
      ssx(iip1, j, l, ntra) = ssx(1, j, l, ntra)
      sy(iip1, j, l, ntra) = sy(1, j, l, ntra)
      sz(iip1, j, l, ntra) = sz(1, j, l, ntra)
    END DO
  END DO

  ! ----------- qqtite totale de traceur dans tte l'atmosphere
  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        sqf = sqf + s0(i, j, l, ntra)
      END DO
    END DO
  END DO

  PRINT *, '------ DIAG DANS ADVX2 - SORTIE -----'
  PRINT *, 'sqf=', sqf
  ! -------------------------------------------------------------
  RETURN
END SUBROUTINE advxp
