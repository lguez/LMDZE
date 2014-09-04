
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/pentes_ini.F,v 1.1.1.1 2004/05/19
! 12:53:07 lmdzadmin Exp $

SUBROUTINE pentes_ini(q, w, masse, pbaru, pbarv, mode)
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE comgeom
  USE nr_util, ONLY: pi
  IMPLICIT NONE

  ! =======================================================================
  ! Adaptation LMDZ:  A.Armengaud (LGGE)
  ! ----------------

  ! ********************************************************************
  ! Transport des traceurs par la methode des pentes
  ! ********************************************************************
  ! Reference possible : Russel. G.L., Lerner J.A.:
  ! A new Finite-Differencing Scheme for Traceur Transport
  ! Equation , Journal of Applied Meteorology, pp 1483-1498,dec. 81
  ! ********************************************************************
  ! q,w,masse,pbaru et pbarv
  ! sont des arguments d'entree  pour le s-pg ....

  ! =======================================================================



  ! Arguments:
  ! ----------
  INTEGER mode
  REAL, INTENT (IN) :: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
  REAL q(iip1, jjp1, llm, 0:3)
  REAL w(ip1jmp1, llm)
  REAL masse(iip1, jjp1, llm)
  ! Local:
  ! ------
  LOGICAL limit
  REAL sm(iip1, jjp1, llm)
  REAL s0(iip1, jjp1, llm), sx(iip1, jjp1, llm)
  REAL sy(iip1, jjp1, llm), sz(iip1, jjp1, llm)
  REAL masn, mass, zz
  INTEGER i, j, l, iq

  ! modif Fred 24 03 96

  REAL sinlon(iip1), sinlondlon(iip1)
  REAL coslon(iip1), coslondlon(iip1)
  SAVE sinlon, coslon, sinlondlon, coslondlon
  REAL dyn1, dyn2, dys1, dys2
  REAL qpn, qps, dqzpn, dqzps
  REAL smn, sms, s0n, s0s, sxn(iip1), sxs(iip1)
  REAL qmin, pente_max

  REAL ssum
  INTEGER ismax, ismin, lati, latf
  EXTERNAL ssum, convflu, ismin, ismax
  LOGICAL first
  SAVE first
  ! fin modif

  ! EXTERNAL masskg
  EXTERNAL advx
  EXTERNAL advy
  EXTERNAL advz

  ! modif Fred 24 03 96
  DATA first/.TRUE./

  limit = .TRUE.
  pente_max = 2
  ! if (mode.eq.1.or.mode.eq.3) then
  ! if (mode.eq.1) then
  IF (mode>=1) THEN
    lati = 2
    latf = jjm
  ELSE
    lati = 1
    latf = jjp1
  END IF

  qmin = 0.4995
  qmin = 0.
  IF (first) THEN
    PRINT *, 'SCHEMA AMONT NOUVEAU'
    first = .FALSE.
    DO i = 2, iip1
      coslon(i) = cos(rlonv(i))
      sinlon(i) = sin(rlonv(i))
      coslondlon(i) = coslon(i)*(rlonu(i)-rlonu(i-1))/pi
      sinlondlon(i) = sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
      PRINT *, coslondlon(i), sinlondlon(i)
    END DO
    coslon(1) = coslon(iip1)
    coslondlon(1) = coslondlon(iip1)
    sinlon(1) = sinlon(iip1)
    sinlondlon(1) = sinlondlon(iip1)
    PRINT *, 'sum sinlondlon ', ssum(iim, sinlondlon, 1)/sinlondlon(1)
    PRINT *, 'sum coslondlon ', ssum(iim, coslondlon, 1)/coslondlon(1)
    DO l = 1, llm
      DO j = 1, jjp1
        DO i = 1, iip1
          q(i, j, l, 1) = 0.
          q(i, j, l, 2) = 0.
          q(i, j, l, 3) = 0.
        END DO
      END DO
    END DO

  END IF

  ! *** q contient les qqtes de traceur avant l'advection

  ! *** Affectation des tableaux S a partir de Q
  ! *** Rem : utilisation de SCOPY ulterieurement

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        s0(i, j, llm+1-l) = q(i, j, l, 0)
        sx(i, j, llm+1-l) = q(i, j, l, 1)
        sy(i, j, llm+1-l) = q(i, j, l, 2)
        sz(i, j, llm+1-l) = q(i, j, l, 3)
      END DO
    END DO
  END DO

  ! *** On calcule la masse d'air en kg

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        sm(i, j, llm+1-l) = masse(i, j, l)
      END DO
    END DO
  END DO

  ! *** On converti les champs S en atome (resp. kg)
  ! *** Les routines d'advection traitent les champs
  ! *** a advecter si ces derniers sont en atome (resp. kg)
  ! *** A optimiser !!!

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        s0(i, j, l) = s0(i, j, l)*sm(i, j, l)
        sx(i, j, l) = sx(i, j, l)*sm(i, j, l)
        sy(i, j, l) = sy(i, j, l)*sm(i, j, l)
        sz(i, j, l) = sz(i, j, l)*sm(i, j, l)
      END DO
    END DO
  END DO

  ! *** Appel des subroutines d'advection en X, en Y et en Z
  ! *** Advection avec "time-splitting"

  IF (mode==2) THEN
    DO l = 1, llm
      s0s = 0.
      s0n = 0.
      dyn1 = 0.
      dys1 = 0.
      dyn2 = 0.
      dys2 = 0.
      smn = 0.
      sms = 0.
      DO i = 1, iim
        smn = smn + sm(i, 1, l)
        sms = sms + sm(i, jjp1, l)
        s0n = s0n + s0(i, 1, l)
        s0s = s0s + s0(i, jjp1, l)
        zz = sy(i, 1, l)/sm(i, 1, l)
        dyn1 = dyn1 + sinlondlon(i)*zz
        dyn2 = dyn2 + coslondlon(i)*zz
        zz = sy(i, jjp1, l)/sm(i, jjp1, l)
        dys1 = dys1 + sinlondlon(i)*zz
        dys2 = dys2 + coslondlon(i)*zz
      END DO
      DO i = 1, iim
        sy(i, 1, l) = dyn1*sinlon(i) + dyn2*coslon(i)
        sy(i, jjp1, l) = dys1*sinlon(i) + dys2*coslon(i)
      END DO
      DO i = 1, iim
        s0(i, 1, l) = s0n/smn + sy(i, 1, l)
        s0(i, jjp1, l) = s0s/sms - sy(i, jjp1, l)
      END DO

      s0(iip1, 1, l) = s0(1, 1, l)
      s0(iip1, jjp1, l) = s0(1, jjp1, l)

      DO i = 1, iim
        sxn(i) = s0(i+1, 1, l) - s0(i, 1, l)
        sxs(i) = s0(i+1, jjp1, l) - s0(i, jjp1, l)
        ! on rerentre les masses
      END DO
      DO i = 1, iim
        sy(i, 1, l) = sy(i, 1, l)*sm(i, 1, l)
        sy(i, jjp1, l) = sy(i, jjp1, l)*sm(i, jjp1, l)
        s0(i, 1, l) = s0(i, 1, l)*sm(i, 1, l)
        s0(i, jjp1, l) = s0(i, jjp1, l)*sm(i, jjp1, l)
      END DO
      sxn(iip1) = sxn(1)
      sxs(iip1) = sxs(1)
      DO i = 1, iim
        sx(i+1, 1, l) = 0.25*(sxn(i)+sxn(i+1))*sm(i+1, 1, l)
        sx(i+1, jjp1, l) = 0.25*(sxs(i)+sxs(i+1))*sm(i+1, jjp1, l)
      END DO
      s0(iip1, 1, l) = s0(1, 1, l)
      s0(iip1, jjp1, l) = s0(1, jjp1, l)
      sy(iip1, 1, l) = sy(1, 1, l)
      sy(iip1, jjp1, l) = sy(1, jjp1, l)
      sx(1, 1, l) = sx(iip1, 1, l)
      sx(1, jjp1, l) = sx(iip1, jjp1, l)
    END DO
  END IF

  IF (mode==4) THEN
    DO l = 1, llm
      DO i = 1, iip1
        sx(i, 1, l) = 0.
        sx(i, jjp1, l) = 0.
        sy(i, 1, l) = 0.
        sy(i, jjp1, l) = 0.
      END DO
    END DO
  END IF
  CALL limx(s0, sx, sm, pente_max)
  CALL advx(limit, .5*dtvr, pbaru, sm, s0, sx, sy, sz, lati, latf)
  IF (mode==4) THEN
    DO l = 1, llm
      DO i = 1, iip1
        sx(i, 1, l) = 0.
        sx(i, jjp1, l) = 0.
        sy(i, 1, l) = 0.
        sy(i, jjp1, l) = 0.
      END DO
    END DO
  END IF
  CALL limy(s0, sy, sm, pente_max)
  CALL advy(limit, .5*dtvr, pbarv, sm, s0, sx, sy, sz)
  DO j = 1, jjp1
    DO i = 1, iip1
      sz(i, j, 1) = 0.
      sz(i, j, llm) = 0.
    END DO
  END DO
  CALL limz(s0, sz, sm, pente_max)
  CALL advz(limit, dtvr, w, sm, s0, sx, sy, sz)
  IF (mode==4) THEN
    DO l = 1, llm
      DO i = 1, iip1
        sx(i, 1, l) = 0.
        sx(i, jjp1, l) = 0.
        sy(i, 1, l) = 0.
        sy(i, jjp1, l) = 0.
      END DO
    END DO
  END IF
  CALL limy(s0, sy, sm, pente_max)
  CALL advy(limit, .5*dtvr, pbarv, sm, s0, sx, sy, sz)
  DO l = 1, llm
    DO j = 1, jjp1
      sm(iip1, j, l) = sm(1, j, l)
      s0(iip1, j, l) = s0(1, j, l)
      sx(iip1, j, l) = sx(1, j, l)
      sy(iip1, j, l) = sy(1, j, l)
      sz(iip1, j, l) = sz(1, j, l)
    END DO
  END DO


  IF (mode==4) THEN
    DO l = 1, llm
      DO i = 1, iip1
        sx(i, 1, l) = 0.
        sx(i, jjp1, l) = 0.
        sy(i, 1, l) = 0.
        sy(i, jjp1, l) = 0.
      END DO
    END DO
  END IF
  CALL limx(s0, sx, sm, pente_max)
  CALL advx(limit, .5*dtvr, pbaru, sm, s0, sx, sy, sz, lati, latf)
  ! ***   On repasse les S dans la variable q directement 14/10/94
  ! On revient a des rapports de melange en divisant par la masse

  ! En dehors des poles:

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        q(i, j, llm+1-l, 0) = s0(i, j, l)/sm(i, j, l)
        q(i, j, llm+1-l, 1) = sx(i, j, l)/sm(i, j, l)
        q(i, j, llm+1-l, 2) = sy(i, j, l)/sm(i, j, l)
        q(i, j, llm+1-l, 3) = sz(i, j, l)/sm(i, j, l)
      END DO
    END DO
  END DO

  ! Traitements specifiques au pole

  IF (mode>=1) THEN
    DO l = 1, llm
      ! filtrages aux poles
      masn = ssum(iim, sm(1,1,l), 1)
      mass = ssum(iim, sm(1,jjp1,l), 1)
      qpn = ssum(iim, s0(1,1,l), 1)/masn
      qps = ssum(iim, s0(1,jjp1,l), 1)/mass
      dqzpn = ssum(iim, sz(1,1,l), 1)/masn
      dqzps = ssum(iim, sz(1,jjp1,l), 1)/mass
      DO i = 1, iip1
        q(i, 1, llm+1-l, 3) = dqzpn
        q(i, jjp1, llm+1-l, 3) = dqzps
        q(i, 1, llm+1-l, 0) = qpn
        q(i, jjp1, llm+1-l, 0) = qps
      END DO
      IF (mode==3) THEN
        dyn1 = 0.
        dys1 = 0.
        dyn2 = 0.
        dys2 = 0.
        DO i = 1, iim
          dyn1 = dyn1 + sinlondlon(i)*sy(i, 1, l)/sm(i, 1, l)
          dyn2 = dyn2 + coslondlon(i)*sy(i, 1, l)/sm(i, 1, l)
          dys1 = dys1 + sinlondlon(i)*sy(i, jjp1, l)/sm(i, jjp1, l)
          dys2 = dys2 + coslondlon(i)*sy(i, jjp1, l)/sm(i, jjp1, l)
        END DO
        DO i = 1, iim
          q(i, 1, llm+1-l, 2) = (sinlon(i)*dyn1+coslon(i)*dyn2)
          q(i, 1, llm+1-l, 0) = q(i, 1, llm+1-l, 0) + q(i, 1, llm+1-l, 2)
          q(i, jjp1, llm+1-l, 2) = (sinlon(i)*dys1+coslon(i)*dys2)
          q(i, jjp1, llm+1-l, 0) = q(i, jjp1, llm+1-l, 0) - &
            q(i, jjp1, llm+1-l, 2)
        END DO
      END IF
      IF (mode==1) THEN
        ! on filtre les valeurs au bord de la "grande maille pole"
        dyn1 = 0.
        dys1 = 0.
        dyn2 = 0.
        dys2 = 0.
        DO i = 1, iim
          zz = s0(i, 2, l)/sm(i, 2, l) - q(i, 1, llm+1-l, 0)
          dyn1 = dyn1 + sinlondlon(i)*zz
          dyn2 = dyn2 + coslondlon(i)*zz
          zz = q(i, jjp1, llm+1-l, 0) - s0(i, jjm, l)/sm(i, jjm, l)
          dys1 = dys1 + sinlondlon(i)*zz
          dys2 = dys2 + coslondlon(i)*zz
        END DO
        DO i = 1, iim
          q(i, 1, llm+1-l, 2) = (sinlon(i)*dyn1+coslon(i)*dyn2)/2.
          q(i, 1, llm+1-l, 0) = q(i, 1, llm+1-l, 0) + q(i, 1, llm+1-l, 2)
          q(i, jjp1, llm+1-l, 2) = (sinlon(i)*dys1+coslon(i)*dys2)/2.
          q(i, jjp1, llm+1-l, 0) = q(i, jjp1, llm+1-l, 0) - &
            q(i, jjp1, llm+1-l, 2)
        END DO
        q(iip1, 1, llm+1-l, 0) = q(1, 1, llm+1-l, 0)
        q(iip1, jjp1, llm+1-l, 0) = q(1, jjp1, llm+1-l, 0)

        DO i = 1, iim
          sxn(i) = q(i+1, 1, llm+1-l, 0) - q(i, 1, llm+1-l, 0)
          sxs(i) = q(i+1, jjp1, llm+1-l, 0) - q(i, jjp1, llm+1-l, 0)
        END DO
        sxn(iip1) = sxn(1)
        sxs(iip1) = sxs(1)
        DO i = 1, iim
          q(i+1, 1, llm+1-l, 1) = 0.25*(sxn(i)+sxn(i+1))
          q(i+1, jjp1, llm+1-l, 1) = 0.25*(sxs(i)+sxs(i+1))
        END DO
        q(1, 1, llm+1-l, 1) = q(iip1, 1, llm+1-l, 1)
        q(1, jjp1, llm+1-l, 1) = q(iip1, jjp1, llm+1-l, 1)

      END IF

    END DO
  END IF

  ! bouclage en longitude
  DO iq = 0, 3
    DO l = 1, llm
      DO j = 1, jjp1
        q(iip1, j, l, iq) = q(1, j, l, iq)
      END DO
    END DO
  END DO

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        IF (q(i,j,l,0)<0.) THEN
          q(i, j, l, 0) = 0.
        END IF
      END DO
    END DO
  END DO

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        IF (q(i,j,l,0)<qmin) PRINT *, 'apres pentes, s0(', i, ',', j, ',', l, &
          ')=', q(i, j, l, 0)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE pentes_ini
