
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advn.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE advn(q, masse, w, pbaru, pbarv, pdt, mode)

  ! Auteur : F. Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....

  ! pbaru,pbarv,w flux de masse en u ,v ,w
  ! pdt pas de temps

  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  USE comgeom
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  INTEGER mode
  REAL masse(ip1jmp1, llm)
  REAL, INTENT (IN) :: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)
  REAL q(ip1jmp1, llm)
  REAL w(ip1jmp1, llm), pdt

  ! Local
  ! ---------

  INTEGER i, ij, l, j, ii
  INTEGER ijlqmin, iqmin, jqmin, lqmin
  INTEGER ismin

  REAL zm(ip1jmp1, llm), newmasse
  REAL mu(ip1jmp1, llm)
  REAL mv(ip1jm, llm)
  REAL mw(ip1jmp1, llm+1)
  REAL zq(ip1jmp1, llm), zz, qpn, qps
  REAL zqg(ip1jmp1, llm), zqd(ip1jmp1, llm)
  REAL zqs(ip1jmp1, llm), zqn(ip1jmp1, llm)
  REAL zqh(ip1jmp1, llm), zqb(ip1jmp1, llm)
  REAL temps0, temps1, temps2, temps3
  REAL ztemps1, ztemps2, ztemps3, ssum
  LOGICAL testcpu
  SAVE testcpu
  SAVE temps1, temps2, temps3
  REAL zzpbar, zzw

  REAL qmin, qmax
  DATA qmin, qmax/0., 1./
  DATA testcpu/.FALSE./
  DATA temps1, temps2, temps3/0., 0., 0./

  zzpbar = 0.5*pdt
  zzw = pdt

  DO l = 1, llm
    DO ij = iip2, ip1jm
      mu(ij, l) = pbaru(ij, l)*zzpbar
    END DO
    DO ij = 1, ip1jm
      mv(ij, l) = pbarv(ij, l)*zzpbar
    END DO
    DO ij = 1, ip1jmp1
      mw(ij, l) = w(ij, l)*zzw
    END DO
  END DO

  DO ij = 1, ip1jmp1
    mw(ij, llm+1) = 0.
  END DO

  DO l = 1, llm
    qpn = 0.
    qps = 0.
    DO ij = 1, iim
      qpn = qpn + q(ij, l)*masse(ij, l)
      qps = qps + q(ip1jm+ij, l)*masse(ip1jm+ij, l)
    END DO
    qpn = qpn/ssum(iim, masse(1,l), 1)
    qps = qps/ssum(iim, masse(ip1jm+1,l), 1)
    DO ij = 1, iip1
      q(ij, l) = qpn
      q(ip1jm+ij, l) = qps
    END DO
  END DO

  DO ij = 1, ip1jmp1
    mw(ij, llm+1) = 0.
  END DO
  DO l = 1, llm
    DO ij = 1, ip1jmp1
      zq(ij, l) = q(ij, l)
      zm(ij, l) = masse(ij, l)
    END DO
  END DO

  ! call minmaxq(zq,qmin,qmax,'avant vlx     ')
  CALL advnqx(zq, zqg, zqd)
  CALL advnx(zq, zqg, zqd, zm, mu, mode)
  CALL advnqy(zq, zqs, zqn)
  CALL advny(zq, zqs, zqn, zm, mv)
  CALL advnqz(zq, zqh, zqb)
  CALL advnz(zq, zqh, zqb, zm, mw)
  ! call vlz(zq,0.,zm,mw)
  CALL advnqy(zq, zqs, zqn)
  CALL advny(zq, zqs, zqn, zm, mv)
  CALL advnqx(zq, zqg, zqd)
  CALL advnx(zq, zqg, zqd, zm, mu, mode)
  ! call minmaxq(zq,qmin,qmax,'apres vlx     ')

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      q(ij, l) = zq(ij, l)
    END DO
    DO ij = 1, ip1jm + 1, iip1
      q(ij+iim, l) = q(ij, l)
    END DO
  END DO

  RETURN
END SUBROUTINE advn

SUBROUTINE advnqx(q, qg, qd)

  ! Auteurs:   Calcul des valeurs de q aux point u.

  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL q(ip1jmp1, llm), qg(ip1jmp1, llm), qd(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL dxqu(ip1jmp1), zqu(ip1jmp1)
  REAL zqmax(ip1jmp1), zqmin(ip1jmp1)
  LOGICAL extremum(ip1jmp1)

  INTEGER mode
  SAVE mode
  DATA mode/1/

  ! calcul des pentes en u:
  ! -----------------------
  IF (mode==0) THEN
    DO l = 1, llm
      DO ij = 1, ip1jm
        qd(ij, l) = q(ij, l)
        qg(ij, l) = q(ij, l)
      END DO
    END DO
  ELSE
    DO l = 1, llm
      DO ij = iip2, ip1jm - 1
        dxqu(ij) = q(ij+1, l) - q(ij, l)
        zqu(ij) = 0.5*(q(ij+1,l)+q(ij,l))
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        dxqu(ij) = dxqu(ij-iim)
        zqu(ij) = zqu(ij-iim)
      END DO
      DO ij = iip2, ip1jm - 1
        zqu(ij) = zqu(ij) - dxqu(ij+1)/12.
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        zqu(ij) = zqu(ij-iim)
      END DO
      DO ij = iip2 + 1, ip1jm
        zqu(ij) = zqu(ij) + dxqu(ij-1)/12.
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        zqu(ij-iim) = zqu(ij)
      END DO

      ! calcul des valeurs max et min acceptees aux interfaces

      DO ij = iip2, ip1jm - 1
        zqmax(ij) = max(q(ij+1,l), q(ij,l))
        zqmin(ij) = min(q(ij+1,l), q(ij,l))
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        zqmax(ij) = zqmax(ij-iim)
        zqmin(ij) = zqmin(ij-iim)
      END DO
      DO ij = iip2 + 1, ip1jm
        extremum(ij) = dxqu(ij)*dxqu(ij-1) <= 0.
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        extremum(ij-iim) = extremum(ij)
      END DO
      DO ij = iip2, ip1jm
        zqu(ij) = min(max(zqmin(ij),zqu(ij)), zqmax(ij))
      END DO
      DO ij = iip2 + 1, ip1jm
        IF (extremum(ij)) THEN
          qg(ij, l) = q(ij, l)
          qd(ij, l) = q(ij, l)
        ELSE
          qd(ij, l) = zqu(ij)
          qg(ij, l) = zqu(ij-1)
        END IF
      END DO
      DO ij = iip1 + iip1, ip1jm, iip1
        qd(ij-iim, l) = qd(ij, l)
        qg(ij-iim, l) = qg(ij, l)
      END DO

      GO TO 8888

      DO ij = iip2 + 1, ip1jm
        IF (extremum(ij) .AND. .NOT. extremum(ij-1)) qd(ij-1, l) = q(ij, l)
      END DO

      DO ij = iip1 + iip1, ip1jm, iip1
        qd(ij-iim, l) = qd(ij, l)
      END DO
      DO ij = iip2, ip1jm - 1
        IF (extremum(ij) .AND. .NOT. extremum(ij+1)) qg(ij+1, l) = q(ij, l)
      END DO

      DO ij = iip1 + iip1, ip1jm, iip1
        qg(ij, l) = qg(ij-iim, l)
      END DO
8888  CONTINUE
    END DO
  END IF
  RETURN
END SUBROUTINE advnqx
SUBROUTINE advnqy(q, qs, qn)

  ! Auteurs:   Calcul des valeurs de q aux point v.

  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL q(ip1jmp1, llm), qs(ip1jmp1, llm), qn(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL dyqv(ip1jm), zqv(ip1jm, llm)
  REAL zqmax(ip1jm), zqmin(ip1jm)
  LOGICAL extremum(ip1jmp1)

  INTEGER mode
  SAVE mode
  DATA mode/1/

  IF (mode==0) THEN
    DO l = 1, llm
      DO ij = 1, ip1jmp1
        qn(ij, l) = q(ij, l)
        qs(ij, l) = q(ij, l)
      END DO
    END DO
  ELSE

    ! calcul des pentes en u:
    ! -----------------------
    DO l = 1, llm
      DO ij = 1, ip1jm
        dyqv(ij) = q(ij, l) - q(ij+iip1, l)
      END DO

      DO ij = iip2, ip1jm - iip1
        zqv(ij, l) = 0.5*(q(ij+iip1,l)+q(ij,l))
        zqv(ij, l) = zqv(ij, l) + (dyqv(ij+iip1)-dyqv(ij-iip1))/12.
      END DO

      DO ij = iip2, ip1jm
        extremum(ij) = dyqv(ij)*dyqv(ij-iip1) <= 0.
      END DO

      ! Pas de pentes aux poles
      DO ij = 1, iip1
        zqv(ij, l) = q(ij, l)
        zqv(ip1jm-iip1+ij, l) = q(ip1jm+ij, l)
        extremum(ij) = .TRUE.
        extremum(ip1jmp1-iip1+ij) = .TRUE.
      END DO

      ! calcul des valeurs max et min acceptees aux interfaces
      DO ij = 1, ip1jm
        zqmax(ij) = max(q(ij+iip1,l), q(ij,l))
        zqmin(ij) = min(q(ij+iip1,l), q(ij,l))
      END DO

      DO ij = 1, ip1jm
        zqv(ij, l) = min(max(zqmin(ij),zqv(ij,l)), zqmax(ij))
      END DO

      DO ij = iip2, ip1jm
        IF (extremum(ij)) THEN
          qs(ij, l) = q(ij, l)
          qn(ij, l) = q(ij, l)
          ! if (.not.extremum(ij-iip1)) qs(ij-iip1,l)=q(ij,l)
          ! if (.not.extremum(ij+iip1)) qn(ij+iip1,l)=q(ij,l)
        ELSE
          qs(ij, l) = zqv(ij, l)
          qn(ij, l) = zqv(ij-iip1, l)
        END IF
      END DO

      DO ij = 1, iip1
        qs(ij, l) = q(ij, l)
        qn(ij, l) = q(ij, l)
        qs(ip1jm+ij, l) = q(ip1jm+ij, l)
        qn(ip1jm+ij, l) = q(ip1jm+ij, l)
      END DO

    END DO
  END IF
  RETURN
END SUBROUTINE advnqy

SUBROUTINE advnqz(q, qh, qb)

  ! Auteurs:   Calcul des valeurs de q aux point v.

  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL q(ip1jmp1, llm), qh(ip1jmp1, llm), qb(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL dzqw(ip1jmp1, llm+1), zqw(ip1jmp1, llm+1)
  REAL zqmax(ip1jmp1, llm), zqmin(ip1jmp1, llm)
  LOGICAL extremum(ip1jmp1, llm)

  INTEGER mode
  SAVE mode

  DATA mode/1/

  ! calcul des pentes en u:
  ! -----------------------

  IF (mode==0) THEN
    DO l = 1, llm
      DO ij = 1, ip1jmp1
        qb(ij, l) = q(ij, l)
        qh(ij, l) = q(ij, l)
      END DO
    END DO
  ELSE
    DO l = 2, llm
      DO ij = 1, ip1jmp1
        dzqw(ij, l) = q(ij, l-1) - q(ij, l)
        zqw(ij, l) = 0.5*(q(ij,l-1)+q(ij,l))
      END DO
    END DO
    DO ij = 1, ip1jmp1
      dzqw(ij, 1) = 0.
      dzqw(ij, llm+1) = 0.
    END DO
    DO l = 2, llm
      DO ij = 1, ip1jmp1
        zqw(ij, l) = zqw(ij, l) + (dzqw(ij,l+1)-dzqw(ij,l-1))/12.
      END DO
    END DO
    DO l = 2, llm - 1
      DO ij = 1, ip1jmp1
        extremum(ij, l) = dzqw(ij, l)*dzqw(ij, l+1) <= 0.
      END DO
    END DO

    ! Pas de pentes en bas et en haut
    DO ij = 1, ip1jmp1
      zqw(ij, 2) = q(ij, 1)
      zqw(ij, llm) = q(ij, llm)
      extremum(ij, 1) = .TRUE.
      extremum(ij, llm) = .TRUE.
    END DO

    ! calcul des valeurs max et min acceptees aux interfaces
    DO l = 2, llm
      DO ij = 1, ip1jmp1
        zqmax(ij, l) = max(q(ij,l-1), q(ij,l))
        zqmin(ij, l) = min(q(ij,l-1), q(ij,l))
      END DO
    END DO

    DO l = 2, llm
      DO ij = 1, ip1jmp1
        zqw(ij, l) = min(max(zqmin(ij,l),zqw(ij,l)), zqmax(ij,l))
      END DO
    END DO

    DO l = 2, llm - 1
      DO ij = 1, ip1jmp1
        IF (extremum(ij,l)) THEN
          qh(ij, l) = q(ij, l)
          qb(ij, l) = q(ij, l)
        ELSE
          qh(ij, l) = zqw(ij, l+1)
          qb(ij, l) = zqw(ij, l)
        END IF
      END DO
    END DO
    ! do l=2,llm-1
    ! do ij=1,ip1jmp1
    ! if(extremum(ij,l)) then
    ! if (.not.extremum(ij,l-1)) qh(ij,l-1)=q(ij,l)
    ! if (.not.extremum(ij,l+1)) qb(ij,l+1)=q(ij,l)
    ! endif
    ! enddo
    ! enddo

    DO ij = 1, ip1jmp1
      qb(ij, 1) = q(ij, 1)
      qh(ij, 1) = q(ij, 1)
      qb(ij, llm) = q(ij, llm)
      qh(ij, llm) = q(ij, llm)
    END DO

  END IF

  RETURN
END SUBROUTINE advnqz

SUBROUTINE advnx(q, qg, qd, masse, u_m, mode)

  ! Auteur : F. Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....


  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  INTEGER mode
  REAL masse(ip1jmp1, llm)
  REAL u_m(ip1jmp1, llm)
  REAL q(ip1jmp1, llm), qd(ip1jmp1, llm), qg(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER i, j, ij, l, indu(ip1jmp1), niju, iju, ijq
  INTEGER n0, nl(llm)

  REAL new_m, zu_m, zdq, zz
  REAL zsigg(ip1jmp1, llm), zsigd(ip1jmp1, llm), zsig
  REAL u_mq(ip1jmp1, llm)

  REAL zm, zq, zsigm, zsigp, zqm, zqp, zu

  LOGICAL ladvplus(ip1jmp1, llm)

  REAL prec
  SAVE prec

  DATA prec/1.E-15/

  DO l = 1, llm
    DO ij = iip2, ip1jm
      zdq = qd(ij, l) - qg(ij, l)
      IF (abs(zdq)>prec) THEN
        zsigd(ij, l) = (q(ij,l)-qg(ij,l))/zdq
        zsigg(ij, l) = 1. - zsigd(ij, l)
      ELSE
        zsigd(ij, l) = 0.5
        zsigg(ij, l) = 0.5
        qd(ij, l) = q(ij, l)
        qg(ij, l) = q(ij, l)
      END IF
    END DO
  END DO

  ! calcul de la pente maximum dans la maille en valeur absolue

  DO l = 1, llm
    DO ij = iip2, ip1jm - 1
      IF (u_m(ij,l)>=0.) THEN
        zsigp = zsigd(ij, l)
        zsigm = zsigg(ij, l)
        zqp = qd(ij, l)
        zqm = qg(ij, l)
        zm = masse(ij, l)
        zq = q(ij, l)
      ELSE
        zsigm = zsigd(ij+1, l)
        zsigp = zsigg(ij+1, l)
        zqm = qd(ij+1, l)
        zqp = qg(ij+1, l)
        zm = masse(ij+1, l)
        zq = q(ij+1, l)
      END IF
      zu = abs(u_m(ij,l))
      ladvplus(ij, l) = zu > zm
      zsig = zu/zm
      IF (zsig==0.) zsigp = 0.1
      IF (mode==1) THEN
        IF (zsig<=zsigp) THEN
          u_mq(ij, l) = u_m(ij, l)*zqp
        ELSE IF (mode==1) THEN
          u_mq(ij, l) = sign(zm, u_m(ij,l))*(zsigp*zqp+(zsig-zsigp)*zqm)
        END IF
      ELSE
        IF (zsig<=zsigp) THEN
          u_mq(ij, l) = u_m(ij, l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
        ELSE
          zz = 0.5*(zsig-zsigp)/zsigm
          u_mq(ij, l) = sign(zm, u_m(ij,l))*(0.5*(zq+zqp)*zsigp+(zsig-zsigp)* &
            (zq+zz*(zqm-zq)))
        END IF
      END IF
    END DO
  END DO

  DO l = 1, llm
    DO ij = iip1 + iip1, ip1jm, iip1
      u_mq(ij, l) = u_mq(ij-iim, l)
      ladvplus(ij, l) = ladvplus(ij-iim, l)
    END DO
  END DO

  ! =================================================================
  ! SCHEMA SEMI-LAGRAGIEN EN X DANS LES REGIONS POLAIRES
  ! =================================================================
  ! tris des regions a traiter
  n0 = 0
  DO l = 1, llm
    nl(l) = 0
    DO ij = iip2, ip1jm
      IF (ladvplus(ij,l)) THEN
        nl(l) = nl(l) + 1
        u_mq(ij, l) = 0.
      END IF
    END DO
    n0 = n0 + nl(l)
  END DO

  IF (n0>1) THEN
    IF (prt_level>9) PRINT *, &
      'Nombre de points pour lesquels on advect plus que le', &
      'contenu de la maille : ', n0

    DO l = 1, llm
      IF (nl(l)>0) THEN
        iju = 0
        ! indicage des mailles concernees par le traitement special
        DO ij = iip2, ip1jm
          IF (ladvplus(ij,l) .AND. mod(ij,iip1)/=0) THEN
            iju = iju + 1
            indu(iju) = ij
          END IF
        END DO
        niju = iju

        ! traitement des mailles
        DO iju = 1, niju
          ij = indu(iju)
          j = (ij-1)/iip1 + 1
          zu_m = u_m(ij, l)
          u_mq(ij, l) = 0.
          IF (zu_m>0.) THEN
            ijq = ij
            i = ijq - (j-1)*iip1
            ! accumulation pour les mailles completements advectees
            DO WHILE (zu_m>masse(ijq,l))
              u_mq(ij, l) = u_mq(ij, l) + q(ijq, l)*masse(ijq, l)
              zu_m = zu_m - masse(ijq, l)
              i = mod(i-2+iim, iim) + 1
              ijq = (j-1)*iip1 + i
            END DO
            ! MODIFS SPECIFIQUES DU SCHEMA
            ! ajout de la maille non completement advectee
            zsig = zu_m/masse(ijq, l)
            IF (zsig<=zsigd(ijq,l)) THEN
              u_mq(ij, l) = u_mq(ij, l) + zu_m*(qd(ijq,l)-0.5*zsig/zsigd(ijq, &
                l)*(qd(ijq,l)-q(ijq,l)))
            ELSE
              ! u_mq(ij,l)=u_mq(ij,l)+zu_m*q(ijq,l)
              ! goto 8888
              zz = 0.5*(zsig-zsigd(ijq,l))/zsigg(ijq, l)
              IF (.NOT. (zz>0. .AND. zz<=0.5)) THEN
                PRINT *, 'probleme2 au point ij=', ij, '  l=', l
                PRINT *, 'zz=', zz
                STOP
              END IF
              u_mq(ij, l) = u_mq(ij, l) + masse(ijq, l)*(0.5*(q(ijq, &
                l)+qd(ijq,l))*zsigd(ijq,l)+(zsig-zsigd(ijq,l))*(q(ijq, &
                l)+zz*(qg(ijq,l)-q(ijq,l))))
            END IF
          ELSE
            ijq = ij + 1
            i = ijq - (j-1)*iip1
            ! accumulation pour les mailles completements advectees
            DO WHILE (-zu_m>masse(ijq,l))
              u_mq(ij, l) = u_mq(ij, l) - q(ijq, l)*masse(ijq, l)
              zu_m = zu_m + masse(ijq, l)
              i = mod(i, iim) + 1
              ijq = (j-1)*iip1 + i
            END DO
            ! ajout de la maille non completement advectee
            ! 2eme MODIF SPECIFIQUE
            zsig = -zu_m/masse(ij+1, l)
            IF (zsig<=zsigg(ijq,l)) THEN
              u_mq(ij, l) = u_mq(ij, l) + zu_m*(qg(ijq,l)-0.5*zsig/zsigg(ijq, &
                l)*(qg(ijq,l)-q(ijq,l)))
            ELSE
              ! u_mq(ij,l)=u_mq(ij,l)+zu_m*q(ijq,l)
              ! goto 9999
              zz = 0.5*(zsig-zsigg(ijq,l))/zsigd(ijq, l)
              IF (.NOT. (zz>0. .AND. zz<=0.5)) THEN
                PRINT *, 'probleme22 au point ij=', ij, '  l=', l
                PRINT *, 'zz=', zz
                STOP
              END IF
              u_mq(ij, l) = u_mq(ij, l) - masse(ijq, l)*(0.5*(q(ijq, &
                l)+qg(ijq,l))*zsigg(ijq,l)+(zsig-zsigg(ijq,l))*(q(ijq, &
                l)+zz*(qd(ijq,l)-q(ijq,l))))
            END IF
            ! fin de la modif
          END IF
        END DO
      END IF
    END DO
  END IF ! n0.gt.0

  ! bouclage en latitude
  DO l = 1, llm
    DO ij = iip1 + iip1, ip1jm, iip1
      u_mq(ij, l) = u_mq(ij-iim, l)
    END DO
  END DO

  ! =================================================================
  ! CALCUL DE LA CONVERGENCE DES FLUX
  ! =================================================================

  DO l = 1, llm
    DO ij = iip2 + 1, ip1jm
      new_m = masse(ij, l) + u_m(ij-1, l) - u_m(ij, l)
      q(ij, l) = (q(ij,l)*masse(ij,l)+u_mq(ij-1,l)-u_mq(ij,l))/new_m
      masse(ij, l) = new_m
    END DO
    ! Modif Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
    DO ij = iip1 + iip1, ip1jm, iip1
      q(ij-iim, l) = q(ij, l)
      masse(ij-iim, l) = masse(ij, l)
    END DO
  END DO

  RETURN
END SUBROUTINE advnx
SUBROUTINE advny(q, qs, qn, masse, v_m)

  ! Auteur : F. Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! ********************************************************************
  ! nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....


  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE comgeom
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL masse(ip1jmp1, llm)
  REAL v_m(ip1jm, llm)
  REAL q(ip1jmp1, llm), qn(ip1jmp1, llm), qs(ip1jmp1, llm)

  ! Local
  ! ---------

  INTEGER ij, l

  REAL new_m, zdq, zz
  REAL zsigs(ip1jmp1), zsign(ip1jmp1), zsig
  REAL v_mq(ip1jm, llm)
  REAL convpn, convps, convmpn, convmps, massen, masses
  REAL zm, zq, zsigm, zsigp, zqm, zqp
  REAL ssum
  REAL prec
  SAVE prec

  DATA prec/1.E-15/

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      zdq = qn(ij, l) - qs(ij, l)
      IF (abs(zdq)>prec) THEN
        zsign(ij) = (q(ij,l)-qs(ij,l))/zdq
        zsigs(ij) = 1. - zsign(ij)
      ELSE
        zsign(ij) = 0.5
        zsigs(ij) = 0.5
      END IF
    END DO

    ! calcul de la pente maximum dans la maille en valeur absolue

    DO ij = 1, ip1jm
      IF (v_m(ij,l)>=0.) THEN
        zsigp = zsign(ij+iip1)
        zsigm = zsigs(ij+iip1)
        zqp = qn(ij+iip1, l)
        zqm = qs(ij+iip1, l)
        zm = masse(ij+iip1, l)
        zq = q(ij+iip1, l)
      ELSE
        zsigm = zsign(ij)
        zsigp = zsigs(ij)
        zqm = qn(ij, l)
        zqp = qs(ij, l)
        zm = masse(ij, l)
        zq = q(ij, l)
      END IF
      zsig = abs(v_m(ij,l))/zm
      IF (zsig==0.) zsigp = 0.1
      IF (zsig<=zsigp) THEN
        v_mq(ij, l) = v_m(ij, l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
      ELSE
        zz = 0.5*(zsig-zsigp)/zsigm
        v_mq(ij, l) = sign(zm, v_m(ij,l))*(0.5*(zq+zqp)*zsigp+(zsig-zsigp)*( &
          zq+zz*(zqm-zq)))
      END IF
    END DO
  END DO

  DO l = 1, llm
    DO ij = iip2, ip1jm
      new_m = masse(ij, l) + v_m(ij, l) - v_m(ij-iip1, l)
      q(ij, l) = (q(ij,l)*masse(ij,l)+v_mq(ij,l)-v_mq(ij-iip1,l))/new_m
      masse(ij, l) = new_m
    END DO
    ! .-. ancienne version
    convpn = ssum(iim, v_mq(1,l), 1)
    convmpn = ssum(iim, v_m(1,l), 1)
    massen = ssum(iim, masse(1,l), 1)
    new_m = massen + convmpn
    q(1, l) = (q(1,l)*massen+convpn)/new_m
    DO ij = 1, iip1
      q(ij, l) = q(1, l)
      masse(ij, l) = new_m*aire(ij)/apoln
    END DO

    convps = -ssum(iim, v_mq(ip1jm-iim,l), 1)
    convmps = -ssum(iim, v_m(ip1jm-iim,l), 1)
    masses = ssum(iim, masse(ip1jm+1,l), 1)
    new_m = masses + convmps
    q(ip1jm+1, l) = (q(ip1jm+1,l)*masses+convps)/new_m
    DO ij = ip1jm + 1, ip1jmp1
      q(ij, l) = q(ip1jm+1, l)
      masse(ij, l) = new_m*aire(ij)/apols
    END DO
  END DO

  RETURN
END SUBROUTINE advny
SUBROUTINE advnz(q, qh, qb, masse, w_m)

  ! Auteurs:   F.Hourdin

  ! ********************************************************************
  ! Shema  d'advection " pseudo amont " .
  ! b designe le bas et h le haut
  ! il y a une correspondance entre le b en z et le d en x
  ! ********************************************************************


  ! --------------------------------------------------------------------
  USE dimens_m
  USE paramet_m
  USE comgeom
  USE conf_gcm_m
  IMPLICIT NONE



  ! Arguments:
  ! ----------
  REAL masse(ip1jmp1, llm)
  REAL w_m(ip1jmp1, llm+1)
  REAL q(ip1jmp1, llm), qb(ip1jmp1, llm), qh(ip1jmp1, llm)


  ! Local
  ! ---------

  INTEGER ij, l

  REAL new_m, zdq, zz
  REAL zsigh(ip1jmp1, llm), zsigb(ip1jmp1, llm), zsig
  REAL w_mq(ip1jmp1, llm+1)
  REAL zm, zq, zsigm, zsigp, zqm, zqp
  REAL prec
  SAVE prec

  DATA prec/1.E-13/

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      zdq = qb(ij, l) - qh(ij, l)
      IF (abs(zdq)>prec) THEN
        zsigb(ij, l) = (q(ij,l)-qh(ij,l))/zdq
        zsigh(ij, l) = 1. - zsigb(ij, l)
        zsigb(ij, l) = min(max(zsigb(ij,l),0.), 1.)
      ELSE
        zsigb(ij, l) = 0.5
        zsigh(ij, l) = 0.5
      END IF
    END DO
  END DO

  ! calcul de la pente maximum dans la maille en valeur absolue
  DO l = 2, llm
    DO ij = 1, ip1jmp1
      IF (w_m(ij,l)>=0.) THEN
        zsigp = zsigb(ij, l)
        zsigm = zsigh(ij, l)
        zqp = qb(ij, l)
        zqm = qh(ij, l)
        zm = masse(ij, l)
        zq = q(ij, l)
      ELSE
        zsigm = zsigb(ij, l-1)
        zsigp = zsigh(ij, l-1)
        zqm = qb(ij, l-1)
        zqp = qh(ij, l-1)
        zm = masse(ij, l-1)
        zq = q(ij, l-1)
      END IF
      zsig = abs(w_m(ij,l))/zm
      IF (zsig==0.) zsigp = 0.1
      IF (zsig<=zsigp) THEN
        w_mq(ij, l) = w_m(ij, l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
      ELSE
        zz = 0.5*(zsig-zsigp)/zsigm
        w_mq(ij, l) = sign(zm, w_m(ij,l))*(0.5*(zq+zqp)*zsigp+(zsig-zsigp)*( &
          zq+zz*(zqm-zq)))
      END IF
    END DO
  END DO

  DO ij = 1, ip1jmp1
    w_mq(ij, llm+1) = 0.
    w_mq(ij, 1) = 0.
  END DO

  DO l = 1, llm
    DO ij = 1, ip1jmp1
      new_m = masse(ij, l) + w_m(ij, l+1) - w_m(ij, l)
      q(ij, l) = (q(ij,l)*masse(ij,l)+w_mq(ij,l+1)-w_mq(ij,l))/new_m
      masse(ij, l) = new_m
    END DO
  END DO

END SUBROUTINE advnz
