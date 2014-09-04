
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advzp.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE advzp(limit, dtz, w, sm, s0, ssx, sy, sz, ssxx, ssxy, ssxz, syy, &
    syz, szz, ntra)

  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE comgeom
  IMPLICIT NONE

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! second-order moments (SOM) advection of tracer in Z direction  C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! Source : Pascal Simon ( Meteo, CNRM )                          C
  ! Adaptation : A.A. (LGGE)                                       C
  ! Derniere Modif : 19/11/95 LAST                                 C
  ! C
  ! sont les arguments d'entree pour le s-pg                       C
  ! C
  ! argument de sortie du s-pg                                     C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Rem : Probleme aux poles il faut reecrire ce cas specifique
  ! Attention au sens de l'indexation



  ! parametres principaux du modele


  ! Arguments :
  ! ----------
  ! dty : frequence fictive d'appel du transport
  ! parbu,pbarv : flux de masse en x et y en Pa.m2.s-1

  INTEGER lon, lat, niv
  INTEGER i, j, jv, k, l, lp
  INTEGER ntra
  ! PARAMETER (ntra = 1)

  REAL dtz
  REAL w(iip1, jjp1, llm)

  ! moments: SM  total mass in each grid box
  ! S0  mass of tracer in each grid box
  ! Si  1rst order moment in i direction

  REAL sm(iip1, jjp1, llm), s0(iip1, jjp1, llm, ntra)
  REAL ssx(iip1, jjp1, llm, ntra), sy(iip1, jjp1, llm, ntra), &
    sz(iip1, jjp1, llm, ntra), ssxx(iip1, jjp1, llm, ntra), &
    ssxy(iip1, jjp1, llm, ntra), ssxz(iip1, jjp1, llm, ntra), &
    syy(iip1, jjp1, llm, ntra), syz(iip1, jjp1, llm, ntra), &
    szz(iip1, jjp1, llm, ntra)

  ! Local :
  ! -------

  ! mass fluxes across the boundaries (UGRI,VGRI,WGRI)
  ! mass fluxes in kg
  ! declaration :

  REAL wgri(iip1, jjp1, 0:llm)

  ! Rem : UGRI et VGRI ne sont pas utilises dans
  ! cette subroutine ( advection en z uniquement )
  ! Rem 2 :le dimensionnement de VGRI depend de celui de pbarv
  ! attention a celui de WGRI

  ! the moments F are similarly defined and used as temporary
  ! storage for portions of the grid boxes in transit

  ! the moments Fij are used as temporary storage for
  ! portions of the grid boxes in transit at the current level

  ! work arrays


  REAL f0(iim, llm, ntra), fm(iim, llm)
  REAL fx(iim, llm, ntra), fy(iim, llm, ntra)
  REAL fz(iim, llm, ntra)
  REAL fxx(iim, llm, ntra), fxy(iim, llm, ntra)
  REAL fxz(iim, llm, ntra), fyy(iim, llm, ntra)
  REAL fyz(iim, llm, ntra), fzz(iim, llm, ntra)

  ! work arrays

  REAL alf(iim), alf1(iim)
  REAL alfq(iim), alf1q(iim)
  REAL alf2(iim), alf3(iim)
  REAL alf4(iim)
  REAL temptm ! Just temporal variable
  REAL slpmax, s1max, s1new, s2new

  REAL sqi, sqf
  LOGICAL limit

  lon = iim ! rem : Il est possible qu'un pbl. arrive ici
  lat = jjp1 ! a cause des dim. differentes entre les
  niv = llm !       tab. S et VGRI

  ! -----------------------------------------------------------------
  ! *** Test : diag de la qtite totale de traceur dans
  ! l'atmosphere avant l'advection en Y

  sqi = 0.
  sqf = 0.

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        sqi = sqi + s0(i, j, l, ntra)
      END DO
    END DO
  END DO
  PRINT *, '---------- DIAG DANS ADVZP - ENTREE --------'
  PRINT *, 'sqi=', sqi

  ! -----------------------------------------------------------------
  ! Interface : adaptation nouveau modele
  ! -------------------------------------

  ! Conversion des flux de masses en kg

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        wgri(i, j, llm+1-l) = w(i, j, l)
      END DO
    END DO
  END DO
  DO j = 1, jjp1
    DO i = 1, iip1
      wgri(i, j, 0) = 0.
    END DO
  END DO

  ! AA rem : Je ne suis pas sur du signe
  ! AA       Je ne suis pas sur pour le 0:llm

  ! -----------------------------------------------------------------
  ! ---------------------- START HERE -------------------------------

  ! boucle sur les latitudes

  DO k = 1, lat

    ! place limits on appropriate moments before transport
    ! (if flux-limiting is to be applied)

    IF (.NOT. limit) GO TO 101

    DO jv = 1, ntra
      DO l = 1, niv
        DO i = 1, lon
          IF (s0(i,k,l,jv)>0.) THEN
            slpmax = s0(i, k, l, jv)
            s1max = 1.5*slpmax
            s1new = amin1(s1max, amax1(-s1max,sz(i,k,l,jv)))
            s2new = amin1(2.*slpmax-abs(s1new)/3., amax1(abs( &
              s1new)-slpmax,szz(i,k,l,jv)))
            sz(i, k, l, jv) = s1new
            szz(i, k, l, jv) = s2new
            ssxz(i, k, l, jv) = amin1(slpmax, amax1(-slpmax,ssxz(i,k,l,jv)))
            syz(i, k, l, jv) = amin1(slpmax, amax1(-slpmax,syz(i,k,l,jv)))
          ELSE
            sz(i, k, l, jv) = 0.
            szz(i, k, l, jv) = 0.
            ssxz(i, k, l, jv) = 0.
            syz(i, k, l, jv) = 0.
          END IF
        END DO
      END DO
    END DO

101 CONTINUE

    ! boucle sur les niveaux intercouches de 1 a NIV-1
    ! (flux nul au sommet L=0 et a la base L=NIV)

    ! calculate flux and moments between adjacent boxes
    ! (flux from LP to L if WGRI(L).lt.0, from L to LP if WGRI(L).gt.0)
    ! 1- create temporary moments/masses for partial boxes in transit
    ! 2- reajusts moments remaining in the box

    DO l = 1, niv - 1
      lp = l + 1

      DO i = 1, lon

        IF (wgri(i,k,l)<0.) THEN
          fm(i, l) = -wgri(i, k, l)*dtz
          alf(i) = fm(i, l)/sm(i, k, lp)
          sm(i, k, lp) = sm(i, k, lp) - fm(i, l)
        ELSE
          fm(i, l) = wgri(i, k, l)*dtz
          alf(i) = fm(i, l)/sm(i, k, l)
          sm(i, k, l) = sm(i, k, l) - fm(i, l)
        END IF

        alfq(i) = alf(i)*alf(i)
        alf1(i) = 1. - alf(i)
        alf1q(i) = alf1(i)*alf1(i)
        alf2(i) = alf1(i) - alf(i)
        alf3(i) = alf(i)*alfq(i)
        alf4(i) = alf1(i)*alf1q(i)

      END DO

      DO jv = 1, ntra
        DO i = 1, lon

          IF (wgri(i,k,l)<0.) THEN

            f0(i, l, jv) = alf(i)*(s0(i,k,lp,jv)-alf1(i)*(sz(i,k,lp, &
              jv)-alf2(i)*szz(i,k,lp,jv)))
            fz(i, l, jv) = alfq(i)*(sz(i,k,lp,jv)-3.*alf1(i)*szz(i,k,lp,jv))
            fzz(i, l, jv) = alf3(i)*szz(i, k, lp, jv)
            fxz(i, l, jv) = alfq(i)*ssxz(i, k, lp, jv)
            fyz(i, l, jv) = alfq(i)*syz(i, k, lp, jv)
            fx(i, l, jv) = alf(i)*(ssx(i,k,lp,jv)-alf1(i)*ssxz(i,k,lp,jv))
            fy(i, l, jv) = alf(i)*(sy(i,k,lp,jv)-alf1(i)*syz(i,k,lp,jv))
            fxx(i, l, jv) = alf(i)*ssxx(i, k, lp, jv)
            fxy(i, l, jv) = alf(i)*ssxy(i, k, lp, jv)
            fyy(i, l, jv) = alf(i)*syy(i, k, lp, jv)

            s0(i, k, lp, jv) = s0(i, k, lp, jv) - f0(i, l, jv)
            sz(i, k, lp, jv) = alf1q(i)*(sz(i,k,lp,jv)+3.*alf(i)*szz(i,k,lp, &
              jv))
            szz(i, k, lp, jv) = alf4(i)*szz(i, k, lp, jv)
            ssxz(i, k, lp, jv) = alf1q(i)*ssxz(i, k, lp, jv)
            syz(i, k, lp, jv) = alf1q(i)*syz(i, k, lp, jv)
            ssx(i, k, lp, jv) = ssx(i, k, lp, jv) - fx(i, l, jv)
            sy(i, k, lp, jv) = sy(i, k, lp, jv) - fy(i, l, jv)
            ssxx(i, k, lp, jv) = ssxx(i, k, lp, jv) - fxx(i, l, jv)
            ssxy(i, k, lp, jv) = ssxy(i, k, lp, jv) - fxy(i, l, jv)
            syy(i, k, lp, jv) = syy(i, k, lp, jv) - fyy(i, l, jv)

          ELSE

            f0(i, l, jv) = alf(i)*(s0(i,k,l,jv)+alf1(i)*(sz(i,k,l, &
              jv)+alf2(i)*szz(i,k,l,jv)))
            fz(i, l, jv) = alfq(i)*(sz(i,k,l,jv)+3.*alf1(i)*szz(i,k,l,jv))
            fzz(i, l, jv) = alf3(i)*szz(i, k, l, jv)
            fxz(i, l, jv) = alfq(i)*ssxz(i, k, l, jv)
            fyz(i, l, jv) = alfq(i)*syz(i, k, l, jv)
            fx(i, l, jv) = alf(i)*(ssx(i,k,l,jv)+alf1(i)*ssxz(i,k,l,jv))
            fy(i, l, jv) = alf(i)*(sy(i,k,l,jv)+alf1(i)*syz(i,k,l,jv))
            fxx(i, l, jv) = alf(i)*ssxx(i, k, l, jv)
            fxy(i, l, jv) = alf(i)*ssxy(i, k, l, jv)
            fyy(i, l, jv) = alf(i)*syy(i, k, l, jv)

            s0(i, k, l, jv) = s0(i, k, l, jv) - f0(i, l, jv)
            sz(i, k, l, jv) = alf1q(i)*(sz(i,k,l,jv)-3.*alf(i)*szz(i,k,l,jv))
            szz(i, k, l, jv) = alf4(i)*szz(i, k, l, jv)
            ssxz(i, k, l, jv) = alf1q(i)*ssxz(i, k, l, jv)
            syz(i, k, l, jv) = alf1q(i)*syz(i, k, l, jv)
            ssx(i, k, l, jv) = ssx(i, k, l, jv) - fx(i, l, jv)
            sy(i, k, l, jv) = sy(i, k, l, jv) - fy(i, l, jv)
            ssxx(i, k, l, jv) = ssxx(i, k, l, jv) - fxx(i, l, jv)
            ssxy(i, k, l, jv) = ssxy(i, k, l, jv) - fxy(i, l, jv)
            syy(i, k, l, jv) = syy(i, k, l, jv) - fyy(i, l, jv)

          END IF

        END DO
      END DO

    END DO

    ! puts the temporary moments Fi into appropriate neighboring boxes

    DO l = 1, niv - 1
      lp = l + 1

      DO i = 1, lon

        IF (wgri(i,k,l)<0.) THEN
          sm(i, k, l) = sm(i, k, l) + fm(i, l)
          alf(i) = fm(i, l)/sm(i, k, l)
        ELSE
          sm(i, k, lp) = sm(i, k, lp) + fm(i, l)
          alf(i) = fm(i, l)/sm(i, k, lp)
        END IF

        alf1(i) = 1. - alf(i)
        alfq(i) = alf(i)*alf(i)
        alf1q(i) = alf1(i)*alf1(i)
        alf2(i) = alf(i)*alf1(i)
        alf3(i) = alf1(i) - alf(i)

      END DO

      DO jv = 1, ntra
        DO i = 1, lon

          IF (wgri(i,k,l)<0.) THEN

            temptm = -alf(i)*s0(i, k, l, jv) + alf1(i)*f0(i, l, jv)
            s0(i, k, l, jv) = s0(i, k, l, jv) + f0(i, l, jv)
            szz(i, k, l, jv) = alfq(i)*fzz(i, l, jv) + &
              alf1q(i)*szz(i, k, l, jv) + 5.*(alf2(i)*(fz(i,l,jv)-sz(i,k,l, &
              jv))+alf3(i)*temptm)
            sz(i, k, l, jv) = alf(i)*fz(i, l, jv) + alf1(i)*sz(i, k, l, jv) + &
              3.*temptm
            ssxz(i, k, l, jv) = alf(i)*fxz(i, l, jv) + &
              alf1(i)*ssxz(i, k, l, jv) + 3.*(alf1(i)*fx(i,l,jv)-alf(i)*ssx(i &
              ,k,l,jv))
            syz(i, k, l, jv) = alf(i)*fyz(i, l, jv) + &
              alf1(i)*syz(i, k, l, jv) + 3.*(alf1(i)*fy(i,l,jv)-alf(i)*sy(i,k &
              ,l,jv))
            ssx(i, k, l, jv) = ssx(i, k, l, jv) + fx(i, l, jv)
            sy(i, k, l, jv) = sy(i, k, l, jv) + fy(i, l, jv)
            ssxx(i, k, l, jv) = ssxx(i, k, l, jv) + fxx(i, l, jv)
            ssxy(i, k, l, jv) = ssxy(i, k, l, jv) + fxy(i, l, jv)
            syy(i, k, l, jv) = syy(i, k, l, jv) + fyy(i, l, jv)

          ELSE

            temptm = alf(i)*s0(i, k, lp, jv) - alf1(i)*f0(i, l, jv)
            s0(i, k, lp, jv) = s0(i, k, lp, jv) + f0(i, l, jv)
            szz(i, k, lp, jv) = alfq(i)*fzz(i, l, jv) + &
              alf1q(i)*szz(i, k, lp, jv) + 5.*(alf2(i)*(sz(i,k,lp,jv)-fz(i,l, &
              jv))-alf3(i)*temptm)
            sz(i, k, lp, jv) = alf(i)*fz(i, l, jv) + &
              alf1(i)*sz(i, k, lp, jv) + 3.*temptm
            ssxz(i, k, lp, jv) = alf(i)*fxz(i, l, jv) + &
              alf1(i)*ssxz(i, k, lp, jv) + 3.*(alf(i)*ssx(i,k,lp,jv)-alf1(i)* &
              fx(i,l,jv))
            syz(i, k, lp, jv) = alf(i)*fyz(i, l, jv) + &
              alf1(i)*syz(i, k, lp, jv) + 3.*(alf(i)*sy(i,k,lp,jv)-alf1(i)*fy &
              (i,l,jv))
            ssx(i, k, lp, jv) = ssx(i, k, lp, jv) + fx(i, l, jv)
            sy(i, k, lp, jv) = sy(i, k, lp, jv) + fy(i, l, jv)
            ssxx(i, k, lp, jv) = ssxx(i, k, lp, jv) + fxx(i, l, jv)
            ssxy(i, k, lp, jv) = ssxy(i, k, lp, jv) + fxy(i, l, jv)
            syy(i, k, lp, jv) = syy(i, k, lp, jv) + fyy(i, l, jv)

          END IF

        END DO
      END DO

    END DO

    ! fin de la boucle principale sur les latitudes

  END DO

  DO l = 1, llm
    DO j = 1, jjp1
      sm(iip1, j, l) = sm(1, j, l)
      s0(iip1, j, l, ntra) = s0(1, j, l, ntra)
      ssx(iip1, j, l, ntra) = ssx(1, j, l, ntra)
      sy(iip1, j, l, ntra) = sy(1, j, l, ntra)
      sz(iip1, j, l, ntra) = sz(1, j, l, ntra)
    END DO
  END DO
  ! C-------------------------------------------------------------
  ! *** Test : diag de la qqtite totale de tarceur
  ! dans l'atmosphere avant l'advection en z
  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        sqf = sqf + s0(i, j, l, ntra)
      END DO
    END DO
  END DO
  PRINT *, '-------- DIAG DANS ADVZ - SORTIE ---------'
  PRINT *, 'sqf=', sqf

  RETURN
END SUBROUTINE advzp
