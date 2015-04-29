
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advyp.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE advyp(limit, dty, pbarv, sm, s0, ssx, sy, sz, ssxx, ssxy, ssxz, &
    syy, syz, szz, ntra)
  USE dimens_m
  USE comconst
  USE paramet_m
  USE disvert_m
  USE comgeom
  IMPLICIT NONE
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! second-order moments (SOM) advection of tracer in Y direction  C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! Source : Pascal Simon ( Meteo, CNRM )			 C
  ! Adaptation : A.A. (LGGE) 					 C
  ! Derniere Modif : 19/10/95 LAST
  ! C
  ! sont les arguments d'entree pour le s-pg			 C
  ! C
  ! argument de sortie du s-pg					 C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Rem : Probleme aux poles il faut reecrire ce cas specifique
  ! Attention au sens de l'indexation

  ! parametres principaux du modele



  ! Arguments :
  ! ----------
  ! dty : frequence fictive d'appel du transport
  ! parbu,pbarv : flux de masse en x et y en Pa.m2.s-1

  INTEGER lon, lat, niv
  INTEGER i, j, jv, k, kp, l
  INTEGER ntra
  ! PARAMETER (ntra = 1)

  REAL dty
  REAL, INTENT (IN) :: pbarv(iip1, jjm, llm)

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

  REAL vgri(iip1, 0:jjp1, llm)

  ! Rem : UGRI et WGRI ne sont pas utilises dans
  ! cette subroutine ( advection en y uniquement )
  ! Rem 2 :le dimensionnement de VGRI depend de celui de pbarv

  ! the moments F are similarly defined and used as temporary
  ! storage for portions of the grid boxes in transit

  ! the moments Fij are used as temporary storage for
  ! portions of the grid boxes in transit at the current level

  ! work arrays


  REAL f0(iim, 0:jjp1, ntra), fm(iim, 0:jjp1)
  REAL fx(iim, jjm, ntra), fy(iim, jjm, ntra)
  REAL fz(iim, jjm, ntra)
  REAL fxx(iim, jjm, ntra), fxy(iim, jjm, ntra)
  REAL fxz(iim, jjm, ntra), fyy(iim, jjm, ntra)
  REAL fyz(iim, jjm, ntra), fzz(iim, jjm, ntra)
  REAL s00(ntra)
  REAL sm0 ! Just temporal variable

  ! work arrays

  REAL alf(iim, 0:jjp1), alf1(iim, 0:jjp1)
  REAL alfq(iim, 0:jjp1), alf1q(iim, 0:jjp1)
  REAL alf2(iim, 0:jjp1), alf3(iim, 0:jjp1)
  REAL alf4(iim, 0:jjp1)
  REAL temptm ! Just temporal variable
  REAL slpmax, s1max, s1new, s2new

  ! Special pour poles

  REAL sbms, sfms, sfzs, sbmn, sfmn, sfzn
  REAL ssum
  EXTERNAL ssum

  REAL sqi, sqf
  LOGICAL limit

  lon = iim ! rem : Il est possible qu'un pbl. arrive ici
  lat = jjp1 ! a cause des dim. differentes entre les
  niv = llm !       tab. S et VGRI

  ! -----------------------------------------------------------------
  ! initialisations

  sbms = 0.
  sfms = 0.
  sfzs = 0.
  sbmn = 0.
  sfmn = 0.
  sfzn = 0.

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
  PRINT *, '---------- DIAG DANS ADVY - ENTREE --------'
  PRINT *, 'sqi=', sqi

  ! -----------------------------------------------------------------
  ! Interface : adaptation nouveau modele
  ! -------------------------------------

  ! Conversion des flux de masses en kg
  ! -AA 20/10/94  le signe -1 est necessaire car indexation opposee

  DO l = 1, llm
    DO j = 1, jjm
      DO i = 1, iip1
        vgri(i, j, llm+1-l) = -1.*pbarv(i, j, l)
      END DO
    END DO
  END DO

  ! AA Initialisation de flux fictifs aux bords sup. des boites pol.

  DO l = 1, llm
    DO i = 1, iip1
      vgri(i, 0, l) = 0.
      vgri(i, jjp1, l) = 0.
    END DO
  END DO

  ! ----------------- START HERE -----------------------
  ! boucle sur les niveaux

  DO l = 1, niv

    ! place limits on appropriate moments before transport
    ! (if flux-limiting is to be applied)

    IF (.NOT. limit) GO TO 11

    DO jv = 1, ntra
      DO k = 1, lat
        DO i = 1, lon
          IF (s0(i,k,l,jv)>0.) THEN
            slpmax = amax1(s0(i,k,l,jv), 0.)
            s1max = 1.5*slpmax
            s1new = amin1(s1max, amax1(-s1max,sy(i,k,l,jv)))
            s2new = amin1(2.*slpmax-abs(s1new)/3., amax1(abs( &
              s1new)-slpmax,syy(i,k,l,jv)))
            sy(i, k, l, jv) = s1new
            syy(i, k, l, jv) = s2new
            ssxy(i, k, l, jv) = amin1(slpmax, amax1(-slpmax,ssxy(i,k,l,jv)))
            syz(i, k, l, jv) = amin1(slpmax, amax1(-slpmax,syz(i,k,l,jv)))
          ELSE
            sy(i, k, l, jv) = 0.
            syy(i, k, l, jv) = 0.
            ssxy(i, k, l, jv) = 0.
            syz(i, k, l, jv) = 0.
          END IF
        END DO
      END DO
    END DO

11  CONTINUE

    ! le flux a travers le pole Nord est traite separement

    sm0 = 0.
    DO jv = 1, ntra
      s00(jv) = 0.
    END DO

    DO i = 1, lon

      IF (vgri(i,0,l)<=0.) THEN
        fm(i, 0) = -vgri(i, 0, l)*dty
        alf(i, 0) = fm(i, 0)/sm(i, 1, l)
        sm(i, 1, l) = sm(i, 1, l) - fm(i, 0)
        sm0 = sm0 + fm(i, 0)
      END IF

      alfq(i, 0) = alf(i, 0)*alf(i, 0)
      alf1(i, 0) = 1. - alf(i, 0)
      alf1q(i, 0) = alf1(i, 0)*alf1(i, 0)
      alf2(i, 0) = alf1(i, 0) - alf(i, 0)
      alf3(i, 0) = alf(i, 0)*alfq(i, 0)
      alf4(i, 0) = alf1(i, 0)*alf1q(i, 0)

    END DO
    ! print*,'ADVYP 21'

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,0,l)<=0.) THEN

          f0(i, 0, jv) = alf(i, 0)*(s0(i,1,l,jv)-alf1(i,0)*(sy(i,1,l, &
            jv)-alf2(i,0)*syy(i,1,l,jv)))

          s00(jv) = s00(jv) + f0(i, 0, jv)
          s0(i, 1, l, jv) = s0(i, 1, l, jv) - f0(i, 0, jv)
          sy(i, 1, l, jv) = alf1q(i, 0)*(sy(i,1,l,jv)+3.*alf(i,0)*syy(i,1,l, &
            jv))
          syy(i, 1, l, jv) = alf4(i, 0)*syy(i, 1, l, jv)
          ssx(i, 1, l, jv) = alf1(i, 0)*(ssx(i,1,l,jv)+alf(i,0)*ssxy(i,1,l,jv &
            ))
          sz(i, 1, l, jv) = alf1(i, 0)*(sz(i,1,l,jv)+alf(i,0)*ssxz(i,1,l,jv))
          ssxx(i, 1, l, jv) = alf1(i, 0)*ssxx(i, 1, l, jv)
          ssxz(i, 1, l, jv) = alf1(i, 0)*ssxz(i, 1, l, jv)
          szz(i, 1, l, jv) = alf1(i, 0)*szz(i, 1, l, jv)
          ssxy(i, 1, l, jv) = alf1q(i, 0)*ssxy(i, 1, l, jv)
          syz(i, 1, l, jv) = alf1q(i, 0)*syz(i, 1, l, jv)

        END IF

      END DO
    END DO

    DO i = 1, lon
      IF (vgri(i,0,l)>0.) THEN
        fm(i, 0) = vgri(i, 0, l)*dty
        alf(i, 0) = fm(i, 0)/sm0
      END IF
    END DO

    DO jv = 1, ntra
      DO i = 1, lon
        IF (vgri(i,0,l)>0.) THEN
          f0(i, 0, jv) = alf(i, 0)*s00(jv)
        END IF
      END DO
    END DO

    ! puts the temporary moments Fi into appropriate neighboring boxes

    ! print*,'av ADVYP 25'
    DO i = 1, lon

      IF (vgri(i,0,l)>0.) THEN
        sm(i, 1, l) = sm(i, 1, l) + fm(i, 0)
        alf(i, 0) = fm(i, 0)/sm(i, 1, l)
      END IF

      alfq(i, 0) = alf(i, 0)*alf(i, 0)
      alf1(i, 0) = 1. - alf(i, 0)
      alf1q(i, 0) = alf1(i, 0)*alf1(i, 0)
      alf2(i, 0) = alf1(i, 0) - alf(i, 0)
      alf3(i, 0) = alf1(i, 0)*alf(i, 0)

    END DO
    ! print*,'av ADVYP 25'

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,0,l)>0.) THEN

          temptm = alf(i, 0)*s0(i, 1, l, jv) - alf1(i, 0)*f0(i, 0, jv)
          s0(i, 1, l, jv) = s0(i, 1, l, jv) + f0(i, 0, jv)
          syy(i, 1, l, jv) = alf1q(i, 0)*syy(i, 1, l, jv) + &
            5.*(alf3(i,0)*sy(i,1,l,jv)-alf2(i,0)*temptm)
          sy(i, 1, l, jv) = alf1(i, 0)*sy(i, 1, l, jv) + 3.*temptm
          ssxy(i, 1, l, jv) = alf1(i, 0)*ssxy(i, 1, l, jv) + &
            3.*alf(i, 0)*ssx(i, 1, l, jv)
          syz(i, 1, l, jv) = alf1(i, 0)*syz(i, 1, l, jv) + &
            3.*alf(i, 0)*sz(i, 1, l, jv)

        END IF

      END DO
    END DO

    ! calculate flux and moments between adjacent boxes
    ! 1- create temporary moments/masses for partial boxes in transit
    ! 2- reajusts moments remaining in the box

    ! flux from KP to K if V(K).lt.0 and from K to KP if V(K).gt.0

    ! print*,'av ADVYP 30'
    DO k = 1, lat - 1
      kp = k + 1
      DO i = 1, lon

        IF (vgri(i,k,l)<0.) THEN
          fm(i, k) = -vgri(i, k, l)*dty
          alf(i, k) = fm(i, k)/sm(i, kp, l)
          sm(i, kp, l) = sm(i, kp, l) - fm(i, k)
        ELSE
          fm(i, k) = vgri(i, k, l)*dty
          alf(i, k) = fm(i, k)/sm(i, k, l)
          sm(i, k, l) = sm(i, k, l) - fm(i, k)
        END IF

        alfq(i, k) = alf(i, k)*alf(i, k)
        alf1(i, k) = 1. - alf(i, k)
        alf1q(i, k) = alf1(i, k)*alf1(i, k)
        alf2(i, k) = alf1(i, k) - alf(i, k)
        alf3(i, k) = alf(i, k)*alfq(i, k)
        alf4(i, k) = alf1(i, k)*alf1q(i, k)

      END DO
    END DO
    ! print*,'ap ADVYP 30'

    DO jv = 1, ntra
      DO k = 1, lat - 1
        kp = k + 1
        DO i = 1, lon

          IF (vgri(i,k,l)<0.) THEN

            f0(i, k, jv) = alf(i, k)*(s0(i,kp,l,jv)-alf1(i,k)*(sy(i,kp,l, &
              jv)-alf2(i,k)*syy(i,kp,l,jv)))
            fy(i, k, jv) = alfq(i, k)*(sy(i,kp,l,jv)-3.*alf1(i,k)*syy(i,kp,l, &
              jv))
            fyy(i, k, jv) = alf3(i, k)*syy(i, kp, l, jv)
            fx(i, k, jv) = alf(i, k)*(ssx(i,kp,l,jv)-alf1(i,k)*ssxy(i,kp,l,jv &
              ))
            fz(i, k, jv) = alf(i, k)*(sz(i,kp,l,jv)-alf1(i,k)*syz(i,kp,l,jv))
            fxy(i, k, jv) = alfq(i, k)*ssxy(i, kp, l, jv)
            fyz(i, k, jv) = alfq(i, k)*syz(i, kp, l, jv)
            fxx(i, k, jv) = alf(i, k)*ssxx(i, kp, l, jv)
            fxz(i, k, jv) = alf(i, k)*ssxz(i, kp, l, jv)
            fzz(i, k, jv) = alf(i, k)*szz(i, kp, l, jv)

            s0(i, kp, l, jv) = s0(i, kp, l, jv) - f0(i, k, jv)
            sy(i, kp, l, jv) = alf1q(i, k)*(sy(i,kp,l,jv)+3.*alf(i,k)*syy(i, &
              kp,l,jv))
            syy(i, kp, l, jv) = alf4(i, k)*syy(i, kp, l, jv)
            ssx(i, kp, l, jv) = ssx(i, kp, l, jv) - fx(i, k, jv)
            sz(i, kp, l, jv) = sz(i, kp, l, jv) - fz(i, k, jv)
            ssxx(i, kp, l, jv) = ssxx(i, kp, l, jv) - fxx(i, k, jv)
            ssxz(i, kp, l, jv) = ssxz(i, kp, l, jv) - fxz(i, k, jv)
            szz(i, kp, l, jv) = szz(i, kp, l, jv) - fzz(i, k, jv)
            ssxy(i, kp, l, jv) = alf1q(i, k)*ssxy(i, kp, l, jv)
            syz(i, kp, l, jv) = alf1q(i, k)*syz(i, kp, l, jv)

          ELSE

            f0(i, k, jv) = alf(i, k)*(s0(i,k,l,jv)+alf1(i,k)*(sy(i,k,l, &
              jv)+alf2(i,k)*syy(i,k,l,jv)))
            fy(i, k, jv) = alfq(i, k)*(sy(i,k,l,jv)+3.*alf1(i,k)*syy(i,k,l,jv &
              ))
            fyy(i, k, jv) = alf3(i, k)*syy(i, k, l, jv)
            fx(i, k, jv) = alf(i, k)*(ssx(i,k,l,jv)+alf1(i,k)*ssxy(i,k,l,jv))
            fz(i, k, jv) = alf(i, k)*(sz(i,k,l,jv)+alf1(i,k)*syz(i,k,l,jv))
            fxy(i, k, jv) = alfq(i, k)*ssxy(i, k, l, jv)
            fyz(i, k, jv) = alfq(i, k)*syz(i, k, l, jv)
            fxx(i, k, jv) = alf(i, k)*ssxx(i, k, l, jv)
            fxz(i, k, jv) = alf(i, k)*ssxz(i, k, l, jv)
            fzz(i, k, jv) = alf(i, k)*szz(i, k, l, jv)

            s0(i, k, l, jv) = s0(i, k, l, jv) - f0(i, k, jv)
            sy(i, k, l, jv) = alf1q(i, k)*(sy(i,k,l,jv)-3.*alf(i,k)*syy(i,k,l &
              ,jv))
            syy(i, k, l, jv) = alf4(i, k)*syy(i, k, l, jv)
            ssx(i, k, l, jv) = ssx(i, k, l, jv) - fx(i, k, jv)
            sz(i, k, l, jv) = sz(i, k, l, jv) - fz(i, k, jv)
            ssxx(i, k, l, jv) = ssxx(i, k, l, jv) - fxx(i, k, jv)
            ssxz(i, k, l, jv) = ssxz(i, k, l, jv) - fxz(i, k, jv)
            szz(i, k, l, jv) = szz(i, k, l, jv) - fzz(i, k, jv)
            ssxy(i, k, l, jv) = alf1q(i, k)*ssxy(i, k, l, jv)
            syz(i, k, l, jv) = alf1q(i, k)*syz(i, k, l, jv)

          END IF

        END DO
      END DO
    END DO
    ! print*,'ap ADVYP 31'

    ! puts the temporary moments Fi into appropriate neighboring boxes

    DO k = 1, lat - 1
      kp = k + 1
      DO i = 1, lon

        IF (vgri(i,k,l)<0.) THEN
          sm(i, k, l) = sm(i, k, l) + fm(i, k)
          alf(i, k) = fm(i, k)/sm(i, k, l)
        ELSE
          sm(i, kp, l) = sm(i, kp, l) + fm(i, k)
          alf(i, k) = fm(i, k)/sm(i, kp, l)
        END IF

        alfq(i, k) = alf(i, k)*alf(i, k)
        alf1(i, k) = 1. - alf(i, k)
        alf1q(i, k) = alf1(i, k)*alf1(i, k)
        alf2(i, k) = alf1(i, k) - alf(i, k)
        alf3(i, k) = alf1(i, k)*alf(i, k)

      END DO
    END DO
    ! print*,'ap ADVYP 32'

    DO jv = 1, ntra
      DO k = 1, lat - 1
        kp = k + 1
        DO i = 1, lon

          IF (vgri(i,k,l)<0.) THEN

            temptm = -alf(i, k)*s0(i, k, l, jv) + alf1(i, k)*f0(i, k, jv)
            s0(i, k, l, jv) = s0(i, k, l, jv) + f0(i, k, jv)
            syy(i, k, l, jv) = alfq(i, k)*fyy(i, k, jv) + &
              alf1q(i, k)*syy(i, k, l, jv) + 5.*(alf3(i,k)*(fy(i,k,jv)-sy(i, &
              k,l,jv))+alf2(i,k)*temptm)
            sy(i, k, l, jv) = alf(i, k)*fy(i, k, jv) + &
              alf1(i, k)*sy(i, k, l, jv) + 3.*temptm
            ssxy(i, k, l, jv) = alf(i, k)*fxy(i, k, jv) + &
              alf1(i, k)*ssxy(i, k, l, jv) + 3.*(alf1(i,k)*fx(i,k,jv)-alf(i,k &
              )*ssx(i,k,l,jv))
            syz(i, k, l, jv) = alf(i, k)*fyz(i, k, jv) + &
              alf1(i, k)*syz(i, k, l, jv) + 3.*(alf1(i,k)*fz(i,k,jv)-alf(i,k) &
              *sz(i,k,l,jv))
            ssx(i, k, l, jv) = ssx(i, k, l, jv) + fx(i, k, jv)
            sz(i, k, l, jv) = sz(i, k, l, jv) + fz(i, k, jv)
            ssxx(i, k, l, jv) = ssxx(i, k, l, jv) + fxx(i, k, jv)
            ssxz(i, k, l, jv) = ssxz(i, k, l, jv) + fxz(i, k, jv)
            szz(i, k, l, jv) = szz(i, k, l, jv) + fzz(i, k, jv)

          ELSE

            temptm = alf(i, k)*s0(i, kp, l, jv) - alf1(i, k)*f0(i, k, jv)
            s0(i, kp, l, jv) = s0(i, kp, l, jv) + f0(i, k, jv)
            syy(i, kp, l, jv) = alfq(i, k)*fyy(i, k, jv) + &
              alf1q(i, k)*syy(i, kp, l, jv) + 5.*(alf3(i,k)*(sy(i,kp,l, &
              jv)-fy(i,k,jv))-alf2(i,k)*temptm)
            sy(i, kp, l, jv) = alf(i, k)*fy(i, k, jv) + &
              alf1(i, k)*sy(i, kp, l, jv) + 3.*temptm
            ssxy(i, kp, l, jv) = alf(i, k)*fxy(i, k, jv) + &
              alf1(i, k)*ssxy(i, kp, l, jv) + 3.*(alf(i,k)*ssx(i,kp,l,jv)- &
              alf1(i,k)*fx(i,k,jv))
            syz(i, kp, l, jv) = alf(i, k)*fyz(i, k, jv) + &
              alf1(i, k)*syz(i, kp, l, jv) + 3.*(alf(i,k)*sz(i,kp,l,jv)-alf1( &
              i,k)*fz(i,k,jv))
            ssx(i, kp, l, jv) = ssx(i, kp, l, jv) + fx(i, k, jv)
            sz(i, kp, l, jv) = sz(i, kp, l, jv) + fz(i, k, jv)
            ssxx(i, kp, l, jv) = ssxx(i, kp, l, jv) + fxx(i, k, jv)
            ssxz(i, kp, l, jv) = ssxz(i, kp, l, jv) + fxz(i, k, jv)
            szz(i, kp, l, jv) = szz(i, kp, l, jv) + fzz(i, k, jv)

          END IF

        END DO
      END DO
    END DO
    ! print*,'ap ADVYP 33'

    ! traitement special pour le pole Sud (idem pole Nord)

    k = lat

    sm0 = 0.
    DO jv = 1, ntra
      s00(jv) = 0.
    END DO

    DO i = 1, lon

      IF (vgri(i,k,l)>=0.) THEN
        fm(i, k) = vgri(i, k, l)*dty
        alf(i, k) = fm(i, k)/sm(i, k, l)
        sm(i, k, l) = sm(i, k, l) - fm(i, k)
        sm0 = sm0 + fm(i, k)
      END IF

      alfq(i, k) = alf(i, k)*alf(i, k)
      alf1(i, k) = 1. - alf(i, k)
      alf1q(i, k) = alf1(i, k)*alf1(i, k)
      alf2(i, k) = alf1(i, k) - alf(i, k)
      alf3(i, k) = alf(i, k)*alfq(i, k)
      alf4(i, k) = alf1(i, k)*alf1q(i, k)

    END DO
    ! print*,'ap ADVYP 41'

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,k,l)>=0.) THEN
          f0(i, k, jv) = alf(i, k)*(s0(i,k,l,jv)+alf1(i,k)*(sy(i,k,l, &
            jv)+alf2(i,k)*syy(i,k,l,jv)))
          s00(jv) = s00(jv) + f0(i, k, jv)

          s0(i, k, l, jv) = s0(i, k, l, jv) - f0(i, k, jv)
          sy(i, k, l, jv) = alf1q(i, k)*(sy(i,k,l,jv)-3.*alf(i,k)*syy(i,k,l, &
            jv))
          syy(i, k, l, jv) = alf4(i, k)*syy(i, k, l, jv)
          ssx(i, k, l, jv) = alf1(i, k)*(ssx(i,k,l,jv)-alf(i,k)*ssxy(i,k,l,jv &
            ))
          sz(i, k, l, jv) = alf1(i, k)*(sz(i,k,l,jv)-alf(i,k)*syz(i,k,l,jv))
          ssxx(i, k, l, jv) = alf1(i, k)*ssxx(i, k, l, jv)
          ssxz(i, k, l, jv) = alf1(i, k)*ssxz(i, k, l, jv)
          szz(i, k, l, jv) = alf1(i, k)*szz(i, k, l, jv)
          ssxy(i, k, l, jv) = alf1q(i, k)*ssxy(i, k, l, jv)
          syz(i, k, l, jv) = alf1q(i, k)*syz(i, k, l, jv)
        END IF

      END DO
    END DO
    ! print*,'ap ADVYP 42'

    DO i = 1, lon
      IF (vgri(i,k,l)<0.) THEN
        fm(i, k) = -vgri(i, k, l)*dty
        alf(i, k) = fm(i, k)/sm0
      END IF
    END DO
    ! print*,'ap ADVYP 43'

    DO jv = 1, ntra
      DO i = 1, lon
        IF (vgri(i,k,l)<0.) THEN
          f0(i, k, jv) = alf(i, k)*s00(jv)
        END IF
      END DO
    END DO

    ! puts the temporary moments Fi into appropriate neighboring boxes

    DO i = 1, lon

      IF (vgri(i,k,l)<0.) THEN
        sm(i, k, l) = sm(i, k, l) + fm(i, k)
        alf(i, k) = fm(i, k)/sm(i, k, l)
      END IF

      alfq(i, k) = alf(i, k)*alf(i, k)
      alf1(i, k) = 1. - alf(i, k)
      alf1q(i, k) = alf1(i, k)*alf1(i, k)
      alf2(i, k) = alf1(i, k) - alf(i, k)
      alf3(i, k) = alf1(i, k)*alf(i, k)

    END DO
    ! print*,'ap ADVYP 45'

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,k,l)<0.) THEN

          temptm = -alf(i, k)*s0(i, k, l, jv) + alf1(i, k)*f0(i, k, jv)
          s0(i, k, l, jv) = s0(i, k, l, jv) + f0(i, k, jv)
          syy(i, k, l, jv) = alf1q(i, k)*syy(i, k, l, jv) + &
            5.*(-alf3(i,k)*sy(i,k,l,jv)+alf2(i,k)*temptm)
          sy(i, k, l, jv) = alf1(i, k)*sy(i, k, l, jv) + 3.*temptm
          ssxy(i, k, l, jv) = alf1(i, k)*ssxy(i, k, l, jv) - &
            3.*alf(i, k)*ssx(i, k, l, jv)
          syz(i, k, l, jv) = alf1(i, k)*syz(i, k, l, jv) - &
            3.*alf(i, k)*sz(i, k, l, jv)

        END IF

      END DO
    END DO
    ! print*,'ap ADVYP 46'

  END DO

  ! --------------------------------------------------
  ! bouclage cyclique horizontal .

  DO l = 1, llm
    DO jv = 1, ntra
      DO j = 1, jjp1
        sm(iip1, j, l) = sm(1, j, l)
        s0(iip1, j, l, jv) = s0(1, j, l, jv)
        ssx(iip1, j, l, jv) = ssx(1, j, l, jv)
        sy(iip1, j, l, jv) = sy(1, j, l, jv)
        sz(iip1, j, l, jv) = sz(1, j, l, jv)
      END DO
    END DO
  END DO

  ! -------------------------------------------------------------------
  ! *** Test  negativite:

  ! DO jv = 1,ntra
  ! DO l = 1,llm
  ! DO j = 1,jjp1
  ! DO i = 1,iip1
  ! IF (s0( i,j,l,jv ).lt.0.) THEN
  ! PRINT*, '------ S0 < 0 en FIN ADVYP ---'
  ! PRINT*, 'S0(',i,j,l,jv,')=', S0(i,j,l,jv)
  ! c                 STOP
  ! ENDIF
  ! ENDDO
  ! ENDDO
  ! ENDDO
  ! ENDDO


  ! -------------------------------------------------------------------
  ! *** Test : diag de la qtite totale de traceur dans
  ! l'atmosphere avant l'advection en Y

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        sqf = sqf + s0(i, j, l, ntra)
      END DO
    END DO
  END DO
  PRINT *, '---------- DIAG DANS ADVY - SORTIE --------'
  PRINT *, 'sqf=', sqf
  ! print*,'ap ADVYP fin'

  ! -----------------------------------------------------------------

  RETURN
END SUBROUTINE advyp












