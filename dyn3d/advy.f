
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advy.F,v 1.1.1.1 2004/05/19
! 12:53:06 lmdzadmin Exp $

SUBROUTINE advy(limit, dty, pbarv, sm, s0, sx, sy, sz)
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE comgeom
  IMPLICIT NONE

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! first-order moments (SOM) advection of tracer in Y direction  C
  ! C
  ! Source : Pascal Simon ( Meteo, CNRM )			 C
  ! Adaptation : A.A. (LGGE) 					 C
  ! Derniere Modif : 15/12/94 LAST
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
  PARAMETER (ntra=1)

  REAL dty
  REAL, INTENT (IN) :: pbarv(iip1, jjm, llm)

  ! moments: SM  total mass in each grid box
  ! S0  mass of tracer in each grid box
  ! Si  1rst order moment in i direction

  REAL sm(iip1, jjp1, llm), s0(iip1, jjp1, llm, ntra)
  REAL sx(iip1, jjp1, llm, ntra), sy(iip1, jjp1, llm, ntra), &
    sz(iip1, jjp1, llm, ntra)


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

  REAL f0(iim, 0:jjp1, ntra), fm(iim, 0:jjp1)
  REAL fx(iim, jjm, ntra), fy(iim, jjm, ntra)
  REAL fz(iim, jjm, ntra)
  REAL s00(ntra)
  REAL sm0 ! Just temporal variable

  ! work arrays

  REAL alf(iim, 0:jjp1), alf1(iim, 0:jjp1)
  REAL alfq(iim, 0:jjp1), alf1q(iim, 0:jjp1)
  REAL temptm ! Just temporal variable

  ! Special pour poles

  REAL sbms, sfms, sfzs, sbmn, sfmn, sfzn
  REAL sns0(ntra), snsz(ntra), snsm
  REAL s1v(llm), slatv(llm)
  REAL qy1(iim, llm, ntra), qylat(iim, llm, ntra)
  REAL cx1(llm, ntra), cxlat(llm, ntra)
  REAL cy1(llm, ntra), cylat(llm, ntra)
  REAL z1(iim), zcos(iim), zsin(iim)
  REAL smpn, smps, s0pn, s0ps
  REAL ssum
  EXTERNAL ssum

  REAL sqi, sqf
  LOGICAL limit

  lon = iim ! rem : Il est possible qu'un pbl. arrive ici
  lat = jjp1 ! a cause des dim. differentes entre les
  niv = llm


  ! the moments Fi are used as temporary storage for
  ! portions of the grid boxes in transit at the current level

  ! work arrays


  DO l = 1, llm
    DO j = 1, jjm
      DO i = 1, iip1
        vgri(i, j, llm+1-l) = -1.*pbarv(i, j, l)
      END DO
    END DO
    DO i = 1, iip1
      vgri(i, 0, l) = 0.
      vgri(i, jjp1, l) = 0.
    END DO
  END DO

  DO l = 1, niv

    ! place limits on appropriate moments before transport
    ! (if flux-limiting is to be applied)

    IF (.NOT. limit) GO TO 11

    DO jv = 1, ntra
      DO k = 1, lat
        DO i = 1, lon
          sy(i, k, l, jv) = sign(amin1(amax1(s0(i,k,l,jv), &
            0.),abs(sy(i,k,l,jv))), sy(i,k,l,jv))
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

    END DO

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,0,l)<=0.) THEN

          f0(i, 0, jv) = alf(i, 0)*(s0(i,1,l,jv)-alf1(i,0)*sy(i,1,l,jv))

          s00(jv) = s00(jv) + f0(i, 0, jv)
          s0(i, 1, l, jv) = s0(i, 1, l, jv) - f0(i, 0, jv)
          sy(i, 1, l, jv) = alf1q(i, 0)*sy(i, 1, l, jv)
          sx(i, 1, l, jv) = alf1(i, 0)*sx(i, 1, l, jv)
          sz(i, 1, l, jv) = alf1(i, 0)*sz(i, 1, l, jv)

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

    DO i = 1, lon

      IF (vgri(i,0,l)>0.) THEN
        sm(i, 1, l) = sm(i, 1, l) + fm(i, 0)
        alf(i, 0) = fm(i, 0)/sm(i, 1, l)
      END IF

      alf1(i, 0) = 1. - alf(i, 0)

    END DO

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,0,l)>0.) THEN

          temptm = alf(i, 0)*s0(i, 1, l, jv) - alf1(i, 0)*f0(i, 0, jv)
          s0(i, 1, l, jv) = s0(i, 1, l, jv) + f0(i, 0, jv)
          sy(i, 1, l, jv) = alf1(i, 0)*sy(i, 1, l, jv) + 3.*temptm

        END IF

      END DO
    END DO

    ! calculate flux and moments between adjacent boxes
    ! 1- create temporary moments/masses for partial boxes in transit
    ! 2- reajusts moments remaining in the box

    ! flux from KP to K if V(K).lt.0 and from K to KP if V(K).gt.0

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

      END DO
    END DO

    DO jv = 1, ntra
      DO k = 1, lat - 1
        kp = k + 1
        DO i = 1, lon

          IF (vgri(i,k,l)<0.) THEN

            f0(i, k, jv) = alf(i, k)*(s0(i,kp,l,jv)-alf1(i,k)*sy(i,kp,l,jv))
            fy(i, k, jv) = alfq(i, k)*sy(i, kp, l, jv)
            fx(i, k, jv) = alf(i, k)*sx(i, kp, l, jv)
            fz(i, k, jv) = alf(i, k)*sz(i, kp, l, jv)

            s0(i, kp, l, jv) = s0(i, kp, l, jv) - f0(i, k, jv)
            sy(i, kp, l, jv) = alf1q(i, k)*sy(i, kp, l, jv)
            sx(i, kp, l, jv) = sx(i, kp, l, jv) - fx(i, k, jv)
            sz(i, kp, l, jv) = sz(i, kp, l, jv) - fz(i, k, jv)

          ELSE

            f0(i, k, jv) = alf(i, k)*(s0(i,k,l,jv)+alf1(i,k)*sy(i,k,l,jv))
            fy(i, k, jv) = alfq(i, k)*sy(i, k, l, jv)
            fx(i, k, jv) = alf(i, k)*sx(i, k, l, jv)
            fz(i, k, jv) = alf(i, k)*sz(i, k, l, jv)

            s0(i, k, l, jv) = s0(i, k, l, jv) - f0(i, k, jv)
            sy(i, k, l, jv) = alf1q(i, k)*sy(i, k, l, jv)
            sx(i, k, l, jv) = sx(i, k, l, jv) - fx(i, k, jv)
            sz(i, k, l, jv) = sz(i, k, l, jv) - fz(i, k, jv)

          END IF

        END DO
      END DO
    END DO

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

        alf1(i, k) = 1. - alf(i, k)

      END DO
    END DO

    DO jv = 1, ntra
      DO k = 1, lat - 1
        kp = k + 1
        DO i = 1, lon

          IF (vgri(i,k,l)<0.) THEN

            temptm = -alf(i, k)*s0(i, k, l, jv) + alf1(i, k)*f0(i, k, jv)
            s0(i, k, l, jv) = s0(i, k, l, jv) + f0(i, k, jv)
            sy(i, k, l, jv) = alf(i, k)*fy(i, k, jv) + &
              alf1(i, k)*sy(i, k, l, jv) + 3.*temptm
            sx(i, k, l, jv) = sx(i, k, l, jv) + fx(i, k, jv)
            sz(i, k, l, jv) = sz(i, k, l, jv) + fz(i, k, jv)

          ELSE

            temptm = alf(i, k)*s0(i, kp, l, jv) - alf1(i, k)*f0(i, k, jv)
            s0(i, kp, l, jv) = s0(i, kp, l, jv) + f0(i, k, jv)
            sy(i, kp, l, jv) = alf(i, k)*fy(i, k, jv) + &
              alf1(i, k)*sy(i, kp, l, jv) + 3.*temptm
            sx(i, kp, l, jv) = sx(i, kp, l, jv) + fx(i, k, jv)
            sz(i, kp, l, jv) = sz(i, kp, l, jv) + fz(i, k, jv)

          END IF

        END DO
      END DO
    END DO

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

    END DO

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,k,l)>=0.) THEN
          f0(i, k, jv) = alf(i, k)*(s0(i,k,l,jv)+alf1(i,k)*sy(i,k,l,jv))
          s00(jv) = s00(jv) + f0(i, k, jv)

          s0(i, k, l, jv) = s0(i, k, l, jv) - f0(i, k, jv)
          sy(i, k, l, jv) = alf1q(i, k)*sy(i, k, l, jv)
          sx(i, k, l, jv) = alf1(i, k)*sx(i, k, l, jv)
          sz(i, k, l, jv) = alf1(i, k)*sz(i, k, l, jv)
        END IF

      END DO
    END DO

    DO i = 1, lon
      IF (vgri(i,k,l)<0.) THEN
        fm(i, k) = -vgri(i, k, l)*dty
        alf(i, k) = fm(i, k)/sm0
      END IF
    END DO

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

      alf1(i, k) = 1. - alf(i, k)

    END DO

    DO jv = 1, ntra
      DO i = 1, lon

        IF (vgri(i,k,l)<0.) THEN

          temptm = -alf(i, k)*s0(i, k, l, jv) + alf1(i, k)*f0(i, k, jv)
          s0(i, k, l, jv) = s0(i, k, l, jv) + f0(i, k, jv)
          sy(i, k, l, jv) = alf1(i, k)*sy(i, k, l, jv) + 3.*temptm

        END IF

      END DO
    END DO

  END DO

  RETURN
END SUBROUTINE advy

