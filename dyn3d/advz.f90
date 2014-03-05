
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advz.F,v 1.2 2005/05/25 13:10:09
! fairhead Exp $

SUBROUTINE advz(limit, dtz, w, sm, s0, sx, sy, sz)
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  IMPLICIT NONE

  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! C
  ! first-order moments (FOM) advection of tracer in Z direction  C
  ! C
  ! Source : Pascal Simon (Meteo,CNRM)                            C
  ! Adaptation : A.Armengaud (LGGE) juin 94                       C
  ! C
  ! C
  ! sont des arguments d'entree pour le s-pg...                   C
  ! C
  ! dq est l'argument de sortie pour le s-pg                      C
  ! C
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! parametres principaux du modele


  ! Arguments :
  ! -----------
  ! dtz : frequence fictive d'appel du transport
  ! w : flux de masse en z en Pa.m2.s-1

  INTEGER ntra
  PARAMETER (ntra=1)

  REAL, INTENT (IN) :: dtz
  REAL w(iip1, jjp1, llm)

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

  REAL wgri(iip1, jjp1, 0:llm)


  ! the moments F are used as temporary  storage for
  ! portions of grid boxes in transit at the current latitude

  REAL fm(iim, llm)
  REAL f0(iim, llm, ntra), fx(iim, llm, ntra)
  REAL fy(iim, llm, ntra), fz(iim, llm, ntra)

  ! work arrays

  REAL alf(iim), alf1(iim), alfq(iim), alf1q(iim)
  REAL temptm ! Just temporal variable
  REAL sqi, sqf

  LOGICAL limit
  INTEGER lon, lat, niv
  INTEGER i, j, jv, k, l, lp

  lon = iim
  lat = jjp1
  niv = llm

  ! *** Test : diag de la qqtite totale de traceur
  ! dans l'atmosphere avant l'advection en z
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
  PRINT *, '-------- DIAG DANS ADVZ - ENTREE ---------'
  PRINT *, 'sqi=', sqi

  ! -----------------------------------------------------------------
  ! Interface : adaptation nouveau modele
  ! -------------------------------------

  ! Conversion du flux de masse en kg.s-1

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        ! wgri (i,j,llm+1-l) =  w (i,j,l) / g
        wgri(i, j, llm+1-l) = w(i, j, l)
        ! wgri (i,j,0) = 0.                ! a detruire ult.
        ! wgri (i,j,l) = 0.1               !    w (i,j,l)
        ! wgri (i,j,llm) = 0.              ! a detruire ult.
      END DO
    END DO
  END DO
  DO j = 1, jjp1
    DO i = 1, iip1
      wgri(i, j, 0) = 0.
    END DO
  END DO

  ! -----------------------------------------------------------------

  ! start here
  ! boucle sur les latitudes

  DO k = 1, lat

    ! place limits on appropriate moments before transport
    ! (if flux-limiting is to be applied)

    IF (.NOT. limit) GO TO 101

    DO jv = 1, ntra
      DO l = 1, niv
        DO i = 1, lon
          sz(i, k, l, jv) = sign(amin1(amax1(s0(i,k,l,jv), &
            0.),abs(sz(i,k,l,jv))), sz(i,k,l,jv))
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

      END DO

      DO jv = 1, ntra
        DO i = 1, lon

          IF (wgri(i,k,l)<0.) THEN

            f0(i, l, jv) = alf(i)*(s0(i,k,lp,jv)-alf1(i)*sz(i,k,lp,jv))
            fz(i, l, jv) = alfq(i)*sz(i, k, lp, jv)
            fx(i, l, jv) = alf(i)*sx(i, k, lp, jv)
            fy(i, l, jv) = alf(i)*sy(i, k, lp, jv)

            s0(i, k, lp, jv) = s0(i, k, lp, jv) - f0(i, l, jv)
            sz(i, k, lp, jv) = alf1q(i)*sz(i, k, lp, jv)
            sx(i, k, lp, jv) = sx(i, k, lp, jv) - fx(i, l, jv)
            sy(i, k, lp, jv) = sy(i, k, lp, jv) - fy(i, l, jv)

          ELSE

            f0(i, l, jv) = alf(i)*(s0(i,k,l,jv)+alf1(i)*sz(i,k,l,jv))
            fz(i, l, jv) = alfq(i)*sz(i, k, l, jv)
            fx(i, l, jv) = alf(i)*sx(i, k, l, jv)
            fy(i, l, jv) = alf(i)*sy(i, k, l, jv)

            s0(i, k, l, jv) = s0(i, k, l, jv) - f0(i, l, jv)
            sz(i, k, l, jv) = alf1q(i)*sz(i, k, l, jv)
            sx(i, k, l, jv) = sx(i, k, l, jv) - fx(i, l, jv)
            sy(i, k, l, jv) = sy(i, k, l, jv) - fy(i, l, jv)

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

      END DO

      DO jv = 1, ntra
        DO i = 1, lon

          IF (wgri(i,k,l)<0.) THEN

            temptm = -alf(i)*s0(i, k, l, jv) + alf1(i)*f0(i, l, jv)
            s0(i, k, l, jv) = s0(i, k, l, jv) + f0(i, l, jv)
            sz(i, k, l, jv) = alf(i)*fz(i, l, jv) + alf1(i)*sz(i, k, l, jv) + &
              3.*temptm
            sx(i, k, l, jv) = sx(i, k, l, jv) + fx(i, l, jv)
            sy(i, k, l, jv) = sy(i, k, l, jv) + fy(i, l, jv)

          ELSE

            temptm = alf(i)*s0(i, k, lp, jv) - alf1(i)*f0(i, l, jv)
            s0(i, k, lp, jv) = s0(i, k, lp, jv) + f0(i, l, jv)
            sz(i, k, lp, jv) = alf(i)*fz(i, l, jv) + &
              alf1(i)*sz(i, k, lp, jv) + 3.*temptm
            sx(i, k, lp, jv) = sx(i, k, lp, jv) + fx(i, l, jv)
            sy(i, k, lp, jv) = sy(i, k, lp, jv) + fy(i, l, jv)

          END IF

        END DO
      END DO

    END DO

    ! fin de la boucle principale sur les latitudes

  END DO

  ! *** ------------------- bouclage cyclique  en X ------------

  ! DO l = 1,llm
  ! DO j = 1,jjp1
  ! SM(iip1,j,l) = SM(1,j,l)
  ! S0(iip1,j,l,ntra) = S0(1,j,l,ntra)
  ! sx(iip1,j,l,ntra) = sx(1,j,l,ntra)
  ! sy(iip1,j,l,ntra) = sy(1,j,l,ntra)
  ! sz(iip1,j,l,ntra) = sz(1,j,l,ntra)
  ! ENDDO
  ! ENDDO

  ! -------------------------------------------------------------
  ! *** Test : diag de la qqtite totale de traceur
  ! dans l'atmosphere avant l'advection en z
  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iim
        ! IM 240305            sqf = sqf + S0(i,j,l,9)
        sqf = sqf + s0(i, j, l, ntra)
      END DO
    END DO
  END DO
  PRINT *, '-------- DIAG DANS ADVZ - SORTIE ---------'
  PRINT *, 'sqf=', sqf

  ! -------------------------------------------------------------
  RETURN
END SUBROUTINE advz
! _______________________________________________________________
! _______________________________________________________________
