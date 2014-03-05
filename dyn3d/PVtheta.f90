SUBROUTINE pvtheta(ilon, ilev, pucov, pvcov, pteta, ztfi, zplay, zplev, &
    nbteta, theta, pvteta)
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  USE comgeom
  USE tourabs_m, ONLY: tourabs
  IMPLICIT NONE

  ! =======================================================================

  ! Auteur:  I. Musat
  ! -------

  ! Objet:
  ! ------

  ! *******************************************************************
  ! Calcul de la vorticite potentielle PVteta sur des iso-theta selon
  ! la methodologie du NCEP/NCAR :
  ! 1) on calcule la stabilite statique N**2=g/T*(dT/dz+g/cp) sur les
  ! niveaux du modele => N2
  ! 2) on interpole les vents, la temperature et le N**2 sur des isentropes
  ! (en fait sur des iso-theta) lineairement en log(theta) =>
  ! ucovteta, vcovteta, N2teta
  ! 3) on calcule la vorticite absolue sur des iso-theta => vorateta
  ! 4) on calcule la densite rho sur des iso-theta => rhoteta

  ! rhoteta = (T/theta)**(cp/R)*p0/(R*T)

  ! 5) on calcule la vorticite potentielle sur des iso-theta => PVteta

  ! PVteta = (vorateta * N2 * theta)/(g * rhoteta) ! en PVU

  ! NB: 1PVU=10**(-6) K*m**2/(s * kg)

  ! PVteta =  vorateta * N2/(g**2 * rhoteta) ! en 1/(Pa*s)


  ! *******************************************************************


  ! Variables d'entree :
  ! ilon,ilev,pucov,pvcov,pteta,ztfi,zplay,zplev,nbteta,theta
  ! -> sur la grille dynamique
  ! Variable de sortie : PVteta
  ! -> sur la grille physique
  ! =======================================================================


  ! variables Input

  INTEGER ilon
  INTEGER, INTENT (IN) :: ilev
  REAL, INTENT (IN) :: pvcov(iip1, jjm, ilev)
  REAL, INTENT (IN) :: pucov(iip1, jjp1, ilev)
  REAL, INTENT (IN) :: pteta(iip1, jjp1, ilev)
  REAL ztfi(ilon, ilev)
  REAL, INTENT (IN) :: zplay(ilon, ilev), zplev(ilon, ilev+1)
  INTEGER nbteta
  REAL theta(nbteta)

  ! variable Output

  REAL pvteta(ilon, nbteta)

  ! variables locales

  INTEGER i, j, l, ig0
  REAL ssum
  REAL teta(ilon, ilev)
  REAL ptetau(ip1jmp1, ilev), ptetav(ip1jm, ilev)
  REAL ucovteta(ip1jmp1, ilev), vcovteta(ip1jm, ilev)
  REAL n2(ilon, ilev-1), n2teta(ilon, nbteta)
  REAL ztfiteta(ilon, nbteta)
  REAL rhoteta(ilon, nbteta)
  REAL vorateta(iip1, jjm, nbteta)
  REAL voratetafi(ilon, nbteta), vorpol(iim)


  ! projection teta sur la grille physique

  DO l = 1, llm
    teta(1, l) = pteta(1, 1, l)
    ig0 = 2
    DO j = 2, jjm
      DO i = 1, iim
        teta(ig0, l) = pteta(i, j, l)
        ig0 = ig0 + 1
      END DO
    END DO
    teta(ig0, l) = pteta(1, jjp1, l)
  END DO

  ! calcul pteta sur les grilles U et V

  DO l = 1, llm
    DO j = 1, jjp1
      DO i = 1, iip1
        ig0 = i + (j-1)*iip1
        ptetau(ig0, l) = pteta(i, j, l)
      END DO !i
    END DO !j
    DO j = 1, jjm
      DO i = 1, iip1
        ig0 = i + (j-1)*iip1
        ptetav(ig0, l) = 0.5*(pteta(i,j,l)+pteta(i,j+1,l))
      END DO !i
    END DO !j
  END DO !l

  ! projection pucov, pvcov sur une surface de theta constante

  DO l = 1, nbteta
    ! IM 1rout CALL tetaleveli1j1(ip1jmp1,llm,.true.,ptetau,theta(l),
    CALL tetalevel(ip1jmp1, llm, .TRUE., ptetau, theta(l), pucov, &
      ucovteta(:,l))
    ! IM 1rout CALL tetaleveli1j(ip1jm,llm,.true.,ptetav,theta(l),
    CALL tetalevel(ip1jm, llm, .TRUE., ptetav, theta(l), pvcov, &
      vcovteta(:,l))
  END DO !l

  ! calcul vorticite absolue sur une iso-theta : vorateta

  CALL tourabs(nbteta, vcovteta, ucovteta, vorateta)

  ! projection vorateta sur la grille physique => voratetafi

  DO l = 1, nbteta
    DO j = 2, jjm
      ig0 = 1 + (j-2)*iim
      DO i = 1, iim
        voratetafi(ig0+i+1, l) = vorateta(i, j-1, l)*alpha4_2d(i+1, j) + &
          vorateta(i+1, j-1, l)*alpha1_2d(i+1, j) + &
          vorateta(i, j, l)*alpha3_2d(i+1, j) + vorateta(i+1, j, l)*alpha2_2d &
          (i+1, j)
      END DO
      voratetafi(ig0+1, l) = voratetafi(ig0+1+iim, l)
    END DO
  END DO

  DO l = 1, nbteta
    DO i = 1, iim
      vorpol(i) = vorateta(i, 1, l)*aire_2d(i, 1)
    END DO
    voratetafi(1, l) = ssum(iim, vorpol, 1)/apoln
  END DO

  DO l = 1, nbteta
    DO i = 1, iim
      vorpol(i) = vorateta(i, jjm, l)*aire_2d(i, jjm+1)
    END DO
    voratetafi(ilon, l) = ssum(iim, vorpol, 1)/apols
  END DO

  ! calcul N**2 sur la grille physique => N2

  DO l = 1, llm - 1
    DO i = 1, ilon
      n2(i, l) = (g**2*zplay(i,l)*(ztfi(i,l+1)-ztfi(i, &
        l)))/(r*ztfi(i,l)*ztfi(i,l)*(zplev(i,l)-zplev(i, &
        l+1))) + (g**2)/(ztfi(i,l)*cpp)
    END DO !i
  END DO !l

  ! calcul N2 sur une iso-theta => N2teta

  DO l = 1, nbteta
    CALL tetalevel(ilon, llm-1, .TRUE., teta, theta(l), n2, n2teta(:,l))
    CALL tetalevel(ilon, llm, .TRUE., teta, theta(l), ztfi, ztfiteta(:,l))
  END DO !l=1, nbteta

  ! calcul rho et PV sur une iso-theta : rhoteta, PVteta

  DO l = 1, nbteta
    DO i = 1, ilon
      rhoteta(i, l) = (ztfiteta(i,l)/theta(l))**(cpp/r)*(preff/(r*ztfiteta(i, &
        l)))

      ! PVteta en PVU

      pvteta(i, l) = (theta(l)*g*voratetafi(i,l)*n2teta(i,l))/ &
        (g**2*rhoteta(i,l))

      ! PVteta en 1/(Pa*s)

      pvteta(i, l) = (voratetafi(i,l)*n2teta(i,l))/(g**2*rhoteta(i,l))
    END DO !i
  END DO !l

  RETURN
END SUBROUTINE pvtheta
