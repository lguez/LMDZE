      SUBROUTINE PVtheta(ilon,ilev,pucov,pvcov,pteta,
     &     ztfi,zplay,zplev,
     &     nbteta,theta,PVteta)
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      IMPLICIT none

c=======================================================================
c
c   Auteur:  I. Musat
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c    Calcul de la vorticite potentielle PVteta sur des iso-theta selon
c    la methodologie du NCEP/NCAR :
c    1) on calcule la stabilite statique N**2=g/T*(dT/dz+g/cp) sur les
c       niveaux du modele => N2
c    2) on interpole les vents, la temperature et le N**2 sur des isentropes
c       (en fait sur des iso-theta) lineairement en log(theta) =>
c       ucovteta, vcovteta, N2teta
c    3) on calcule la vorticite absolue sur des iso-theta => vorateta
c    4) on calcule la densite rho sur des iso-theta => rhoteta 
c
c       rhoteta = (T/theta)**(cp/R)*p0/(R*T)
c
c    5) on calcule la vorticite potentielle sur des iso-theta => PVteta
c
c       PVteta = (vorateta * N2 * theta)/(g * rhoteta) ! en PVU
c
c       NB: 1PVU=10**(-6) K*m**2/(s * kg)
c
c       PVteta =  vorateta * N2/(g**2 * rhoteta) ! en 1/(Pa*s)
c
c
c    *******************************************************************
c
c
c     Variables d'entree : ilon,ilev,pucov,pvcov,pteta,ztfi,zplay,zplev,nbteta,theta
c                       -> sur la grille dynamique
c     Variable de sortie : PVteta
c                       -> sur la grille physique 
c=======================================================================

c
c variables Input
c
      INTEGER ilon
      integer, intent(in):: ilev
      REAL pvcov(iip1,jjm,ilev)
      REAL pucov(iip1,jjp1,ilev)
      REAL pteta(iip1,jjp1,ilev)
      REAL ztfi(ilon,ilev)
      REAL zplay(ilon,ilev), zplev(ilon,ilev+1)
      INTEGER nbteta
      REAL theta(nbteta)
c
c variable Output
c
      REAL PVteta(ilon,nbteta)
c
c variables locales
c
      INTEGER i, j, l, ig0
      REAL SSUM
      REAL teta(ilon, ilev)
      REAL ptetau(ip1jmp1, ilev), ptetav(ip1jm, ilev)
      REAL ucovteta(ip1jmp1,ilev), vcovteta(ip1jm,ilev)
      REAL N2(ilon,ilev-1), N2teta(ilon,nbteta)
      REAL ztfiteta(ilon,nbteta)
      REAL rhoteta(ilon,nbteta)
      REAL vorateta(iip1,jjm,nbteta)
      REAL voratetafi(ilon,nbteta), vorpol(iim)
c
c
c projection teta sur la grille physique
c
      DO l=1,llm
       teta(1,l)   =  pteta(1,1,l)
       ig0         = 2
       DO j = 2, jjm
        DO i = 1, iim
         teta(ig0,l)    = pteta(i,j,l)
         ig0            = ig0 + 1
        ENDDO
       ENDDO
       teta(ig0,l)    = pteta(1,jjp1,l)
      ENDDO
c
c calcul pteta sur les grilles U et V
c
      DO l=1, llm
       DO j=1, jjp1
        DO i=1, iip1
         ig0=i+(j-1)*iip1
         ptetau(ig0,l)=pteta(i,j,l)
        ENDDO !i
       ENDDO !j
       DO j=1, jjm
        DO i=1, iip1
         ig0=i+(j-1)*iip1
         ptetav(ig0,l)=0.5*(pteta(i,j,l)+pteta(i,j+1,l))
        ENDDO !i
       ENDDO !j
      ENDDO !l
c
c projection pucov, pvcov sur une surface de theta constante
c
      DO l=1, nbteta
cIM 1rout CALL tetaleveli1j1(ip1jmp1,llm,.true.,ptetau,theta(l),
       CALL tetalevel(ip1jmp1,llm,.true.,ptetau,theta(l),
     .                pucov,ucovteta(:,l))
cIM 1rout CALL tetaleveli1j(ip1jm,llm,.true.,ptetav,theta(l),
       CALL tetalevel(ip1jm,llm,.true.,ptetav,theta(l),
     .                pvcov,vcovteta(:,l))
      ENDDO !l
c
c calcul vorticite absolue sur une iso-theta : vorateta
c
      CALL tourabs(nbteta,vcovteta,ucovteta,vorateta)
c
c projection vorateta sur la grille physique => voratetafi
c
      DO l=1,nbteta
       DO j=2,jjm
        ig0=1+(j-2)*iim
        DO i=1,iim
           voratetafi(ig0+i+1,l)
     &          = vorateta( i ,j-1,l) * alpha4_2d(i+1,j) +
     &          vorateta(i+1,j-1,l) * alpha1_2d(i+1,j) +
     &          vorateta(i  ,j  ,l) * alpha3_2d(i+1,j) +
     &          vorateta(i+1,j  ,l) * alpha2_2d(i+1,j)
        ENDDO
        voratetafi(ig0 +1,l) = voratetafi(ig0 +1+ iim,l)
       ENDDO
      ENDDO
c
      DO l=1,nbteta
       DO i=1,iim
        vorpol(i)  = vorateta(i,1,l)*aire_2d(i,1)
       ENDDO
       voratetafi(1,l)= SSUM(iim,vorpol,1)/apoln
      ENDDO
c
      DO l=1,nbteta
       DO i=1,iim
        vorpol(i)  = vorateta(i,jjm,l)* aire_2d(i,jjm +1)
       ENDDO
       voratetafi(ilon,l)= SSUM(iim,vorpol,1)/apols
      ENDDO
c 
c calcul N**2 sur la grille physique => N2
c
      DO l=1, llm-1 
       DO i=1, ilon
        N2(i,l) = (g**2 * zplay(i,l) * 
     &            (ztfi(i,l+1)-ztfi(i,l)) )/
     &            (R*ztfi(i,l)*ztfi(i,l)*
     &            (zplev(i,l)-zplev(i,l+1)) )+
     &            (g**2)/(ztfi(i,l)*CPP)
       ENDDO !i
      ENDDO !l
c
c calcul N2 sur une iso-theta => N2teta 
c
      DO l=1, nbteta
       CALL tetalevel(ilon,llm-1,.true.,teta,theta(l),
     &                N2,N2teta(:,l))
       CALL tetalevel(ilon,llm,.true.,teta,theta(l),
     &                ztfi,ztfiteta(:,l))
      ENDDO !l=1, nbteta
c
c calcul rho et PV sur une iso-theta : rhoteta, PVteta
c
      DO l=1, nbteta
       DO i=1, ilon
        rhoteta(i,l)=(ztfiteta(i,l)/theta(l))**(CPP/R)*
     &  (preff/(R*ztfiteta(i,l)))
c
c PVteta en PVU
c
        PVteta(i,l)=(theta(l)*g*voratetafi(i,l)*N2teta(i,l))/
     &              (g**2*rhoteta(i,l))
c
c PVteta en 1/(Pa*s)
c
        PVteta(i,l)=(voratetafi(i,l)*N2teta(i,l))/
     &              (g**2*rhoteta(i,l))
       ENDDO !i
      ENDDO !l
c
      RETURN
      END 
