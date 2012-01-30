      SUBROUTINE clvent(knon,dtime, u1lay,v1lay,coef,t,ven,
     e                  paprs,pplay,delp,
     s                  d_ven,flux_v)
      use dimens_m
      use dimphy
      use conf_gcm_m
      use SUPHEC_M
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: diffusion vertical de la vitesse "ven"
c======================================================================
c Arguments:
c dtime----input-R- intervalle du temps (en second)
c u1lay----input-R- vent u de la premiere couche (m/s)
c v1lay----input-R- vent v de la premiere couche (m/s)
c coef-----input-R- le coefficient d'echange (m**2/s) multiplie par
c                   le cisaillement du vent (dV/dz); la premiere
c                   valeur indique la valeur de Cdrag (sans unite)
c t--------input-R- temperature (K)
c ven------input-R- vitesse horizontale (m/s)
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c delp-----input-R- epaisseur de couche (Pa)
c
c
c d_ven----output-R- le changement de "ven"
c flux_v---output-R- (diagnostic) flux du vent: (kg m/s)/(m**2 s)
c======================================================================
      INTEGER knon
      REAL, intent(in):: dtime
      REAL u1lay(klon), v1lay(klon)
      REAL coef(klon,klev)
      REAL t(klon,klev), ven(klon,klev)
      REAL paprs(klon,klev+1), pplay(klon,klev), delp(klon,klev)
      REAL d_ven(klon,klev)
      REAL flux_v(klon,klev)
c======================================================================
c======================================================================
      INTEGER i, k
      REAL zx_cv(klon,2:klev)
      REAL zx_dv(klon,2:klev)
      REAL zx_buf(klon)
      REAL zx_coef(klon,klev)
      REAL local_ven(klon,klev)
      REAL zx_alf1(klon), zx_alf2(klon)
c======================================================================
      DO k = 1, klev
      DO i = 1, knon
         local_ven(i,k) = ven(i,k)
      ENDDO
      ENDDO
c======================================================================
      DO i = 1, knon
ccc         zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
         zx_alf1(i) = 1.0
         zx_alf2(i) = 1.0 - zx_alf1(i)
         zx_coef(i,1) = coef(i,1)
     .                 * (1.0+SQRT(u1lay(i)**2+v1lay(i)**2))
     .                 * pplay(i,1)/(RD*t(i,1))
         zx_coef(i,1) = zx_coef(i,1) * dtime*RG
      ENDDO
c======================================================================
      DO k = 2, klev
      DO i = 1, knon
         zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k))
     .                  *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
         zx_coef(i,k) = zx_coef(i,k) * dtime*RG
      ENDDO
      ENDDO
c======================================================================
      DO i = 1, knon
         zx_buf(i) = delp(i,1) + zx_coef(i,1)*zx_alf1(i)+zx_coef(i,2)
         zx_cv(i,2) = local_ven(i,1)*delp(i,1) / zx_buf(i)
         zx_dv(i,2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1))
     .                /zx_buf(i)
      ENDDO
      DO k = 3, klev
      DO i = 1, knon
         zx_buf(i) = delp(i,k-1) + zx_coef(i,k)
     .                         + zx_coef(i,k-1)*(1.-zx_dv(i,k-1))
         zx_cv(i,k) = (local_ven(i,k-1)*delp(i,k-1)
     .                  +zx_coef(i,k-1)*zx_cv(i,k-1) )/zx_buf(i)
         zx_dv(i,k) = zx_coef(i,k)/zx_buf(i)
      ENDDO
      ENDDO
      DO i = 1, knon
         local_ven(i,klev) = ( local_ven(i,klev)*delp(i,klev)
     .                        +zx_coef(i,klev)*zx_cv(i,klev) )
     .                   / ( delp(i,klev) + zx_coef(i,klev)
     .                       -zx_coef(i,klev)*zx_dv(i,klev) )
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, knon
         local_ven(i,k) = zx_cv(i,k+1) + zx_dv(i,k+1)*local_ven(i,k+1)
      ENDDO
      ENDDO
c======================================================================
c== flux_v est le flux de moment angulaire (positif vers bas)
c== dont l'unite est: (kg m/s)/(m**2 s)
      DO i = 1, knon
         flux_v(i,1) = zx_coef(i,1)/(RG*dtime)
     .                 *(local_ven(i,1)*zx_alf1(i)
     .                  +local_ven(i,2)*zx_alf2(i))
      ENDDO
      DO k = 2, klev
      DO i = 1, knon
         flux_v(i,k) = zx_coef(i,k)/(RG*dtime)
     .               * (local_ven(i,k)-local_ven(i,k-1))
      ENDDO
      ENDDO
c
      DO k = 1, klev
      DO i = 1, knon
         d_ven(i,k) = local_ven(i,k) - ven(i,k)
      ENDDO
      ENDDO
c
      RETURN
      END
