      SUBROUTINE coefkz2(nsrf, knon, paprs, pplay,t,
     .                  pcfm, pcfh)
      use dimens_m
      use indicesol
      use dimphy
      use iniprint
      use SUPHEC_M
      IMPLICIT none
c======================================================================
c J'introduit un peu de diffusion sauf dans les endroits
c ou une forte inversion est presente
c On peut dire qu'il represente la convection peu profonde
c
c Arguments:
c nsrf-----input-I- indicateur de la nature du sol
c knon-----input-I- nombre de points a traiter
c paprs----input-R- pression a chaque intercouche (en Pa)
c pplay----input-R- pression au milieu de chaque couche (en Pa)
c t--------input-R- temperature (K)
c
c pcfm-----output-R- coefficients a calculer (vitesse)
c pcfh-----output-R- coefficients a calculer (chaleur et humidite)
c======================================================================
c
c Arguments:
c
      INTEGER knon, nsrf
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL t(klon,klev)
c
      REAL pcfm(klon,klev), pcfh(klon,klev)
c
c Quelques constantes et options:
c
      REAL prandtl
      PARAMETER (prandtl=0.4)
      REAL kstable
      PARAMETER (kstable=0.002)
ccc      PARAMETER (kstable=0.001)
      REAL mixlen ! constante controlant longueur de melange
      PARAMETER (mixlen=35.0)
      REAL seuil ! au-dela l'inversion est consideree trop faible
      PARAMETER (seuil=-0.02)
ccc      PARAMETER (seuil=-0.04)
ccc      PARAMETER (seuil=-0.06)
ccc      PARAMETER (seuil=-0.09)
c
c Variables locales:
c
      INTEGER i, k, invb(knon)
      REAL zl2(knon)
      REAL zdthmin(knon), zdthdp
c
c Initialiser les sorties
c
      DO k = 1, klev
      DO i = 1, knon
         pcfm(i,k) = 0.0
         pcfh(i,k) = 0.0
      ENDDO
      ENDDO
c
c Chercher la zone d'inversion forte
c
      DO i = 1, knon
         invb(i) = klev
         zdthmin(i)=0.0
      ENDDO
      DO k = 2, klev/2-1
      DO i = 1, knon
         zdthdp = (t(i,k)-t(i,k+1))/(pplay(i,k)-pplay(i,k+1))
     .          - RD * 0.5*(t(i,k)+t(i,k+1))/RCPD/paprs(i,k+1)
         zdthdp = zdthdp * 100.0
         IF (pplay(i,k).GT.0.8*paprs(i,1) .AND.
     .       zdthdp.LT.zdthmin(i) ) THEN
            zdthmin(i) = zdthdp
            invb(i) = k
         ENDIF
      ENDDO
      ENDDO
c
c Introduire une diffusion:
c
      DO k = 2, klev
      DO i = 1, knon
cIM cf FH/GK   IF ( (nsrf.NE.is_oce) .OR.  ! si ce n'est pas sur l'ocean
cIM cf FH/GK  .     (invb(i).EQ.klev) .OR. ! s'il n'y a pas d'inversion
      !IM cf JLD/ GKtest TERkz2
      ! IF ( (nsrf.EQ.is_ter) .OR.  ! si on est sur la terre
      ! fin GKtest
      IF ( (nsrf.EQ.is_oce) .AND.  ! si on est sur ocean et si 
     .     ( (invb(i).EQ.klev) .OR.      ! s'il n'y a pas d'inversion
     .     (zdthmin(i).GT.seuil) ) )THEN ! si l'inversion est trop faible
         zl2(i)=(mixlen*MAX(0.0,(paprs(i,k)-paprs(i,klev+1))
     .                       /(paprs(i,2)-paprs(i,klev+1)) ))**2
         pcfm(i,k)= zl2(i)* kstable
         pcfh(i,k) = pcfm(i,k) /prandtl ! h et m different
      ENDIF
      ENDDO
      ENDDO
c
      RETURN
      END
