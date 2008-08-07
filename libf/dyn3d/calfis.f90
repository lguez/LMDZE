module calfis_m

  ! Clean: no C preprocessor directive, no include line

  IMPLICIT NONE

contains

  SUBROUTINE calfis(nq, lafin, rdayvrai, heure, pucov, pvcov, pteta, pq, &
       pmasse, pps, ppk, pphis, pphi, pducov, pdvcov, pdteta, pdq, pw, &
       pdufi, pdvfi, pdhfi, pdqfi, pdpsfi)

    ! From dyn3d/calfis.F,v 1.3 2005/05/25 13:10:09

    ! Auteurs : P. Le Van, F. Hourdin

    !   1. rearrangement des tableaux et transformation
    !      variables dynamiques  >  variables physiques
    !   2. calcul des termes physiques
    !   3. retransformation des tendances physiques en tendances dynamiques

    !   remarques:
    !   ----------

    !    - les vents sont donnes dans la physique par leurs composantes 
    !      naturelles.
    !    - la variable thermodynamique de la physique est une variable
    !      intensive :   T 
    !      pour la dynamique on prend    T * (preff / p(l)) **kappa
    !    - les deux seules variables dependant de la geometrie necessaires
    !      pour la physique sont la latitude pour le rayonnement et 
    !      l'aire de la maille quand on veut integrer une grandeur 
    !      horizontalement.

    !     Input :
    !     -------
    !       pucov           covariant zonal velocity
    !       pvcov           covariant meridional velocity 
    !       pteta           potential temperature
    !       pps             surface pressure
    !       pmasse          masse d'air dans chaque maille
    !       pts             surface temperature  (K)
    !       callrad         clef d'appel au rayonnement

    !    Output :
    !    --------
    !        pdufi          tendency for the natural zonal velocity (ms-1)
    !        pdvfi          tendency for the natural meridional velocity 
    !        pdhfi          tendency for the potential temperature
    !        pdtsfi         tendency for the surface temperature

    !        pdtrad         radiative tendencies  \  both input
    !        pfluxrad       radiative fluxes      /  and output

    use dimens_m, only: iim, jjm, llm, nqmx
    use dimphy, only: klon
    use comconst, only: kappa, cpp, dtphys, g, pi
    use comvert, only: preff, presnivs
    use comgeom, only: apoln, cu_2d, cv_2d, unsaire_2d, apols, rlonu, rlonv
    use iniadvtrac_m, only: niadv
    use grid_change, only: dyn_phy, gr_fi_dyn
    use physiq_m, only: physiq
    use pressure_var, only: p3d, pls

    !    0.  Declarations :

    INTEGER, intent(in):: nq

    !    Arguments :

    LOGICAL, intent(in):: lafin
    REAL, intent(in):: heure ! heure de la journée en fraction de jour

    REAL pvcov(iim + 1,jjm,llm)
    REAL pucov(iim + 1,jjm + 1,llm)
    REAL pteta(iim + 1,jjm + 1,llm)
    REAL pmasse(iim + 1,jjm + 1,llm)

    REAL, intent(in):: pq(iim + 1,jjm + 1,llm,nqmx)
    ! (mass fractions of advected fields)

    REAL pphis(iim + 1,jjm + 1)
    REAL pphi(iim + 1,jjm + 1,llm)

    REAL pdvcov(iim + 1,jjm,llm)
    REAL pducov(iim + 1,jjm + 1,llm)
    REAL pdteta(iim + 1,jjm + 1,llm)
    REAL pdq(iim + 1,jjm + 1,llm,nqmx)

    REAL pw(iim + 1,jjm + 1,llm)

    REAL pps(iim + 1,jjm + 1)
    REAL, intent(in):: ppk(iim + 1,jjm + 1,llm)

    REAL pdvfi(iim + 1,jjm,llm)
    REAL pdufi(iim + 1,jjm + 1,llm)
    REAL pdhfi(iim + 1,jjm + 1,llm)
    REAL pdqfi(iim + 1,jjm + 1,llm,nqmx)
    REAL pdpsfi(iim + 1,jjm + 1)

    INTEGER, PARAMETER:: longcles = 20

    !    Local variables :

    INTEGER i,j,l,ig0,ig,iq,iiq
    REAL zpsrf(klon)
    REAL zplev(klon,llm+1),zplay(klon,llm)
    REAL zphi(klon,llm),zphis(klon)

    REAL zufi(klon,llm), zvfi(klon,llm)
    REAL ztfi(klon,llm) ! temperature
    real zqfi(klon,llm,nqmx) ! mass fractions of advected fields

    REAL pcvgu(klon,llm), pcvgv(klon,llm)
    REAL pcvgt(klon,llm), pcvgq(klon,llm,2)

    REAL pvervel(klon,llm)

    REAL zdufi(klon,llm),zdvfi(klon,llm)
    REAL zdtfi(klon,llm),zdqfi(klon,llm,nqmx)
    REAL zdpsrf(klon)

    REAL zsin(iim),zcos(iim),z1(iim)
    REAL zsinbis(iim),zcosbis(iim),z1bis(iim)
    REAL pksurcp(iim + 1,jjm + 1)

    ! I. Musat: diagnostic PVteta, Amip2
    INTEGER, PARAMETER:: ntetaSTD=3
    REAL:: rtetaSTD(ntetaSTD) = (/350., 380., 405./)
    REAL PVteta(klon,ntetaSTD)

    REAL SSUM

    LOGICAL:: firstcal = .true.
    REAL, intent(in):: rdayvrai

    !-----------------------------------------------------------------------

    !!print *, "Call sequence information: calfis"

    !    1. Initialisations :
    !   latitude, longitude et aires des mailles pour la physique:

    !   40. transformation des variables dynamiques en variables physiques:
    !   41. pressions au sol (en Pascals)

    zpsrf(1) = pps(1,1)

    ig0  = 2
    DO j = 2,jjm
       CALL SCOPY(iim,pps(1,j),1,zpsrf(ig0), 1)
       ig0 = ig0+iim
    ENDDO

    zpsrf(klon) = pps(1,jjm + 1)

    !   42. pression intercouches :

    !     .... zplev  definis aux (llm +1) interfaces des couches  ....
    !     .... zplay  definis aux (llm)    milieux des couches  .... 

    !    ...    Exner = cp * (p(l) / preff) ** kappa     ....

    forall (l = 1: llm+1) zplev(:, l) = pack(p3d(:, :, l), dyn_phy)

    !   43. temperature naturelle (en K) et pressions milieux couches .
    DO l=1,llm
       pksurcp     =  ppk(:, :, l) / cpp
       pls(:, :, l) = preff * pksurcp**(1./ kappa)
       zplay(:, l) = pack(pls(:, :, l), dyn_phy)
       ztfi(:, l) = pack(pteta(:, :, l) * pksurcp, dyn_phy)
       pcvgt(:, l) = pack(pdteta(:, :, l) * pksurcp / pmasse(:, :, l), dyn_phy)
    ENDDO

    !   43.bis traceurs

    DO iq=1,nq
       iiq=niadv(iq) 
       DO l=1,llm
          zqfi(1,l,iq) = pq(1,1,l,iiq)
          ig0          = 2
          DO j=2,jjm
             DO i = 1, iim
                zqfi(ig0,l,iq)  = pq(i,j,l,iiq)
                ig0             = ig0 + 1
             ENDDO
          ENDDO
          zqfi(ig0,l,iq) = pq(1,jjm + 1,l,iiq)
       ENDDO
    ENDDO

    !   convergence dynamique pour les traceurs "EAU"

    DO iq=1,2
       DO l=1,llm
          pcvgq(1,l,iq)= pdq(1,1,l,iq) / pmasse(1,1,l)
          ig0          = 2
          DO j=2,jjm
             DO i = 1, iim
                pcvgq(ig0,l,iq) = pdq(i,j,l,iq) / pmasse(i,j,l)
                ig0             = ig0 + 1
             ENDDO
          ENDDO
          pcvgq(ig0,l,iq)= pdq(1,jjm + 1,l,iq) / pmasse(1,jjm + 1,l)
       ENDDO
    ENDDO

    !   Geopotentiel calcule par rapport a la surface locale:

    forall (l = 1:llm) zphi(:, l) = pack(pphi(:, :, l), dyn_phy)
    zphis = pack(pphis, dyn_phy)
    DO l=1,llm
       DO ig=1,klon
          zphi(ig,l)=zphi(ig,l)-zphis(ig)
       ENDDO
    ENDDO

    !   ....  Calcul de la vitesse  verticale  (en Pa*m*s  ou Kg/s)  ....

    DO l=1,llm
       pvervel(1,l)=pw(1,1,l) * g /apoln
       ig0=2
       DO j=2,jjm
          DO i = 1, iim
             pvervel(ig0,l) = pw(i,j,l) * g * unsaire_2d(i,j)
             ig0 = ig0 + 1
          ENDDO
       ENDDO
       pvervel(ig0,l)=pw(1,jjm + 1,l) * g /apols
    ENDDO

    !   45. champ u:

    DO  l=1,llm

       DO  j=2,jjm
          ig0 = 1+(j-2)*iim
          zufi(ig0+1,l)= 0.5 *  &
               (pucov(iim,j,l)/cu_2d(iim,j) + pucov(1,j,l)/cu_2d(1,j))
          pcvgu(ig0+1,l)= 0.5 *  &
               (pducov(iim,j,l)/cu_2d(iim,j) + pducov(1,j,l)/cu_2d(1,j))
          DO i=2,iim
             zufi(ig0+i,l)= 0.5 * &
                  (pucov(i-1,j,l)/cu_2d(i-1,j) &
                  + pucov(i,j,l)/cu_2d(i,j))
             pcvgu(ig0+i,l)= 0.5 * &
                  (pducov(i-1,j,l)/cu_2d(i-1,j) &
                  + pducov(i,j,l)/cu_2d(i,j))
          end DO
       end DO

    end DO

    !   46.champ v:

    DO l=1,llm
       DO j=2,jjm
          ig0=1+(j-2)*iim
          DO i=1,iim
             zvfi(ig0+i,l)= 0.5 * &
                  (pvcov(i,j-1,l)/cv_2d(i,j-1) &
                  + pvcov(i,j,l)/cv_2d(i,j))
             pcvgv(ig0+i,l)= 0.5 * &
                  (pdvcov(i,j-1,l)/cv_2d(i,j-1) &
                  + pdvcov(i,j,l)/cv_2d(i,j))
          ENDDO
       ENDDO
    ENDDO

    !   47. champs de vents aux pole nord   
    !        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
    !        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

    DO l=1,llm

       z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,1,l)/cv_2d(1,1)
       z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,1,l)/cv_2d(1,1)
       DO i=2,iim
          z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,1,l)/cv_2d(i,1)
          z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,1,l)/cv_2d(i,1)
       ENDDO

       DO i=1,iim
          zcos(i)   = COS(rlonv(i))*z1(i)
          zcosbis(i)= COS(rlonv(i))*z1bis(i)
          zsin(i)   = SIN(rlonv(i))*z1(i)
          zsinbis(i)= SIN(rlonv(i))*z1bis(i)
       ENDDO

       zufi(1,l)  = SSUM(iim,zcos,1)/pi
       pcvgu(1,l) = SSUM(iim,zcosbis,1)/pi
       zvfi(1,l)  = SSUM(iim,zsin,1)/pi
       pcvgv(1,l) = SSUM(iim,zsinbis,1)/pi

    ENDDO

    !   48. champs de vents aux pole sud:
    !        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
    !        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

    DO l=1,llm

       z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,jjm,l) &
            /cv_2d(1,jjm)
       z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,jjm,l) &
            /cv_2d(1,jjm)
       DO i=2,iim
          z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,jjm,l)/cv_2d(i,jjm)
          z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,jjm,l)/cv_2d(i,jjm)
       ENDDO

       DO i=1,iim
          zcos(i)    = COS(rlonv(i))*z1(i)
          zcosbis(i) = COS(rlonv(i))*z1bis(i)
          zsin(i)    = SIN(rlonv(i))*z1(i)
          zsinbis(i) = SIN(rlonv(i))*z1bis(i)
       ENDDO

       zufi(klon,l)  = SSUM(iim,zcos,1)/pi
       pcvgu(klon,l) = SSUM(iim,zcosbis,1)/pi
       zvfi(klon,l)  = SSUM(iim,zsin,1)/pi
       pcvgv(klon,l) = SSUM(iim,zsinbis,1)/pi

    ENDDO

    !IM calcul PV a teta=350, 380, 405K
    CALL PVtheta(klon,llm,pucov,pvcov,pteta, &
         ztfi,zplay,zplev, &
         ntetaSTD,rtetaSTD,PVteta)

    !   Appel de la physique:

    CALL physiq(nq, firstcal, lafin, rdayvrai, heure, dtphys, &
         zplev, zplay, zphi, zphis, presnivs, zufi, zvfi, &
         ztfi, zqfi, pvervel, zdufi, zdvfi, zdtfi, zdqfi, zdpsrf, pducov, &
         PVteta) ! IM diagnostique PVteta, Amip2

    !   transformation des tendances physiques en tendances dynamiques:

    !  tendance sur la pression :

    pdpsfi = gr_fi_dyn(zdpsrf)

    !   62. enthalpie potentielle

    DO l=1,llm

       DO i=1,iim + 1
          pdhfi(i,1,l)    = cpp *  zdtfi(1,l)      / ppk(i, 1  ,l)
          pdhfi(i,jjm + 1,l) = cpp *  zdtfi(klon,l)/ ppk(i,jjm + 1,l)
       ENDDO

       DO j=2,jjm
          ig0=1+(j-2)*iim
          DO i=1,iim
             pdhfi(i,j,l) = cpp * zdtfi(ig0+i,l) / ppk(i,j,l)
          ENDDO
          pdhfi(iim + 1,j,l) =  pdhfi(1,j,l)
       ENDDO

    ENDDO

    !   62. humidite specifique

    DO iq=1,nqmx
       DO l=1,llm
          DO i=1,iim + 1
             pdqfi(i,1,l,iq)    = zdqfi(1,l,iq)
             pdqfi(i,jjm + 1,l,iq) = zdqfi(klon,l,iq)
          ENDDO
          DO j=2,jjm
             ig0=1+(j-2)*iim
             DO i=1,iim
                pdqfi(i,j,l,iq) = zdqfi(ig0+i,l,iq)
             ENDDO
             pdqfi(iim + 1,j,l,iq) = pdqfi(1,j,l,iq)
          ENDDO
       ENDDO
    ENDDO

    !   63. traceurs

    !     initialisation des tendances
    pdqfi=0.

    DO iq=1,nq
       iiq=niadv(iq)
       DO l=1,llm
          DO i=1,iim + 1
             pdqfi(i,1,l,iiq)    = zdqfi(1,l,iq)
             pdqfi(i,jjm + 1,l,iiq) = zdqfi(klon,l,iq)
          ENDDO
          DO j=2,jjm
             ig0=1+(j-2)*iim
             DO i=1,iim
                pdqfi(i,j,l,iiq) = zdqfi(ig0+i,l,iq)
             ENDDO
             pdqfi(iim + 1,j,l,iiq) = pdqfi(1,j,l,iq)
          ENDDO
       ENDDO
    ENDDO

    !   65. champ u:

    DO l=1,llm

       DO i=1,iim + 1
          pdufi(i,1,l)    = 0.
          pdufi(i,jjm + 1,l) = 0.
       ENDDO

       DO j=2,jjm
          ig0=1+(j-2)*iim
          DO i=1,iim-1
             pdufi(i,j,l)= &
                  0.5*(zdufi(ig0+i,l)+zdufi(ig0+i+1,l))*cu_2d(i,j)
          ENDDO
          pdufi(iim,j,l)= &
               0.5*(zdufi(ig0+1,l)+zdufi(ig0+iim,l))*cu_2d(iim,j)
          pdufi(iim + 1,j,l)=pdufi(1,j,l)
       ENDDO

    ENDDO

    !   67. champ v:

    DO l=1,llm

       DO j=2,jjm-1
          ig0=1+(j-2)*iim
          DO i=1,iim
             pdvfi(i,j,l)= &
                  0.5*(zdvfi(ig0+i,l)+zdvfi(ig0+i+iim,l))*cv_2d(i,j)
          ENDDO
          pdvfi(iim + 1,j,l) = pdvfi(1,j,l)
       ENDDO
    ENDDO

    !   68. champ v pres des poles:
    !      v = U * cos(long) + V * SIN(long)

    DO l=1,llm

       DO i=1,iim
          pdvfi(i,1,l)= &
               zdufi(1,l)*COS(rlonv(i))+zdvfi(1,l)*SIN(rlonv(i))
          pdvfi(i,jjm,l)=zdufi(klon,l)*COS(rlonv(i)) &
               +zdvfi(klon,l)*SIN(rlonv(i))
          pdvfi(i,1,l)= &
               0.5*(pdvfi(i,1,l)+zdvfi(i+1,l))*cv_2d(i,1)
          pdvfi(i,jjm,l)= &
               0.5*(pdvfi(i,jjm,l)+zdvfi(klon-iim-1+i,l))*cv_2d(i,jjm)
       ENDDO

       pdvfi(iim + 1,1,l)  = pdvfi(1,1,l)
       pdvfi(iim + 1,jjm,l)= pdvfi(1,jjm,l)

    ENDDO

    firstcal = .FALSE.

  END SUBROUTINE calfis

end module calfis_m
