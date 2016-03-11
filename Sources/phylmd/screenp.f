module screenp_m

  IMPLICIT none

contains

  ! $Header: /home/cvsroot/LMDZ4/libf/phylmd/screenp.F90,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
  !
  SUBROUTINE screenp(klon, knon, &
       &                   speed, tair, qair, &
       &                   ts, qsurf, rugos, lmon, &
       &                   ustar, testar, qstar, zref, &
       &                   delu, delte, delq) 

    ! Objet : calcul "predicteur" des anomalies du vent, de la temperature 
    !         potentielle et de l'humidite relative au niveau de reference zref et 
    !         par rapport au 1er niveau (pour u) ou a la surface (pour theta et q) 
    !         a partir des relations de Dyer-Businger.
    !
    ! Reference : Hess, Colman et McAvaney (1995)
    !
    ! I. Musat, 01.07.2002
    !-------------------------------------------------------------------------
    !
    ! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
    ! knon----input-I- nombre de points pour un type de surface
    ! speed---input-R- module du vent au 1er niveau du modele
    ! tair----input-R- temperature de l'air au 1er niveau du modele
    ! qair----input-R- humidite relative au 1er niveau du modele
    ! ts------input-R- temperature de l'air a la surface
    ! qsurf---input-R- humidite relative a la surface
    ! rugos---input-R- rugosite
    ! lmon----input-R- longueur de Monin-Obukov
    ! ustar---input-R- facteur d'echelle pour le vent
    ! testar--input-R- facteur d'echelle pour la temperature potentielle
    ! qstar---input-R- facteur d'echelle pour l'humidite relative
    ! zref----input-R- altitude de reference
    !
    ! delu----input-R- anomalie du vent par rapport au 1er niveau
    ! delte---input-R- anomalie de la temperature potentielle par rapport a la surface
    ! delq----input-R- anomalie de l'humidite relative par rapport a la surface
    !
    INTEGER, intent(in) :: klon, knon
    REAL, dimension(klon), intent(in) :: speed, tair, qair
    REAL, dimension(klon), intent(in) :: ts, qsurf, rugos
    DOUBLE PRECISION, dimension(klon), intent(in) :: lmon
    REAL, dimension(klon), intent(in) :: ustar, testar, qstar
    REAL, intent(in) :: zref
    !
    REAL, dimension(klon), intent(out) :: delu, delte, delq
    !
    !-------------------------------------------------------------------------
    ! Variables locales et constantes :
    REAL, PARAMETER :: RKAR=0.40
    INTEGER :: i
    REAL :: xtmp, xtmp0
    !-------------------------------------------------------------------------
    DO i = 1, knon
       !
       IF (lmon(i).GE.0.) THEN
          !
          ! STABLE CASE
          !
          IF (speed(i).GT.1.5.AND.lmon(i).LE.1.0) THEN
             delu(i) = (ustar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) + &
                  min(5d0, 5d0 *(zref - rugos(i))/lmon(i)))
             delte(i) = (testar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) + &
                  min(5d0, 5d0 * (zref - rugos(i))/lmon(i)))
             delq(i) = (qstar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) + &
                  min(5d0, 5d0 * (zref - rugos(i))/lmon(i)))
          ELSE
             delu(i)  = 0.1 * speed(i)
             delte(i) = 0.1 * (tair(i) - ts(i) )
             delq(i)  = 0.1 * (max(qair(i),0.0) - max(qsurf(i),0.0))
          ENDIF
       ELSE  
          !
          ! UNSTABLE CASE
          !
          IF (speed(i).GT.5.0.AND.abs(lmon(i)).LE.50.0) THEN
             xtmp = (1. - 16. * (zref/lmon(i)))**(1./4.)
             xtmp0 = (1. - 16. * (rugos(i)/lmon(i)))**(1./4.)
             delu(i) = (ustar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) & 
                  - 2.*log(0.5*(1. + xtmp)) &
                  + 2.*log(0.5*(1. + xtmp0)) &
                  - log(0.5*(1. + xtmp*xtmp)) &
                  + log(0.5*(1. + xtmp0*xtmp0)) &
                  + 2.*atan(xtmp) - 2.*atan(xtmp0))
             delte(i) = (testar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) &
                  - 2.0 * log(0.5*(1. + xtmp*xtmp)) & 
                  + 2.0 * log(0.5*(1. + xtmp0*xtmp0)))
             delq(i)  = (qstar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) &
                  - 2.0 * log(0.5*(1. + xtmp*xtmp)) & 
                  + 2.0 * log(0.5*(1. + xtmp0*xtmp0)))
          ELSE
             delu(i)  = 0.5 * speed(i)
             delte(i) = 0.5 * (tair(i) - ts(i) )
             delq(i)  = 0.5 * (max(qair(i),0.0) - max(qsurf(i),0.0))
          ENDIF
       ENDIF
       !
    ENDDO
  END SUBROUTINE screenp

end module screenp_m
