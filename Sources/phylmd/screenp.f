module screenp_m

  IMPLICIT none

contains

  SUBROUTINE screenp(knon, speed, tair, qair, ts, qsurf, rugos, lmon, ustar, &
       testar, qstar, zref, delu, delte, delq)

    ! From LMDZ4/libf/phylmd/screenp.F90, version 1.1.1.1, 2004/05/19 12:53:09

    ! Objet : calcul "pr\'edicteur" des anomalies du vent, de la
    ! temp\'erature potentielle et de l'humidit\'e relative au niveau
    ! de r\'ef\'erence zref et par rapport au 1er niveau (pour u) ou
    ! \`a la surface (pour theta et q) \`a partir des relations de
    ! Dyer-Businger.

    ! Reference: Hess, Colman et McAvaney (1995)

    ! I. Musat, 01.07.2002

    use dimphy, only: klon

    INTEGER, intent(in):: knon
    ! knon----input-I- nombre de points pour un type de surface
    REAL, dimension(klon), intent(in):: speed, tair, qair
    ! speed---input-R- module du vent au 1er niveau du modele
    ! tair----input-R- temperature de l'air au 1er niveau du modele
    ! qair----input-R- humidite relative au 1er niveau du modele
    REAL, dimension(klon), intent(in):: ts, qsurf, rugos
    ! ts------input-R- temperature de l'air a la surface
    ! qsurf---input-R- humidite relative a la surface
    ! rugos---input-R- rugosite
    DOUBLE PRECISION, dimension(klon), intent(in):: lmon
    ! lmon----input-R- longueur de Monin-Obukov
    REAL, dimension(klon), intent(in):: ustar, testar, qstar
    ! ustar---input-R- facteur d'echelle pour le vent
    ! testar--input-R- facteur d'echelle pour la temperature potentielle
    ! qstar---input-R- facteur d'echelle pour l'humidite relative
    REAL, intent(in):: zref
    ! zref----input-R- altitude de reference

    REAL, dimension(klon), intent(out):: delu
    ! delu----input-R- anomalie du vent par rapport au 1er niveau

    REAL, dimension(klon), intent(out):: delte
    ! anomalie de la temperature potentielle par rapport a la surface

    REAL, dimension(klon), intent(out):: delq
    ! delq----input-R- anomalie de l'humidite relative par rapport a la surface

    ! Local:
    REAL, PARAMETER:: RKAR=0.40
    INTEGER i
    REAL xtmp, xtmp0
    
    !-------------------------------------------------------------------------
    
    DO i = 1, knon
       IF (lmon(i) >= 0.) THEN
          ! STABLE CASE
          IF (speed(i) > 1.5.AND.lmon(i) <= 1.0) THEN
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
             delu(i) = 0.1 * speed(i)
             delte(i) = 0.1 * (tair(i) - ts(i))
             delq(i) = 0.1 * (max(qair(i), 0.0) - max(qsurf(i), 0.0))
          ENDIF
       ELSE
          ! UNSTABLE CASE
          IF (speed(i) > 5.0.AND.abs(lmon(i)) <= 50.0) THEN
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
             delq(i) = (qstar(i)/RKAR)* &
                  (log(zref/(rugos(i))+1.) &
                  - 2.0 * log(0.5*(1. + xtmp*xtmp)) &
                  + 2.0 * log(0.5*(1. + xtmp0*xtmp0)))
          ELSE
             delu(i) = 0.5 * speed(i)
             delte(i) = 0.5 * (tair(i) - ts(i))
             delq(i) = 0.5 * (max(qair(i), 0.0) - max(qsurf(i), 0.0))
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE screenp

end module screenp_m
