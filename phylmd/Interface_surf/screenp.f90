module screenp_m

  IMPLICIT none

contains

  SUBROUTINE screenp(speed, tair, qair, ts, qsurf, rugos, lmon, ustar, testar, &
       qstar, zref, delu, delte, delq)

    ! From LMDZ4/libf/phylmd/screenp.F90, version 1.1.1.1, 2004/05/19 12:53:09

    ! Objet : calcul "pr\'edicteur" des anomalies du vent, de la
    ! temp\'erature potentielle et de l'humidit\'e relative au niveau
    ! de r\'ef\'erence zref et par rapport au 1er niveau (pour u) ou
    ! \`a la surface (pour theta et q) \`a partir des relations de
    ! Dyer-Businger.

    ! Reference: Hess, Colman and McAvaney (1995)

    ! I. Musat, July 2002

    REAL, intent(in):: speed(:) ! (knon) module du vent au 1er niveau du modele

    REAL, intent(in):: tair(:) ! (knon)
    ! temperature de l'air au 1er niveau du modele

    REAL, intent(in):: qair(:) ! (knon)
    ! humidite relative au 1er niveau du modele

    REAL, intent(in):: ts(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidite relative a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosite
    DOUBLE PRECISION, intent(in):: lmon(:) ! (knon) longueur de Monin-Obukov
    REAL, intent(in):: ustar(:) ! (knon) facteur d'\'echelle pour le vent

    REAL, intent(in):: testar(:) ! (knon)
    ! facteur d'echelle pour la temperature potentielle

    REAL, intent(in):: qstar(:) ! (knon)
    ! facteur d'echelle pour l'humidite relative

    REAL, intent(in):: zref ! altitude de reference

    REAL, intent(out):: delu(:) ! (knon)
    ! anomalie du vent par rapport au 1er niveau

    REAL, intent(out):: delte(:) ! (knon)
    ! anomalie de la temperature potentielle par rapport a la surface

    REAL, intent(out):: delq(:) ! (knon)
    ! anomalie de l'humidite relative par rapport a la surface

    ! Local:
    REAL, PARAMETER:: RKAR = 0.4
    INTEGER i
    REAL xtmp, xtmp0

    !-------------------------------------------------------------------------

    DO i = 1, size(speed)
       IF (lmon(i) >= 0.) THEN
          ! Stable case
          IF (speed(i) > 1.5 .AND. lmon(i) <= 1.) THEN
             delu(i) = (ustar(i) / RKAR) * (log(zref / (rugos(i)) + 1.) &
                  + min(5d0, 5d0 * (zref - rugos(i)) / lmon(i)))
             delte(i) = (testar(i) / RKAR) * (log(zref / (rugos(i)) + 1.) &
                  + min(5d0, 5d0 * (zref - rugos(i)) / lmon(i)))
             delq(i) = (qstar(i) / RKAR) * (log(zref / (rugos(i)) + 1.) &
                  + min(5d0, 5d0 * (zref - rugos(i)) / lmon(i)))
          ELSE
             delu(i) = 0.1 * speed(i)
             delte(i) = 0.1 * (tair(i) - ts(i))
             delq(i) = 0.1 * (max(qair(i), 0.) - max(qsurf(i), 0.))
          ENDIF
       ELSE
          ! Unstable case
          IF (speed(i) > 5. .AND. abs(lmon(i)) <= 50.) THEN
             xtmp = (1. - 16. * (zref / lmon(i)))**(1. / 4.)
             xtmp0 = (1. - 16. * (rugos(i) / lmon(i)))**(1. / 4.)
             delu(i) = (ustar(i) / RKAR) * (log(zref / (rugos(i)) + 1.) - 2. &
                  * log(0.5 * (1. + xtmp)) + 2. * log(0.5 * (1. + xtmp0)) &
                  - log(0.5 * (1. + xtmp * xtmp)) + log(0.5 * (1. + xtmp0 &
                  * xtmp0)) + 2. * atan(xtmp) - 2. * atan(xtmp0))
             delte(i) = (testar(i) / RKAR) * (log(zref / (rugos(i)) + 1.) - 2. &
                  * log(0.5 * (1. + xtmp * xtmp)) + 2. * log(0.5 * (1. + xtmp0 &
                  * xtmp0)))
             delq(i) = (qstar(i) / RKAR) * (log(zref / (rugos(i)) + 1.) - 2. &
                  * log(0.5 * (1. + xtmp * xtmp)) + 2. * log(0.5 * (1. + xtmp0 &
                  * xtmp0)))
          ELSE
             delu(i) = 0.5 * speed(i)
             delte(i) = 0.5 * (tair(i) - ts(i))
             delq(i) = 0.5 * (max(qair(i), 0.) - max(qsurf(i), 0.))
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE screenp

end module screenp_m
