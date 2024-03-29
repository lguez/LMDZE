module conf_dat3d_m

  IMPLICIT NONE

contains

  SUBROUTINE conf_dat3d(xd, yd, zd, xf, yf, zf, champd)

    ! From dyn3d/conf_dat3d.F, version 1.1.1.1, 2004/05/19 12:53:05

    ! Author: P. Le Van

    ! Ce sous-programme configure le champ de donn\'ees 3D 'champd' pour
    ! que la longitude varie de - pi \`a pi, la latitude de pi/2 \`a
    ! - pi/2 et pour que la coordonn\'ee pression soit d\'ecroissante.

    use jumble, only: assert_eq, pi

    REAL, intent(in):: xd(:), yd(:) ! longitudes et latitudes initiales, en rad
    REAL, intent(in):: zd(:) ! pressure levels, in Pa or hPa

    REAL, intent(out):: xf(:) ! longitude, in rad, - pi to pi
    REAL, intent(out):: yf(:) ! latitude, in rad, pi/2 to -pi/2
    REAL, intent(out):: zf(:) ! pressure levels, in decreasing order, in Pa
    REAL, intent(inout):: champd(:, :, :)

    ! Variables locales :

    INTEGER lons, lats, levs
    LOGICAL invlon
    REAL rlatmin, rlatmax, oldxd1
    INTEGER i

    !--------------------------------------

    lons = assert_eq(size(xd), size(xf), size(champd, 1), "conf_dat3d lons")
    lats = assert_eq(size(yd), size(yf), size(champd, 2), "conf_dat3d lats")
    levs = assert_eq(size(zd), size(zf), size(champd, 3), "conf_dat3d levs")

    IF (xd(1) >= - pi - 0.5 .AND. xd(lons) <=  pi + 0.5) THEN
       invlon    = .FALSE.
    ELSE IF (xd(1) >= -0.5 .AND. xd(lons) <= 2 * pi+0.5) THEN
       invlon    = .TRUE.
    ELSE
       print *, "Probl\`eme sur les longitudes des donn\'ees"
       stop 1
    ENDIF

    rlatmin = MIN(yd(1), yd(lats))
    rlatmax = MAX(yd(1), yd(lats))

    IF (rlatmin < -pi / 2 - 0.5 .or. rlatmax > pi / 2 + 0.5) THEN
       print *, "Probl\`eme sur les latitudes des donn\'ees"
       stop 1
    ENDIF

    xf(:) = xd(:)
    yf(:) = yd(:)

    IF (invlon) THEN
       ! On tourne les longitudes pour avoir - pi  \`a pi
       DO i=1, lons
          IF (xf(i) >  pi) exit
       ENDDO

       where (xf > pi) xf = xf - 2 * pi
       xf = cshift(xf, i - 1)
       ! On tourne les longitudes pour "champd":
       champd = cshift(champd, i - 1)
    ENDIF

    IF (yd(1) < yd(lats)) THEN
       yf = yf(lats:1:-1)
       champd = champd(:, lats:1:-1, :)
    ENDIF

    oldxd1 = xf(1)
    forall (i = 1: lons-1) xf(i) = 0.5 * (xf(i) + xf(i+1))
    xf(lons) = 0.5 * (xf(lons) + oldxd1 + 2 * pi)
    forall (i = 1: lats-1) yf(i) = 0.5 * (yf(i) + yf(i+1))

    IF (MAX(zd(1), zd(levs)) < 1200.) THEN
       zf(:) = zd(:) * 100. ! convert from hPa to Pa
    else
       zf(:) = zd(:)
    ENDIF

    IF (zd(1) < zd(levs)) THEN
       zf(:) = zf(levs:1:-1)
       champd(:, :, :) = champd(:, :, levs:1:-1)
    ENDIF

  END SUBROUTINE conf_dat3d

end module conf_dat3d_m
