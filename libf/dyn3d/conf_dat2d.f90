module conf_dat2d_m

  ! From conf_dat2d.F, version 1.2 2006/01/27 15:14:22

  IMPLICIT NONE

contains

  SUBROUTINE conf_dat2d(xd, yd, xf, yf, champd, interbar)

    ! Auteur : P. Le Van

    ! Ce sous-programme configure le champ de données 2D 'champd' et
    ! les longitudes et latitudes de telle façon qu'on ait - pi à pi
    ! en longitude et pi/2 à - pi/2 en latitude.
    
    ! This procedure receives a 2D field, with the corresponding
    ! coordinate variables: longitude and latitude.
    ! The procedure converts longitude and latitude to radians, if the
    ! input values are in degrees.
    ! If the input longitudes are between 0 and 2 pi, the procedure
    ! computes the congruent longitudes between -pi and pi, and permutes
    ! them so they stay in increasing order.
    ! If the input latitudes are from south pole to north pole, the
    ! procedure permutes them so they become from north to south.
    ! Any change on longitudes or latitudes induces a change on the 2D field.
    ! If required, the longitudes and latitudes are finally replaced
    ! by their mid-values.

    use nr_util, only: assert_eq
    use comconst, only: pi

    REAL, intent(in):: xd(:)
    ! (longitudes, in degrees or radians, in increasing order, from 0°
    ! to 360° or -180° to 180°)

    REAL, intent(in):: yd(:)
    ! (latitudes, in degrees or radians, in increasing or decreasing
    ! order, from pole to pole)

    LOGICAL, intent(in), optional:: interbar
    REAL, intent(out):: xf(:), yf(:) ! longitudes and latitudes, in rad
    REAL, intent(inout):: champd(:, :)

    ! Variables locales:

    INTEGER  lons, lats
    LOGICAL radianlon ! "xd" is in degrees
    logical invlon ! "xd" contains longitudes between 0 and 2 pi
    logical radianlat ! "yd" is in rad
    REAL rlatmin, rlatmax, old_xf_1
    INTEGER i, j
    logical mid_values

    !------------------------------

    lons = assert_eq(size(xd), size(xf), size(champd, 1), "conf_dat2d lons")
    lats = assert_eq(size(yd), size(yf), size(champd, 2), "conf_dat2d lats")

    IF (xd(1) >= - pi -0.5 .AND. xd(lons) <= pi + 0.5) THEN
       radianlon = .TRUE.
       invlon    = .FALSE.
    ELSE IF (xd(1) >= -0.5 .AND. xd(lons) <= 2 * pi+0.5) THEN
       radianlon = .TRUE.
       invlon    = .TRUE.
    ELSE IF (xd(1) >= -180.5 .AND. xd(lons) <= 180.5) THEN
       radianlon = .FALSE.
       invlon    = .FALSE.
    ELSE IF (xd(1) >= -0.5 .AND. xd(lons) <= 360.5) THEN
       radianlon = .FALSE.
       invlon    = .TRUE.
    ELSE
       stop '"conf_dat2d": problem with longitudes'
    ENDIF

    rlatmin = MIN(yd(1), yd(lats))
    rlatmax = MAX(yd(1), yd(lats))

    IF (rlatmin >= -pi / 2 - 0.5 .AND. rlatmax <= pi / 2 + 0.5)THEN
       radianlat = .TRUE.
    ELSE IF (rlatmin >= -90. - 0.5 .AND. rlatmax <= 90. + 0.5) THEN
       radianlat = .FALSE.
    ELSE
       stop '"conf_dat2d": problem with latitudes'
    ENDIF

    IF (radianlon)  THEN
       xf(:) = xd(:)
    else
       xf(:) = xd(:) * pi / 180. ! convert to rad
    ENDIF

    IF (radianlat)  THEN
       yf(:) = yd(:)
    else
       yf(:) = yd(:) * pi / 180. ! convert to rad
    ENDIF

    IF (invlon) THEN
       ! On tourne les longitudes pour avoir - pi à + pi :

       ! Get the index of the first longitude > pi:
       i = 1
       do while (xf(i) <= pi)
          i = i + 1
       end do

       xf(i:) = xf(i:) - 2 * pi
       xf(:) = cshift(xf, shift=i - 1)
       champd(:, :) = cshift(champd, shift=i - 1)
    ENDIF

    IF (yd(1) < yd(lats)) THEN
       ! "yd" contains latitudes from south pole to north pole,
       ! reverse their order in "yf":
       yf(lats:1:-1) = yf(:)
       champd(:, lats:1:-1) = champd(:, :)
    ENDIF

    if (present(interbar)) then
       mid_values = interbar
    else
       mid_values = .true. ! default
    end if
    if (mid_values) then
       ! Replace longitudes and latitudes by their mid-values:
       old_xf_1 = xf(1)
       forall (i = 1: lons - 1) xf(i) = 0.5 * (xf(i) + xf(i+1))
       xf(lons) = 0.5 * (xf(lons) + old_xf_1 + 2 * pi)
          
       forall (j = 1: lats - 1) yf(j) = 0.5 * (yf(j) + yf(j+1))
    end if

  END SUBROUTINE conf_dat2d

end module conf_dat2d_m
