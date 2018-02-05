module ord_coord_m

  implicit none

contains


  !******************************

  SUBROUTINE ord_coord(xi, xo, decrois)

    ! From dyn3d/ord_coord.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Author : P. Le Van

    ! This procedure receives an array of latitudes.
    ! It converts them to degrees if they are in radians.
    ! If the input latitudes are in decreasing order, the procedure
    ! reverses their order.
    ! Finally, the procedure adds 90° as the last value of the array.

    use nr_util, only: assert_eq, pi


    REAL, intent(in):: xi(:)
    ! (latitude, in degrees or radians, in increasing or decreasing order)
    ! ("xi" should contain latitudes from pole to pole.
    ! "xi" should contain the latitudes of the boundaries of grid
    ! cells, not the centers of grid cells.
    ! So the extreme values should not be 90° and -90°.)

    REAL, intent(out):: xo(:) ! angles in degrees
    LOGICAL, intent(out):: decrois

    ! Variables  local to the procedure:
    INTEGER nmax, i

    !--------------------

    nmax = assert_eq(size(xi), size(xo) - 1, "ord_coord")

    ! Check monotonicity:
    decrois = xi(2) < xi(1)
    DO i = 3, nmax
       IF (decrois .neqv. xi(i) < xi(i-1)) then
          print *, '"ord_coord":  latitudes are not monotonic'
          stop 1
       end IF
    ENDDO

    IF (abs(xi(1)) < pi) then
       ! "xi" contains latitudes in radians
       xo(:nmax) = xi(:) * 180. / pi
    else
       ! "xi" contains latitudes in degrees
       xo(:nmax) = xi(:)
    end IF

    IF (ABS(abs(xo(1)) - 90) < 0.001 .or. ABS(abs(xo(nmax)) - 90) < 0.001) THEN
       print *, "ord_coord"
       PRINT *, '"xi" should contain the latitudes of the boundaries of ' &
            // 'grid cells, not the centers of grid cells.'
       STOP 1
    ENDIF

    IF (decrois) xo(:nmax) = xo(nmax:1:- 1)
    xo(nmax + 1) = 90.

  END SUBROUTINE ord_coord

end module ord_coord_m
