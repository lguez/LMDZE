module inter_barxy_m

  ! From inter_barxy.F, version 1.1.1.1, 2004/05/19 12:53:07

  implicit none

contains

  SUBROUTINE inter_barxy(dlonid, dlatid, champ, rlonimod, rlatimod, champint)

    ! Author: P. Le Van

    use jumble, only: assert_eq, assert

    use comgeom, only: aire_2d, apoln, apols
    use dimensions, only: iim, jjm
    use inter_barx_m, only: inter_barx
    use inter_bary_m, only: inter_bary
    use ord_coord_m, only: ord_coord
    use ord_coordm_m, only: ord_coordm

    REAL, intent(in):: dlonid(:)
    ! longitude from input file, in rad, from -pi to pi

    REAL, intent(in):: dlatid(:), champ(:, :), rlonimod(:)

    REAL, intent(in):: rlatimod(:)
    ! latitude angle, in degrees or rad, in strictly decreasing order

    real, intent(out):: champint(:, :)
    ! Si taille de la seconde dim = jjm + 1, on veut interpoler sur les
    ! jjm+1 latitudes rlatu du modele (latitudes des scalaires et de U)
    ! Si taille de la seconde dim = jjm, on veut interpoler sur les
    ! jjm latitudes rlatv du mod\`ele (latitudes de V)

    ! Local:

    REAL champy(iim, size(champ, 2))
    integer j, i, jnterfd, jmods

    REAL yjmod(size(champint, 2))
    ! (angle, in degrees, in strictly increasing order)

    REAL   yjdat(size(dlatid) + 1) ! angle, in degrees, in increasing order
    LOGICAL decrois ! "dlatid" is in decreasing order

    !-----------------------------------

    jnterfd = assert_eq(size(champ, 2) - 1, size(dlatid), "inter_barxy jnterfd")
    jmods = size(champint, 2)
    call assert(size(champ, 1) == size(dlonid), "inter_barxy size(champ, 1)")
    call assert((/size(rlonimod), size(champint, 1)/) == iim, &
         "inter_barxy iim")
    call assert(any(jmods == (/jjm, jjm + 1/)), 'inter_barxy jmods')
    call assert(size(rlatimod) == jjm, "inter_barxy size(rlatimod)")

    ! Check decreasing order for "rlatimod":
    DO i = 2, jjm
       IF (rlatimod(i) >= rlatimod(i-1)) then
          print *, '"inter_barxy": "rlatimod" should be strictly decreasing'
          stop 1
       end IF
    ENDDO

    yjmod(:jjm) = ord_coordm(rlatimod)
    IF (jmods == jjm + 1) THEN
       IF (90. - yjmod(jjm) < 0.01) then
          print *, '"inter_barxy": with jmods = jjm + 1, ' &
               // 'yjmod(jjm) should be < 90.'
          stop 1
       end IF
    ELSE
       ! jmods = jjm
       IF (ABS(yjmod(jjm) - 90.) > 0.01) then
          print *, '"inter_barxy": with jmods = jjm, yjmod(jjm) should be 90.'
          stop 1
       end IF
    ENDIF

    if (jmods == jjm + 1) yjmod(jjm + 1) = 90.

    DO j = 1, jnterfd + 1
       champy(:, j) = inter_barx(dlonid, champ(:, j), rlonimod)
    ENDDO

    CALL ord_coord(dlatid, yjdat, decrois)
    IF (decrois) champy(:, :) = champy(:, jnterfd + 1:1:-1)
    DO i = 1, iim
       champint(i, :) = inter_bary(yjdat, champy(i, :), yjmod)
    ENDDO
    champint(:, :) = champint(:, jmods:1:-1)

    IF (jmods == jjm + 1) THEN
       ! Valeurs uniques aux poles
       champint(:, 1) = SUM(aire_2d(:iim,  1) * champint(:, 1)) / apoln
       champint(:, jjm + 1) = SUM(aire_2d(:iim, jjm + 1) &
            * champint(:, jjm + 1)) / apols
    ENDIF

  END SUBROUTINE inter_barxy

end module inter_barxy_m
