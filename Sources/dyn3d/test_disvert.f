module test_disvert_m

  implicit none

contains

  subroutine test_disvert

    ! Author: Lionel GUEZ

    ! This procedure tests the order of pressure values at half-levels
    ! and full levels, and writes the vertical coordinates to files
    ! "half_level.csv" and "full_level.csv" for a reference surface
    ! pressure. We arbitrarily choose to test ngrid values of the
    ! surface pressure, which sample possible values on Earth.

    use abort_gcm_m, only: abort_gcm
    use comconst, only: kappa, cpp
    use disvert_m, only: ap, bp, preff, presnivs, s
    use dimens_m, only: llm
    use exner_hyb_m, only: exner_hyb
    use jumble, only: new_unit

    integer i, unit, l

    integer, parameter:: ngrid = 8

    real p(ngrid, 1, llm + 1) ! pressure at half-level, in Pa
    real z(llm) ! pressure-altitude at half-level (km)
    real pks(ngrid, 1) ! exner function at the surface, in J K-1 kg-1 
    real pk(ngrid, 1, llm) ! exner function at full level, in J K-1 kg-1 
    real ps(ngrid, 1) ! surface pressure, in Pa
    real p_lay(ngrid, llm) ! pressure at full level, in Pa
    real, parameter:: delta_ps = 6e4 / (ngrid - 2) ! in Pa

    !---------------------

    print *, "Call sequence information: test_disvert"

    ps(:, 1) = (/(5e4 + delta_ps * i, i = 0, ngrid - 2), preff/)
    forall (l = 1: llm + 1) p(:, 1, l) = ap(l) + bp(l) * ps(:, 1)
    call exner_hyb(ps, p, pks, pk)
    p_lay = preff * (pk(:, 1, :) / cpp)**(1. / kappa)

    ! Write distribution for the reference surface pressure (index ngrid):

    z = 7. * log(preff / p(ngrid, 1, :llm))
    call new_unit(unit)

    open(unit, file="half_level.csv", status="replace", action="write")
    ! Title line:
    write(unit, fmt=*) '"ap (hPa)" "bp" "s" "pressure (hPa)" "z (km)"'
    do l = 1, llm
       write(unit, fmt=*) ap(l) / 100., bp(l), s(l), p(ngrid, 1, l) / 100., z(l)
    end do
    close(unit)
    print *, 'The file "half_level.csv" has been created.'

    open(unit, file="full_level.csv", status="replace", action="write")
    ! Title line:
    write(unit, fmt=*) &
         '"pressure (hPa)" "z (km)" "presnivs (hPa)" "delta z (km)"'
    do l = 1, llm - 1
       write(unit, fmt=*) p_lay(ngrid, l) / 100., &
            7. * log(preff / p_lay(ngrid, l)), presnivs(l) / 100., z(l+1) - z(l)
    end do
    write(unit, fmt=*) p_lay(ngrid, llm) / 100., &
         7. * log(preff / p_lay(ngrid, llm)), presnivs(llm) / 100.
    close(unit)
    print *, 'The file "full_level.csv" has been created.'

    ! Are pressure values in the right order?
    if (any(p(:, 1, :llm) <= p_lay .or. p_lay <= p(:, 1, 2:))) then
       ! List details and stop:
       do l = 1, llm
          do i = 1, ngrid
             if (p(i, 1, l) <= p_lay(i, l)) then
                print 1000, "ps = ", ps(i, 1) / 100., "hPa, p(level ",  l, &
                     ") = ", p(i, 1, l) / 100., " hPa <= p(layer ", l, ") = ", &
                     p_lay(i, l) / 100., " hPa"
             end if
             if (p_lay(i, l) <= p(i, 1, l+1)) then
                print 1000, "ps = ", ps(i, 1) / 100., &
                     "hPa, p(layer ", l, ") = ", p_lay(i, l) / 100., &
                     " hPa <= p(level ", l + 1, ") = ", &
                     p(i, 1, l + 1) / 100., " hPa"
             end if
          end do
       end do
       call abort_gcm("test_disvert", "bad order of pressure values")
    end if

1000 format (3(a, g10.4, a, i0))

  end subroutine test_disvert

end module test_disvert_m
