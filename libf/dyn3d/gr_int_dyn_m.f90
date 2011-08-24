module gr_int_dyn_m

  ! Clean: no C preprocessor directive, no include line

  implicit none

contains

  function gr_int_dyn(champin)

    ! From dyn3d/gr_int_dyn.F,v 1.1.1.1 2004/05/19 12:53:07

    ! Passage d'un champ interpolé à un champ sur grille scalaire

    REAL, intent(in):: champin(:, :)
    REAL gr_int_dyn(size(champin, 1) + 1, size(champin, 2))

    ! Variables local to the procedure:
    integer iim, jp1

    !-----------------------------------------------------------------------

    iim = size(champin, 1)
    jp1 = size(champin, 2)

    gr_int_dyn(:, 1) = sum(champin(:, 1)) / iim ! north pole
    gr_int_dyn(:, jp1) = sum(champin(:, jp1)) / iim ! south pole
    gr_int_dyn(: iim, 2: jp1 - 1) = champin(:, 2: jp1 - 1)
    gr_int_dyn(iim + 1, 2: jp1 - 1) = gr_int_dyn(1, 2: jp1 - 1)

  END function gr_int_dyn

end module gr_int_dyn_m
