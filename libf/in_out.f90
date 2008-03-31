module new_unit_m

  implicit none

contains

  subroutine new_unit(unit)

    integer, intent(out):: unit

    ! Variables local to the procedure:
    logical opened, exist

    !------------------------------------------------------

    unit = 0
    do
       inquire(unit=unit, opened=opened, exist=exist)
       if (exist .and. .not. opened) exit
       unit = unit + 1
    end do

  end subroutine new_unit

end module new_unit_m
