module int2str_m

  implicit none

contains

  function int2str(i)

    integer, parameter :: MaxLen=10
    integer,intent(in) :: i
    character(len=MaxLen) :: int2str

    !---------------------------------------------------

    write(unit = int2str, fmt = *) i
    int2str = adjustl(int2str)

  end function int2str

end module int2str_m
