module getfieldindex_m

  implicit none

contains

  integer function GetFieldIndex(name)

    USE write_field, ONLY: fieldname, nbfield

    character(len=*), intent(in):: name

    ! Local:
    character(len=255) TrueName
    integer i

    !--------------------------------------------------

    TrueName=TRIM(ADJUSTL(name))

    GetFieldIndex=-1

    do i=1,NbField
       if (TrueName==FieldName(i)) then
          GetFieldIndex=i
          exit
       endif
    enddo

  end function GetFieldIndex

end module getfieldindex_m
