module getfieldindex_m

  implicit none

contains

  integer function GetFieldIndex(name)

    USE createnewfield_m, ONLY: fieldname, nbfield

    character(len=*), intent(in):: name

    ! Local:
    character(len(name)) TrueName

    !--------------------------------------------------

    TrueName = ADJUSTL(name)

    if (NbField >= 1) then
       GetFieldIndex = 1

       do while (TrueName /= FieldName(getfieldindex) &
            .and. GetFieldIndex < NbField)
          GetFieldIndex = GetFieldIndex + 1
       end do

       if (TrueName /= FieldName(getfieldindex)) GetFieldIndex = - 1
    else
       GetFieldIndex = - 1
    end if

  end function GetFieldIndex

end module getfieldindex_m
