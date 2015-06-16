module WriteField_m

  use CreateNewField_m, only: CreateNewField, ncid, nbfield, MaxWriteField
  use GetFieldIndex_m, only: GetFieldIndex
  USE netcdf95, ONLY: nf95_put_var, nf95_sync

  implicit none

  integer, save:: Record(MaxWriteField)

  interface WriteField
     module procedure WriteField3d,WriteField2d,WriteField1d
  end interface WriteField

  private
  public WriteField

contains

  subroutine WriteField1d(name,Field)

    character(len=*), intent(in):: name
    real, intent(in):: Field(:)

    ! Local:
    integer index

    !-------------------------------------------

    Index=GetFieldIndex(name)

    if (Index==-1) then
       call CreateNewField(name, shape(field))
       Index=nbfield
       Record(Index) = 1
    else
       Record(Index)=Record(Index)+1
    endif

    call NF95_PUT_VAR(Ncid(Index), Varid = 1, values = Field, &
         start = (/1, Record(Index)/))
    call nf95_sync(Ncid(Index))

  end subroutine WriteField1d

  !****************************************************************

  subroutine WriteField2d(name,Field)

    character(len=*), intent(in):: name
    real, intent(in):: Field(:, :)

    ! Local:
    integer index

    !-------------------------------------------

    Index=GetFieldIndex(name)

    if (Index==-1) then
       call CreateNewField(name, shape(field))
       Index=nbfield
       Record(Index) = 1
    else
       Record(Index)=Record(Index)+1
    endif

    call NF95_PUT_VAR(Ncid(Index), Varid = 1, values = Field, &
         start = (/1, 1, Record(Index)/))
    call nf95_sync(Ncid(Index))

  end subroutine WriteField2d

  !****************************************************************

  subroutine WriteField3d(name,Field)

    character(len=*), intent(in):: name
    real, intent(in):: Field(:, :, :)

    ! Local:
    integer index

    !-------------------------------------------

    Index=GetFieldIndex(name)

    if (Index==-1) then
       call CreateNewField(name, shape(field))
       Index=nbfield
       Record(Index) = 1
    else
       Record(Index)=Record(Index)+1
    endif

    call NF95_PUT_VAR(Ncid(Index), Varid = 1, values = Field, &
         start = (/1, 1, 1, Record(Index)/))
    call nf95_sync(Ncid(Index))

  end subroutine WriteField3d

end module WriteField_m
