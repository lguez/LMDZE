module writefield_gen_m

  implicit none

contains

  subroutine WriteField_gen(name, Field, dimx, dimy, dimz)

    use CreateNewField_m, only: CreateNewField
    use GetFieldIndex_m, only: GetFieldIndex
    USE netcdf95, ONLY: nf95_put_var
    USE write_field, ONLY: fieldid, fieldindex, fieldvarid

    character(len=*) :: name
    integer :: dimx, dimy, dimz
    real, dimension(dimx, dimy, dimz) :: Field
    integer :: index
    integer :: start(4)
    integer :: count(4)

    !------------------------------------------------------------

    Index=GetFieldIndex(name)
    if (Index==-1) then
       call CreateNewField(name, dimx, dimy, dimz)
       Index=GetFieldIndex(name)
    else
       FieldIndex(Index)=FieldIndex(Index)+1
    endif

    start(1)=1
    start(2)=1
    start(3)=1
    start(4)=FieldIndex(Index)

    count(1)=dimx
    count(2)=dimy
    count(3)=dimz
    count(4)=1

    call NF95_PUT_VAR(FieldId(Index), FieldVarId(Index), Field, start, count)

  end subroutine WriteField_gen

end module writefield_gen_m
