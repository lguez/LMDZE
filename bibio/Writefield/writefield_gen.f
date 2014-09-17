module writefield_gen_m

  implicit none

contains

  subroutine WriteField_gen(name, Field, dimx, dimy, dimz)

    use CreateNewField_m, only: CreateNewField
    use GetFieldIndex_m, only: GetFieldIndex
    USE netcdf95, ONLY: nf95_put_var
    USE write_field, ONLY: ncid, record, varid, nbfield

    character(len=*) :: name
    integer :: dimx, dimy, dimz
    real, dimension(dimx, dimy, dimz) :: Field
    integer :: index

    !------------------------------------------------------------

    Index=GetFieldIndex(name)
    if (Index==-1) then
       call CreateNewField(name, (/dimx, dimy, dimz/))
       Index=nbfield
    else
       Record(Index)=Record(Index)+1
    endif

    call NF95_PUT_VAR(Ncid(Index), Varid(Index), Field, &
         start = (/1, 1, 1, Record(Index)/))

  end subroutine WriteField_gen

end module writefield_gen_m
