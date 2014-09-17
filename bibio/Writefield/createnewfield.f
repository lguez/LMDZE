module createnewfield_m

  implicit none

contains

  subroutine CreateNewField(name, my_shape)

    USE netcdf, ONLY: nf90_clobber, nf90_double, nf90_unlimited
    USE netcdf95, ONLY: nf95_create, nf95_def_dim, nf95_def_var, nf95_enddef
    USE write_field, ONLY: ncid, record, fieldname, varid, nbfield

    character(len = *), intent(in):: name
    integer, intent(in):: my_shape(:)

    ! Local:
    integer Dimid(size(my_shape) + 1), n_dim

    !--------------------------------------------------------------------

    NbField = NbField + 1
    FieldName(NbField) = ADJUSTL(name)
    Record(NbField) = 1
    n_dim = size(my_shape)

    call NF95_CREATE(TRIM(FieldName(NbField)) // '.nc', NF90_CLOBBER, &
         Ncid(NbField))
    call NF95_DEF_DIM(Ncid(NbField), 'X', my_shape(1), Dimid(1))
    if (n_dim >= 2) then
       call NF95_DEF_DIM(Ncid(NbField), 'Y', my_shape(2), Dimid(2))
       if (n_dim >= 3) then
          call NF95_DEF_DIM(Ncid(NbField), 'Z', my_shape(3), Dimid(3))
       end if
    end if
    call NF95_DEF_DIM(Ncid(NbField), 'iter', NF90_UNLIMITED, Dimid(n_dim + 1))
    call NF95_DEF_VAR(Ncid(NbField), FieldName(NbField), NF90_DOUBLE, &
         Dimid, Varid(NbField))
    call NF95_ENDDEF(Ncid(NbField))

  end subroutine CreateNewField

end module createnewfield_m
