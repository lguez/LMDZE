module createnewfield_m

  implicit none

contains

  subroutine CreateNewField(name, dimx, dimy, dimz)

    USE netcdf, ONLY: nf90_clobber, nf90_double, nf90_unlimited
    USE netcdf95, ONLY: nf95_create, nf95_def_dim, nf95_def_var, nf95_enddef
    USE write_field, ONLY: fieldid, fieldindex, fieldname, fieldvarid, nbfield

    character(len = *), intent(in):: name
    integer, intent(in):: dimx, dimy, dimz

    ! Local:
    integer TabDim(4)

    !--------------------------------------------------------------------

    NbField = NbField + 1
    FieldName(NbField) = TRIM(ADJUSTL(name))
    FieldIndex(NbField) = 1

    call NF95_CREATE(TRIM(ADJUSTL(name)) // '.nc', NF90_CLOBBER, &
         FieldId(NbField))
    call NF95_DEF_DIM(FieldId(NbField), 'X', dimx, TabDim(1))
    call NF95_DEF_DIM(FieldId(NbField), 'Y', dimy, TabDim(2))
    call NF95_DEF_DIM(FieldId(NbField), 'Z', dimz, TabDim(3))
    call NF95_DEF_DIM(FieldId(NbField), 'iter', NF90_UNLIMITED, TabDim(4))
    call NF95_DEF_VAR(FieldId(NbField), FieldName(NbField), NF90_DOUBLE, &
         TabDim, FieldVarId(NbField))
    call NF95_ENDDEF(FieldId(NbField))

  end subroutine CreateNewField

end module createnewfield_m
