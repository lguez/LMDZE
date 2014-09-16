module write_field

  ! From write_field.F90 1907 2013-11-26 13:10:46Z

  implicit none

  integer, parameter:: MaxWriteField = 100
  integer, dimension(MaxWriteField), save:: FieldId, FieldVarId, FieldIndex
  character(len=255), save:: FieldName(MaxWriteField)
  integer:: NbField = 0

end module write_field
