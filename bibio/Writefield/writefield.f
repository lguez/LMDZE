module WriteField_m

  use writefield_gen_m, only: writefield_gen

  implicit none

  interface WriteField
     module procedure WriteField3d,WriteField2d,WriteField1d
  end interface WriteField

  private
  public WriteField

contains

  subroutine WriteField1d(name,Field)
    character(len=*) :: name
    real, dimension(:) :: Field
    integer, dimension(1) :: Dim

    Dim=shape(Field)
    call WriteField_gen(name,Field,Dim(1),1,1)

  end subroutine WriteField1d

  !****************************************************************

  subroutine WriteField2d(name,Field)
    character(len=*) :: name
    real, dimension(:,:) :: Field
    integer, dimension(2) :: Dim

    Dim=shape(Field)
    call WriteField_gen(name,Field,Dim(1),Dim(2),1)

  end subroutine WriteField2d

  !****************************************************************

  subroutine WriteField3d(name,Field)
    character(len=*) :: name
    real, dimension(:,:,:) :: Field
    integer, dimension(3) :: Dim

    Dim=shape(Field)
    call WriteField_gen(name,Field,Dim(1),Dim(2),Dim(3))

  end subroutine WriteField3d

end module WriteField_m
