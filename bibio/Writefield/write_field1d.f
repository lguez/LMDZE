module write_field1d_m

  implicit none

contains

  subroutine write_field1D(name,Field)

    use int2str_m, only: int2str

    integer, parameter :: MaxDim=1
    character(len=*)   :: name
    real, dimension(:) :: Field
    integer, dimension(MaxDim) :: Dim
    integer :: i,nb
    integer, parameter :: id=10
    integer, parameter :: NbCol=4
    integer :: ColumnSize
    integer :: pos
    character(len=255) :: form
    character(len=255) :: MaxLen


    open(unit=id,file=name//'.field',form='formatted',status='replace')
    write (id,'("----- Field '//name//'",//)')
    Dim=shape(Field)
    MaxLen=int2str(len(trim(int2str(Dim(1)))))
    ColumnSize=20+6+3+len(trim(int2str(Dim(1))))
    Nb=0
    Pos=2
    do i=1,Dim(1)
      nb=nb+1

      if (MOD(nb,NbCol)==0) then
        form='(t'//trim(int2str(pos))// ',i'//trim(MaxLen) //'," ---> ",g22.16,/)'
        Pos=2
      else
        form='(t'//trim(int2str(pos))// ',i'//trim(MaxLen) //'," ---> ",g22.16," | ",)'
        Pos=Pos+ColumnSize
      endif
      write (id,form,advance='no') i,Field(i)
    enddo

    close(id)

  end subroutine write_field1D

end module write_field1d_m
