module write_field2d_m

  implicit none

contains


  subroutine write_field2D(name,Field)

    use int2str_m, only: int2str

    integer, parameter :: MaxDim=2
    character(len=*)   :: name
    real, dimension(:,:) :: Field
    integer, dimension(MaxDim) :: Dim
    integer :: i,j,nb
    integer, parameter :: id=10
    integer, parameter :: NbCol=4
    integer :: ColumnSize
    integer :: pos,offset
    character(len=255) :: form
    character(len=255) :: spacing

    open(unit=id,file=name//'.field',form='formatted',status='replace')
    write (id,'("----- Field '//name//'",//)')

    Dim=shape(Field)
    offset=len(trim(int2str(Dim(1))))+len(trim(int2str(Dim(2))))+3
    ColumnSize=20+6+3+offset

    spacing='(t2,"'//repeat('-',ColumnSize*NbCol)//'")'

    do i=1,Dim(2)
      nb=0
      Pos=2
      do j=1,Dim(1)
        nb=nb+1

        if (MOD(nb,NbCol)==0) then
          form='(t'//trim(int2str(pos))//            &
               ',"('//trim(int2str(j))//','          &
                    //trim(int2str(i))//')",t'       &
                    //trim(int2str(pos+offset))     &
                    //'," ---> ",g22.16,/)'
          Pos=2
        else
          form='(t'//trim(int2str(pos))//            &
               ',"('//trim(int2str(j))//','          &
                    //trim(int2str(i))//')",t'       &
                    //trim(int2str(pos+offset))     &
                    //'," ---> ",g22.16," | ")'
          Pos=Pos+ColumnSize
        endif
        write (id,form,advance='no') Field(j,i)
      enddo
      if (MOD(nb,NbCol)==0) then
        write (id,spacing)
      else
        write (id,'("")')
        write (id,spacing)
      endif
    enddo

  end subroutine write_field2D

end module write_field2d_m
