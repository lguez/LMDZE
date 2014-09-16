module write_field3d_m

  implicit none

contains


  subroutine write_field3D(name,Field)

    use int2str_m, only: int2str

    integer, parameter :: MaxDim=3
    character(len=*)   :: name
    real, dimension(:,:,:) :: Field
    integer, dimension(MaxDim) :: Dim
    integer :: i,j,k,nb
    integer, parameter :: id=10
    integer, parameter :: NbCol=4
    integer :: ColumnSize
    integer :: pos,offset
    character(len=255) :: form
    character(len=255) :: spacing

    open(unit=id,file=name//'.field',form='formatted',status='replace')
    write (id,'("----- Field '//name//'"//)')

    Dim=shape(Field)
    offset=len(trim(int2str(Dim(1))))+len(trim(int2str(Dim(2))))+len(trim(int2str(Dim(3))))+4
    ColumnSize=22+6+3+offset

!    open(unit=id,file=name,form=formatted

    spacing='(t2,"'//repeat('-',ColumnSize*NbCol)//'")'

    do i=1,Dim(3)

      do j=1,Dim(2)
        nb=0
        Pos=2

        do k=1,Dim(1)
        nb=nb+1

          if (MOD(nb,NbCol)==0) then
            form='(t'//trim(int2str(pos))//            &
                 ',"('//trim(int2str(k))//','          &
                      //trim(int2str(j))//','          &
                      //trim(int2str(i))//')",t'       &
                      //trim(int2str(pos+offset))      &
                      //'," ---> ",g22.16,/)'
           Pos=2
          else
            form='(t'//trim(int2str(pos))//            &
                 ',"('//trim(int2str(k))//','          &
                      //trim(int2str(j))//','          &
                      //trim(int2str(i))//')",t'       &
                      //trim(int2str(pos+offset))      &
                      //'," ---> ",g22.16," | ")'
! d�pent de l'impl�mention, sur compaq, c'est necessaire
!            Pos=Pos+ColumnSize
          endif
          write (id,form,advance='no') Field(k,j,i)
        enddo
        if (MOD(nb,NbCol)==0) then
          write (id,spacing)
        else
          write (id,'("")')
          write (id,spacing)
        endif
      enddo
      write (id,spacing)
    enddo

    close(id)

  end subroutine write_field3D

end module write_field3d_m
