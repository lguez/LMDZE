subroutine correctbid(iim,nl,x)

  ! Warning Correction bidon pour palier a un probleme dans la
  ! creation des fichiers nc

  integer iim,nl
  real x(iim+1,nl)
  integer i,l
  real zz

  do l=1,nl
     do i=2,iim-1
        if(abs(x(i,l)).gt.1.e10) then
           zz=0.5*(x(i-1,l)+x(i+1,l))
           !              print*,'correction ',i,l,x(i,l),zz
           x(i,l)=zz
        endif
     enddo
  enddo

end subroutine correctbid
