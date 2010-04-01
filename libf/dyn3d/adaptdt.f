!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/adaptdt.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      subroutine adaptdt(nadv,dtbon,n,pbaru,
     c                   masse)

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use conf_gcm_m
      use logic
      use comgeom
      use temps
      use ener
      IMPLICIT NONE


c----------------------------------------------------------
c     Arguments
c----------------------------------------------------------
      INTEGER n,nadv
      REAL dtbon 
      REAL, intent(in):: pbaru(iip1,jjp1,llm)
      REAL masse(iip1,jjp1,llm)
c----------------------------------------------------------    
c     Local
c----------------------------------------------------------
      INTEGER i,j,l
      REAL CFLmax,aaa,bbb
      
        CFLmax=0.
        do l=1,llm
         do j=2,jjm
          do i=1,iim
             aaa=pbaru(i,j,l)*dtvr/masse(i,j,l)
             CFLmax=max(CFLmax,aaa)
             bbb=-pbaru(i,j,l)*dtvr/masse(i+1,j,l)
             CFLmax=max(CFLmax,bbb)
          enddo
         enddo
        enddo              
        n=int(CFLmax)+1
c pour reproduire cas VL du code qui appele x,y,z,y,x
c        if (nadv.eq.30) n=n/2   ! Pour Prather
        dtbon=dtvr/n
        
       return
       end







