!
      SUBROUTINE cv_closure(nloc,ncum,nd,nk,icb &
                           ,tv,tvp,p,ph,dph,plcl,cpn &
                           ,iflag,cbmf)
            use cvthermo
            use cvparam
      implicit none

! inputs:
      integer ncum, nd, nloc
      integer nk(nloc), icb(nloc)
      real tv(nloc,nd), tvp(nloc,nd), p(nloc,nd), dph(nloc,nd)
      real ph(nloc,nd+1) ! caution nd instead ndp1 to be consistent...
      real plcl(nloc), cpn(nloc,nd)

! outputs:
      integer iflag(nloc)
      real cbmf(nloc) ! also an input

! local variables:
      integer i, k, icbmax
      real dtpbl(nloc), dtmin(nloc), tvpplcl(nloc), tvaplcl(nloc)
      real work(nloc)


!-------------------------------------------------------------------
! Compute icbmax.
!-------------------------------------------------------------------

      icbmax=2
      do 230 i=1,ncum
       icbmax=max(icbmax,icb(i))
 230  continue

!=====================================================================
! ---  CALCULATE CLOUD BASE MASS FLUX
!=====================================================================
!
! tvpplcl = parcel temperature lifted adiabatically from level
!           icb-1 to the LCL.
! tvaplcl = virtual temperature at the LCL.
!
      do 610 i=1,ncum
        dtpbl(i)=0.0
        tvpplcl(i)=tvp(i,icb(i)-1) &
        -rrd*tvp(i,icb(i)-1)*(p(i,icb(i)-1)-plcl(i)) &
        /(cpn(i,icb(i)-1)*p(i,icb(i)-1))
        tvaplcl(i)=tv(i,icb(i)) &
        +(tvp(i,icb(i))-tvp(i,icb(i)+1))*(plcl(i)-p(i,icb(i))) &
        /(p(i,icb(i))-p(i,icb(i)+1))
 610  continue

!-------------------------------------------------------------------
! --- Interpolate difference between lifted parcel and
! --- environmental temperatures to lifted condensation level
!-------------------------------------------------------------------
!
! dtpbl = average of tvp-tv in the PBL (k=nk to icb-1).
!
      do 630 k=minorig,icbmax
        do 620 i=1,ncum
        if((k.ge.nk(i)).and.(k.le.(icb(i)-1)))then
          dtpbl(i)=dtpbl(i)+(tvp(i,k)-tv(i,k))*dph(i,k)
        endif
 620    continue
 630  continue
      do 640 i=1,ncum
        dtpbl(i)=dtpbl(i)/(ph(i,nk(i))-ph(i,icb(i)))
        dtmin(i)=tvpplcl(i)-tvaplcl(i)+dtmax+dtpbl(i)
 640  continue
!
!-------------------------------------------------------------------
! --- Adjust cloud base mass flux
!-------------------------------------------------------------------
!
      do 650 i=1,ncum
       work(i)=cbmf(i)
       cbmf(i)=max(0.0,(1.0-damp)*cbmf(i)+0.1*alpha*dtmin(i))
       if((work(i).eq.0.0).and.(cbmf(i).eq.0.0))then
         iflag(i)=3
       endif
 650  continue

       return
       end
