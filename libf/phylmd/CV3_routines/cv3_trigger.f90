
      SUBROUTINE cv3_trigger(len,nd,icb,plcl,p,th,tv,tvp &
                      ,pbase,buoybase,iflag,sig,w0)
            use cv3_param_m
      implicit none

!-------------------------------------------------------------------
! --- TRIGGERING
!
!     - computes the cloud base
!   - triggering (crude in this version)
!     - relaxation of sig and w0 when no convection
!
!     Caution1: if no convection, we set iflag=4
!              (it used to be 0 in convect3)
!
!     Caution2: at this stage, tvp (and thus buoy) are know up
!             through icb only!
! -> the buoyancy below cloud base not (yet) set to the cloud base buoyancy
!-------------------------------------------------------------------


! input:
      integer len, nd
      integer icb(len)
      real plcl(len), p(len,nd)
      real th(len,nd), tv(len,nd), tvp(len,nd)

! output:
      real pbase(len), buoybase(len)

! input AND output:
      integer iflag(len)
      real, intent(inout):: sig(len,nd), w0(len,nd)

! local variables:
      integer i,k
      real tvpbase, tvbase, tdif, ath, ath1

!
! ***   set cloud base buoyancy at (plcl+dpbase) level buoyancy
!
      do 100 i=1,len
       pbase(i) = plcl(i) + dpbase
       tvpbase = tvp(i,icb(i))*(pbase(i)-p(i,icb(i)+1)) &
                              /(p(i,icb(i))-p(i,icb(i)+1)) &
               + tvp(i,icb(i)+1)*(p(i,icb(i))-pbase(i)) &
                                /(p(i,icb(i))-p(i,icb(i)+1))
       tvbase = tv(i,icb(i))*(pbase(i)-p(i,icb(i)+1)) &
                            /(p(i,icb(i))-p(i,icb(i)+1)) &
              + tv(i,icb(i)+1)*(p(i,icb(i))-pbase(i)) &
                              /(p(i,icb(i))-p(i,icb(i)+1))
       buoybase(i) = tvpbase - tvbase
100   continue

!
!   ***   make sure that column is dry adiabatic between the surface  ***
!   ***    and cloud base, and that lifted air is positively buoyant  ***
!   ***                         at cloud base                         ***
!   ***       if not, return to calling program after resetting       ***
!   ***                        sig(i) and w0(i)                       ***
!

      do  60 k=1,nl
      do 200 i=1,len

       tdif = buoybase(i)
       ath1 = th(i,1)
       ath  = th(i,icb(i)-1) - dttrig

       if (tdif.lt.dtcrit .or. ath.gt.ath1) then
            sig(i,k) = beta*sig(i,k) - 2.*alpha*tdif*tdif
            sig(i,k) = AMAX1(sig(i,k),0.0)
            w0(i,k)  = beta*w0(i,k)
        iflag(i)=4 ! pour version vectorisee
! convect3         iflag(i)=0
       endif

200   continue
 60   continue

! fin oct3 --

      return
      end
