
      SUBROUTINE cv_trigger(len,nd,icb,cbmf,tv,tvp,iflag)
            use cvparam
      implicit none

!-------------------------------------------------------------------
! --- Test for instability.
! --- If there was no convection at last time step and parcel
! --- is stable at icb, then set iflag to 4.
!-------------------------------------------------------------------


! inputs:
       integer, intent(in):: len, nd, icb(len)
       real cbmf(len), tv(len,nd), tvp(len,nd)

! outputs:
       integer iflag(len) ! also an input

! local variables:
       integer i


      do 390 i=1,len
        if((cbmf(i).eq.0.0) .and.(iflag(i).eq.0).and. &
        (tvp(i,icb(i)).le.(tv(i,icb(i))-dtmax)))iflag(i)=4
 390  continue

      return
      end
