
      SUBROUTINE cv3_feed(len,nd,t,q,qs,p,ph,hm,gz &
                        ,nk,icb,icbmax,iflag,tnk,qnk,gznk,plcl)
            use cvparam3
      implicit none

!================================================================
! Purpose: CONVECTIVE FEED
!
! Main differences with cv_feed:
!   - ph added in input
!     - here, nk(i)=minorig
!     - icb defined differently (plcl compared with ph instead of p)
!
! Main differences with convect3:
!     - we do not compute dplcldt and dplcldr of CLIFT anymore
!     - values iflag different (but tests identical)
!   - A,B explicitely defined (!...)
!================================================================


! inputs:
        integer len, nd
      real, intent(in):: t(len,nd)
      real q(len,nd), qs(len,nd), p(len,nd)
      real hm(len,nd), gz(len,nd)
      real ph(len,nd+1)

! outputs:
        integer iflag(len), nk(len), icb(len), icbmax
      real tnk(len), qnk(len), gznk(len), plcl(len)

! local variables:
      integer i, k
      integer ihmin(len)
      real work(len)
      real pnk(len), qsnk(len), rh(len), chi(len)
      real A, B ! convect3
!ym
      plcl=0.0
!@ !-------------------------------------------------------------------
!@ ! --- Find level of minimum moist static energy
!@ ! --- If level of minimum moist static energy coincides with
!@ ! --- or is lower than minimum allowable parcel origin level,
!@ ! --- set iflag to 6.
!@ !-------------------------------------------------------------------
!@
!@       do 180 i=1,len
!@        work(i)=1.0e12
!@        ihmin(i)=nl
!@  180  continue
!@       do 200 k=2,nlp
!@         do 190 i=1,len
!@          if((hm(i,k).lt.work(i)).and.
!@      &      (hm(i,k).lt.hm(i,k-1)))then
!@            work(i)=hm(i,k)
!@            ihmin(i)=k
!@          endif
!@  190    continue
!@  200  continue
!@       do 210 i=1,len
!@         ihmin(i)=min(ihmin(i),nlm)
!@         if(ihmin(i).le.minorig)then
!@           iflag(i)=6
!@         endif
!@  210  continue
!@ c
!@ !-------------------------------------------------------------------
!@ ! --- Find that model level below the level of minimum moist static
!@ ! --- energy that has the maximum value of moist static energy
!@ !-------------------------------------------------------------------
!@
!@       do 220 i=1,len
!@        work(i)=hm(i,minorig)
!@        nk(i)=minorig
!@  220  continue
!@       do 240 k=minorig+1,nl
!@         do 230 i=1,len
!@          if((hm(i,k).gt.work(i)).and.(k.le.ihmin(i)))then
!@            work(i)=hm(i,k)
!@            nk(i)=k
!@          endif
!@  230     continue
!@  240  continue

!-------------------------------------------------------------------
! --- Origin level of ascending parcels for convect3:
!-------------------------------------------------------------------

         do 220 i=1,len
          nk(i)=minorig
  220    continue

!-------------------------------------------------------------------
! --- Check whether parcel level temperature and specific humidity
! --- are reasonable
!-------------------------------------------------------------------
       do 250 i=1,len
       if( (     ( t(i,nk(i)).lt.250.0    ) &
             .or.( q(i,nk(i)).le.0.0      )     ) &
         .and. &
             ( iflag(i).eq.0) ) iflag(i)=7
 250   continue
!-------------------------------------------------------------------
! --- Calculate lifted condensation level of air at parcel origin level
! --- (Within 0.2% of formula of Bolton, MON. WEA. REV.,1980)
!-------------------------------------------------------------------

       A = 1669.0 ! convect3
       B = 122.0  ! convect3

       do 260 i=1,len

        if (iflag(i).ne.7) then ! modif sb Jun7th 2002

        tnk(i)=t(i,nk(i))
        qnk(i)=q(i,nk(i))
        gznk(i)=gz(i,nk(i))
        pnk(i)=p(i,nk(i))
        qsnk(i)=qs(i,nk(i))
!
        rh(i)=qnk(i)/qsnk(i)
! ori        rh(i)=min(1.0,rh(i)) ! removed for convect3
! ori        chi(i)=tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
        chi(i)=tnk(i)/(A-B*rh(i)-tnk(i)) ! convect3
        plcl(i)=pnk(i)*(rh(i)**chi(i))
        if(((plcl(i).lt.200.0).or.(plcl(i).ge.2000.0)) &
         .and.(iflag(i).eq.0))iflag(i)=8

        endif ! iflag=7

 260   continue

!-------------------------------------------------------------------
! --- Calculate first level above lcl (=icb)
!-------------------------------------------------------------------

!@      do 270 i=1,len
!@       icb(i)=nlm
!@ 270  continue
!@c
!@      do 290 k=minorig,nl
!@        do 280 i=1,len
!@          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i)))
!@     &    icb(i)=min(icb(i),k)
!@ 280    continue
!@ 290  continue
!@c
!@      do 300 i=1,len
!@        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
!@ 300  continue

      do 270 i=1,len
       icb(i)=nlm
 270  continue
!
! la modification consiste a comparer plcl a ph et non a p:
! icb est defini par :  ph(icb)<plcl<ph(icb-1)
!@      do 290 k=minorig,nl
      do 290 k=3,nl-1 ! modif pour que icb soit sup/egal a 2
        do 280 i=1,len
          if( ph(i,k).lt.plcl(i) ) icb(i)=min(icb(i),k)
 280    continue
 290  continue
!
      do 300 i=1,len
!@        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
        if((icb(i).eq.nlm).and.(iflag(i).eq.0))iflag(i)=9
 300  continue

      do 400 i=1,len
        icb(i) = icb(i)-1 ! icb sup ou egal a 2
 400  continue
!
! Compute icbmax.
!
      icbmax=2
      do 310 i=1,len
!!        icbmax=max(icbmax,icb(i))
       if (iflag(i).lt.7) icbmax=max(icbmax,icb(i)) ! sb Jun7th02
 310  continue

      return
      end
