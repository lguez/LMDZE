
      SUBROUTINE cv_feed(len,nd,t,q,qs,p,hm,gz &
                        ,nk,icb,icbmax,iflag,tnk,qnk,gznk,plcl)
            use cvparam
      implicit none

!================================================================
! Purpose: CONVECTIVE FEED
!================================================================


! inputs:
        integer len, nd
      real, intent(in):: t(len,nd)
      real q(len,nd), qs(len,nd), p(len,nd)
      real hm(len,nd), gz(len,nd)

! outputs:
        integer iflag(len), nk(len), icb(len), icbmax
      real tnk(len), qnk(len), gznk(len), plcl(len)

! local variables:
      integer i, k
      integer ihmin(len)
      real work(len)
      real pnk(len), qsnk(len), rh(len), chi(len)

!-------------------------------------------------------------------
! --- Find level of minimum moist static energy
! --- If level of minimum moist static energy coincides with
! --- or is lower than minimum allowable parcel origin level,
! --- set iflag to 6.
!-------------------------------------------------------------------

      do 180 i=1,len
       work(i)=1.0e12
       ihmin(i)=nl
 180  continue
      do 200 k=2,nlp
        do 190 i=1,len
         if((hm(i,k).lt.work(i)).and. &
            (hm(i,k).lt.hm(i,k-1)))then
           work(i)=hm(i,k)
           ihmin(i)=k
         endif
 190    continue
 200  continue
      do 210 i=1,len
        ihmin(i)=min(ihmin(i),nlm)
        if(ihmin(i).le.minorig)then
          iflag(i)=6
        endif
 210  continue
!
!-------------------------------------------------------------------
! --- Find that model level below the level of minimum moist static
! --- energy that has the maximum value of moist static energy
!-------------------------------------------------------------------

      do 220 i=1,len
       work(i)=hm(i,minorig)
       nk(i)=minorig
 220  continue
      do 240 k=minorig+1,nl
        do 230 i=1,len
         if((hm(i,k).gt.work(i)).and.(k.le.ihmin(i)))then
           work(i)=hm(i,k)
           nk(i)=k
         endif
 230     continue
 240  continue
!-------------------------------------------------------------------
! --- Check whether parcel level temperature and specific humidity
! --- are reasonable
!-------------------------------------------------------------------
       do 250 i=1,len
       if(((t(i,nk(i)).lt.250.0).or. &
            (q(i,nk(i)).le.0.0).or. &
            (p(i,ihmin(i)).lt.400.0)).and. &
            (iflag(i).eq.0))iflag(i)=7
 250   continue
!-------------------------------------------------------------------
! --- Calculate lifted condensation level of air at parcel origin level
! --- (Within 0.2% of formula of Bolton, MON. WEA. REV.,1980)
!-------------------------------------------------------------------
       do 260 i=1,len
        tnk(i)=t(i,nk(i))
        qnk(i)=q(i,nk(i))
        gznk(i)=gz(i,nk(i))
        pnk(i)=p(i,nk(i))
        qsnk(i)=qs(i,nk(i))
!
        rh(i)=qnk(i)/qsnk(i)
        rh(i)=min(1.0,rh(i))
        chi(i)=tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
        plcl(i)=pnk(i)*(rh(i)**chi(i))
        if(((plcl(i).lt.200.0).or.(plcl(i).ge.2000.0)) &
         .and.(iflag(i).eq.0))iflag(i)=8
 260   continue
!-------------------------------------------------------------------
! --- Calculate first level above lcl (=icb)
!-------------------------------------------------------------------
      do 270 i=1,len
       icb(i)=nlm
 270  continue
!
      do 290 k=minorig,nl
        do 280 i=1,len
          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i))) &
          icb(i)=min(icb(i),k)
 280    continue
 290  continue
!
      do 300 i=1,len
        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
 300  continue
!
! Compute icbmax.
!
      icbmax=2
      do 310 i=1,len
        icbmax=max(icbmax,icb(i))
 310  continue

      return
      end
