
      SUBROUTINE cv3_yield(nloc,ncum,nd,na,ntra  &
                          ,icb,inb,delt &
                          ,t,rr,u,v,tra,gz,p,ph,h,hp,lv,cpn,th &
                          ,ep,clw,m,tp,mp,rp,up,vp,trap &
                          ,wt,water,evap,b &
                          ,ment,qent,uent,vent,nent,elij,traent,sig &
                          ,tv,tvp &
                          ,iflag,precip,VPrecip,ft,fr,fu,fv,ftra &
                          ,upwd,dnwd,dnwd0,ma,mike,tls,tps,qcondc,wd)
      use conema3_m
            use cvparam3
            use cvthermo
            use cvflag
      implicit none


! inputs:
      integer ncum,nd,na,ntra,nloc
      integer icb(nloc), inb(nloc)
      real, intent(in):: delt
      real t(nloc,nd), rr(nloc,nd), u(nloc,nd), v(nloc,nd)
      real tra(nloc,nd,ntra), sig(nloc,nd)
      real gz(nloc,na), ph(nloc,nd+1), h(nloc,na), hp(nloc,na)
      real th(nloc,na), p(nloc,nd), tp(nloc,na)
      real lv(nloc,na), cpn(nloc,na), ep(nloc,na), clw(nloc,na)
      real m(nloc,na), mp(nloc,na), rp(nloc,na), up(nloc,na)
      real vp(nloc,na), wt(nloc,nd), trap(nloc,nd,ntra)
      real water(nloc,na), evap(nloc,na), b(nloc,na)
      real ment(nloc,na,na), qent(nloc,na,na), uent(nloc,na,na)
!ym      real vent(nloc,na,na), nent(nloc,na), elij(nloc,na,na)
      real vent(nloc,na,na), elij(nloc,na,na)
      integer nent(nloc,na)
      real traent(nloc,na,na,ntra)
      real tv(nloc,nd), tvp(nloc,nd)

! input/output:
      integer iflag(nloc)

! outputs:
      real precip(nloc)
      real VPrecip(nloc,nd+1)
      real ft(nloc,nd), fr(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real ftra(nloc,nd,ntra)
      real upwd(nloc,nd), dnwd(nloc,nd), ma(nloc,nd)
      real dnwd0(nloc,nd), mike(nloc,nd)
      real tls(nloc,nd), tps(nloc,nd)
      real qcondc(nloc,nd)                               ! cld
      real wd(nloc)                                      ! gust

! local variables:
      integer i,k,il,n,j,num1
      real rat, awat, delti
      real ax, bx, cx, dx, ex
      real cpinv, rdcp, dpinv
      real lvcp(nloc,na), mke(nloc,na)
      real am(nloc), work(nloc), ad(nloc), amp1(nloc)
!!!      real up1(nloc), dn1(nloc)
      real up1(nloc,nd,nd), dn1(nloc,nd,nd)
      real asum(nloc), bsum(nloc), csum(nloc), dsum(nloc)
      real qcond(nloc,nd), nqcond(nloc,nd), wa(nloc,nd)  ! cld
      real siga(nloc,nd), sax(nloc,nd), mac(nloc,nd)      ! cld


!-------------------------------------------------------------

! initialization:

      delti = 1.0/delt

      do il=1,ncum
       precip(il)=0.0
       wd(il)=0.0     ! gust
       VPrecip(il,nd+1)=0.
      enddo

      do i=1,nd
       do il=1,ncum
         VPrecip(il,i)=0.0
         ft(il,i)=0.0
         fr(il,i)=0.0
         fu(il,i)=0.0
         fv(il,i)=0.0
         qcondc(il,i)=0.0                                ! cld
         qcond(il,i)=0.0                                 ! cld
         nqcond(il,i)=0.0                                ! cld
       enddo
      enddo

!      do j=1,ntra
!       do i=1,nd
!        do il=1,ncum
!          ftra(il,i,j)=0.0
!        enddo
!       enddo
!      enddo

      do i=1,nl
       do il=1,ncum
         lvcp(il,i)=lv(il,i)/cpn(il,i)
       enddo
      enddo


!
!   ***  calculate surface precipitation in mm/day     ***
!
      do il=1,ncum
       if(ep(il,inb(il)).ge.0.0001)then
        if (cvflag_grav) then
         precip(il)=wt(il,1)*sigd*water(il,1)*86400.*1000./(rowl*grav)
        else
         precip(il)=wt(il,1)*sigd*water(il,1)*8640.
        endif
       endif
      enddo

!   ***  CALCULATE VERTICAL PROFILE OF  PRECIPITATIONs IN kg/m2/s  ===
!
! MAF rajout pour lessivage
       do k=1,nl
         do il=1,ncum
          if (k.le.inb(il)) then
            if (cvflag_grav) then
             VPrecip(il,k) = wt(il,k)*sigd*water(il,k)/grav
            else
             VPrecip(il,k) = wt(il,k)*sigd*water(il,k)/10.
            endif
          endif
         end do
       end do
!
!
!   ***  Calculate downdraft velocity scale    ***
!   ***  NE PAS UTILISER POUR L'INSTANT ***
!
!!      do il=1,ncum
!!        wd(il)=betad*abs(mp(il,icb(il)))*0.01*rrd*t(il,icb(il))
!!     :                                  /(sigd*p(il,icb(il)))
!!      enddo

!
!   ***  calculate tendencies of lowest level potential temperature  ***
!   ***                      and mixing ratio                        ***
!
      do il=1,ncum
       work(il)=1.0/(ph(il,1)-ph(il,2))
       am(il)=0.0
      enddo

      do k=2,nl
       do il=1,ncum
        if (k.le.inb(il)) then
         am(il)=am(il)+m(il,k)
        endif
       enddo
      enddo

      do il=1,ncum

! convect3      if((0.1*dpinv*am).ge.delti)iflag(il)=4
      if (cvflag_grav) then
      if((0.01*grav*work(il)*am(il)).ge.delti)iflag(il)=1!consist vect
       ft(il,1)=0.01*grav*work(il)*am(il)*(t(il,2)-t(il,1) &
                  +(gz(il,2)-gz(il,1))/cpn(il,1))
      else
       if((0.1*work(il)*am(il)).ge.delti)iflag(il)=1 !consistency vect
       ft(il,1)=0.1*work(il)*am(il)*(t(il,2)-t(il,1) &
                  +(gz(il,2)-gz(il,1))/cpn(il,1))
      endif

      ft(il,1)=ft(il,1)-0.5*lvcp(il,1)*sigd*(evap(il,1)+evap(il,2))

      if (cvflag_grav) then
       ft(il,1)=ft(il,1)-0.009*grav*sigd*mp(il,2) &
                                   *t(il,1)*b(il,1)*work(il)
      else
       ft(il,1)=ft(il,1)-0.09*sigd*mp(il,2)*t(il,1)*b(il,1)*work(il)
      endif

      ft(il,1)=ft(il,1)+0.01*sigd*wt(il,1)*(cl-cpd)*water(il,2)*(t(il,2) &
      -t(il,1))*work(il)/cpn(il,1)

      if (cvflag_grav) then
!jyg1  Correction pour mieux conserver l'eau (conformite avec CONVECT4.3)
! (sb: pour l'instant, on ne fait que le chgt concernant grav, pas evap)
       fr(il,1)=0.01*grav*mp(il,2)*(rp(il,2)-rr(il,1))*work(il) &
                +sigd*0.5*(evap(il,1)+evap(il,2))
!+tard     :          +sigd*evap(il,1)

       fr(il,1)=fr(il,1)+0.01*grav*am(il)*(rr(il,2)-rr(il,1))*work(il)

       fu(il,1)=fu(il,1)+0.01*grav*work(il)*(mp(il,2)*(up(il,2)-u(il,1)) &
               +am(il)*(u(il,2)-u(il,1)))
       fv(il,1)=fv(il,1)+0.01*grav*work(il)*(mp(il,2)*(vp(il,2)-v(il,1)) &
               +am(il)*(v(il,2)-v(il,1)))
      else  ! cvflag_grav
       fr(il,1)=0.1*mp(il,2)*(rp(il,2)-rr(il,1))*work(il) &
                +sigd*0.5*(evap(il,1)+evap(il,2))
       fr(il,1)=fr(il,1)+0.1*am(il)*(rr(il,2)-rr(il,1))*work(il)
       fu(il,1)=fu(il,1)+0.1*work(il)*(mp(il,2)*(up(il,2)-u(il,1)) &
               +am(il)*(u(il,2)-u(il,1)))
       fv(il,1)=fv(il,1)+0.1*work(il)*(mp(il,2)*(vp(il,2)-v(il,1)) &
               +am(il)*(v(il,2)-v(il,1)))
      endif ! cvflag_grav

      enddo ! il

!      do j=1,ntra
!       do il=1,ncum
!        if (cvflag_grav) then
!         ftra(il,1,j)=ftra(il,1,j)+0.01*grav*work(il)
!     :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
!     :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
!        else
!         ftra(il,1,j)=ftra(il,1,j)+0.1*work(il)
!     :                     *(mp(il,2)*(trap(il,2,j)-tra(il,1,j))
!     :             +am(il)*(tra(il,2,j)-tra(il,1,j)))
!        endif
!       enddo
!      enddo

      do j=2,nl
       do il=1,ncum
        if (j.le.inb(il)) then
         if (cvflag_grav) then
          fr(il,1)=fr(il,1) &
             +0.01*grav*work(il)*ment(il,j,1)*(qent(il,j,1)-rr(il,1))
          fu(il,1)=fu(il,1) &
             +0.01*grav*work(il)*ment(il,j,1)*(uent(il,j,1)-u(il,1))
          fv(il,1)=fv(il,1) &
             +0.01*grav*work(il)*ment(il,j,1)*(vent(il,j,1)-v(il,1))
         else   ! cvflag_grav
          fr(il,1)=fr(il,1) &
               +0.1*work(il)*ment(il,j,1)*(qent(il,j,1)-rr(il,1))
          fu(il,1)=fu(il,1) &
               +0.1*work(il)*ment(il,j,1)*(uent(il,j,1)-u(il,1))
          fv(il,1)=fv(il,1) &
               +0.1*work(il)*ment(il,j,1)*(vent(il,j,1)-v(il,1))
         endif  ! cvflag_grav
        endif ! j
       enddo
      enddo

!      do k=1,ntra
!       do j=2,nl
!        do il=1,ncum
!         if (j.le.inb(il)) then

!          if (cvflag_grav) then
!           ftra(il,1,k)=ftra(il,1,k)+0.01*grav*work(il)*ment(il,j,1)
!     :                *(traent(il,j,1,k)-tra(il,1,k))
!          else
!           ftra(il,1,k)=ftra(il,1,k)+0.1*work(il)*ment(il,j,1)
!     :                *(traent(il,j,1,k)-tra(il,1,k))
!          endif

!         endif
!        enddo
!       enddo
!      enddo

!
!   ***  calculate tendencies of potential temperature and mixing ratio  ***
!   ***               at levels above the lowest level                   ***
!
!   ***  first find the net saturated updraft and downdraft mass fluxes  ***
!   ***                      through each level                          ***
!

      do 500 i=2,nl+1 ! newvecto: mettre nl au lieu nl+1?

       num1=0
       do il=1,ncum
        if(i.le.inb(il))num1=num1+1
       enddo
       if(num1.le.0)go to 500

       call zilch(amp1,ncum)
       call zilch(ad,ncum)

      do 440 k=i+1,nl+1
       do 441 il=1,ncum
        if (i.le.inb(il) .and. k.le.(inb(il)+1)) then
         amp1(il)=amp1(il)+m(il,k)
        endif
 441   continue
 440  continue

      do 450 k=1,i
       do 451 j=i+1,nl+1
        do 452 il=1,ncum
         if (i.le.inb(il) .and. j.le.(inb(il)+1)) then
          amp1(il)=amp1(il)+ment(il,k,j)
         endif
452     continue
451    continue
450   continue

      do 470 k=1,i-1
       do 471 j=i,nl+1 ! newvecto: nl au lieu nl+1?
        do 472 il=1,ncum
        if (i.le.inb(il) .and. j.le.inb(il)) then
         ad(il)=ad(il)+ment(il,j,k)
        endif
472     continue
471    continue
470   continue

      do 1350 il=1,ncum
      if (i.le.inb(il)) then
       dpinv=1.0/(ph(il,i)-ph(il,i+1))
       cpinv=1.0/cpn(il,i)

! convect3      if((0.1*dpinv*amp1).ge.delti)iflag(il)=4
      if (cvflag_grav) then
       if((0.01*grav*dpinv*amp1(il)).ge.delti)iflag(il)=1 ! vecto
      else
       if((0.1*dpinv*amp1(il)).ge.delti)iflag(il)=1 ! vecto
      endif

      if (cvflag_grav) then
       ft(il,i)=0.01*grav*dpinv*(amp1(il)*(t(il,i+1)-t(il,i) &
          +(gz(il,i+1)-gz(il,i))*cpinv) &
          -ad(il)*(t(il,i)-t(il,i-1)+(gz(il,i)-gz(il,i-1))*cpinv)) &
          -0.5*sigd*lvcp(il,i)*(evap(il,i)+evap(il,i+1))
       rat=cpn(il,i-1)*cpinv
       ft(il,i)=ft(il,i)-0.009*grav*sigd*(mp(il,i+1)*t(il,i)*b(il,i) &
         -mp(il,i)*t(il,i-1)*rat*b(il,i-1))*dpinv
       ft(il,i)=ft(il,i)+0.01*grav*dpinv*ment(il,i,i)*(hp(il,i)-h(il,i) &
          +t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
      else  ! cvflag_grav
       ft(il,i)=0.1*dpinv*(amp1(il)*(t(il,i+1)-t(il,i) &
          +(gz(il,i+1)-gz(il,i))*cpinv) &
          -ad(il)*(t(il,i)-t(il,i-1)+(gz(il,i)-gz(il,i-1))*cpinv)) &
          -0.5*sigd*lvcp(il,i)*(evap(il,i)+evap(il,i+1))
       rat=cpn(il,i-1)*cpinv
       ft(il,i)=ft(il,i)-0.09*sigd*(mp(il,i+1)*t(il,i)*b(il,i) &
         -mp(il,i)*t(il,i-1)*rat*b(il,i-1))*dpinv
       ft(il,i)=ft(il,i)+0.1*dpinv*ment(il,i,i)*(hp(il,i)-h(il,i) &
          +t(il,i)*(cpv-cpd)*(rr(il,i)-qent(il,i,i)))*cpinv
      endif ! cvflag_grav


      ft(il,i)=ft(il,i)+0.01*sigd*wt(il,i)*(cl-cpd)*water(il,i+1) &
                 *(t(il,i+1)-t(il,i))*dpinv*cpinv

      if (cvflag_grav) then
       fr(il,i)=0.01*grav*dpinv*(amp1(il)*(rr(il,i+1)-rr(il,i)) &
                 -ad(il)*(rr(il,i)-rr(il,i-1)))
       fu(il,i)=fu(il,i)+0.01*grav*dpinv*(amp1(il)*(u(il,i+1)-u(il,i)) &
                   -ad(il)*(u(il,i)-u(il,i-1)))
       fv(il,i)=fv(il,i)+0.01*grav*dpinv*(amp1(il)*(v(il,i+1)-v(il,i)) &
                   -ad(il)*(v(il,i)-v(il,i-1)))
      else  ! cvflag_grav
       fr(il,i)=0.1*dpinv*(amp1(il)*(rr(il,i+1)-rr(il,i)) &
                 -ad(il)*(rr(il,i)-rr(il,i-1)))
       fu(il,i)=fu(il,i)+0.1*dpinv*(amp1(il)*(u(il,i+1)-u(il,i)) &
                   -ad(il)*(u(il,i)-u(il,i-1)))
       fv(il,i)=fv(il,i)+0.1*dpinv*(amp1(il)*(v(il,i+1)-v(il,i)) &
                   -ad(il)*(v(il,i)-v(il,i-1)))
      endif ! cvflag_grav

      endif ! i
1350  continue

!      do k=1,ntra
!       do il=1,ncum
!        if (i.le.inb(il)) then
!         dpinv=1.0/(ph(il,i)-ph(il,i+1))
!         cpinv=1.0/cpn(il,i)
!         if (cvflag_grav) then
!           ftra(il,i,k)=ftra(il,i,k)+0.01*grav*dpinv
!     :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
!     :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
!         else
!           ftra(il,i,k)=ftra(il,i,k)+0.1*dpinv
!     :         *(amp1(il)*(tra(il,i+1,k)-tra(il,i,k))
!     :           -ad(il)*(tra(il,i,k)-tra(il,i-1,k)))
!         endif
!        endif
!       enddo
!      enddo

      do 480 k=1,i-1
       do 1370 il=1,ncum
        if (i.le.inb(il)) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)

      awat=elij(il,k,i)-(1.-ep(il,i))*clw(il,i)
      awat=amax1(awat,0.0)

      if (cvflag_grav) then
      fr(il,i)=fr(il,i) &
         +0.01*grav*dpinv*ment(il,k,i)*(qent(il,k,i)-awat-rr(il,i))
      fu(il,i)=fu(il,i) &
               +0.01*grav*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
      fv(il,i)=fv(il,i) &
               +0.01*grav*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
      else  ! cvflag_grav
      fr(il,i)=fr(il,i) &
         +0.1*dpinv*ment(il,k,i)*(qent(il,k,i)-awat-rr(il,i))
      fu(il,i)=fu(il,i) &
         +0.01*grav*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
      fv(il,i)=fv(il,i) &
         +0.1*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
      endif ! cvflag_grav

! (saturated updrafts resulting from mixing)        ! cld
        qcond(il,i)=qcond(il,i)+(elij(il,k,i)-awat) ! cld
        nqcond(il,i)=nqcond(il,i)+1.                ! cld
      endif ! i
1370  continue
480   continue

!      do j=1,ntra
!       do k=1,i-1
!        do il=1,ncum
!         if (i.le.inb(il)) then
!          dpinv=1.0/(ph(il,i)-ph(il,i+1))
!          cpinv=1.0/cpn(il,i)
!          if (cvflag_grav) then
!           ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
!     :        *(traent(il,k,i,j)-tra(il,i,j))
!          else
!           ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
!     :        *(traent(il,k,i,j)-tra(il,i,j))
!          endif
!         endif
!        enddo
!       enddo
!      enddo

      do 490 k=i,nl+1
       do 1380 il=1,ncum
        if (i.le.inb(il) .and. k.le.inb(il)) then
         dpinv=1.0/(ph(il,i)-ph(il,i+1))
         cpinv=1.0/cpn(il,i)

         if (cvflag_grav) then
         fr(il,i)=fr(il,i) &
               +0.01*grav*dpinv*ment(il,k,i)*(qent(il,k,i)-rr(il,i))
         fu(il,i)=fu(il,i) &
               +0.01*grav*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
         fv(il,i)=fv(il,i) &
               +0.01*grav*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
         else  ! cvflag_grav
         fr(il,i)=fr(il,i) &
               +0.1*dpinv*ment(il,k,i)*(qent(il,k,i)-rr(il,i))
         fu(il,i)=fu(il,i) &
               +0.1*dpinv*ment(il,k,i)*(uent(il,k,i)-u(il,i))
         fv(il,i)=fv(il,i) &
               +0.1*dpinv*ment(il,k,i)*(vent(il,k,i)-v(il,i))
         endif ! cvflag_grav
        endif ! i and k
1380   continue
490   continue

!      do j=1,ntra
!       do k=i,nl+1
!        do il=1,ncum
!         if (i.le.inb(il) .and. k.le.inb(il)) then
!          dpinv=1.0/(ph(il,i)-ph(il,i+1))
!          cpinv=1.0/cpn(il,i)
!          if (cvflag_grav) then
!           ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv*ment(il,k,i)
!     :         *(traent(il,k,i,j)-tra(il,i,j))
!          else
!           ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv*ment(il,k,i)
!     :             *(traent(il,k,i,j)-tra(il,i,j))
!          endif
!         endif ! i and k
!        enddo
!       enddo
!      enddo

      do 1400 il=1,ncum
       if (i.le.inb(il)) then
        dpinv=1.0/(ph(il,i)-ph(il,i+1))
        cpinv=1.0/cpn(il,i)

        if (cvflag_grav) then
! sb: on ne fait pas encore la correction permettant de mieux
! conserver l'eau:
         fr(il,i)=fr(il,i)+0.5*sigd*(evap(il,i)+evap(il,i+1)) &
              +0.01*grav*(mp(il,i+1)*(rp(il,i+1)-rr(il,i))-mp(il,i) &
                     *(rp(il,i)-rr(il,i-1)))*dpinv

         fu(il,i)=fu(il,i)+0.01*grav*(mp(il,i+1)*(up(il,i+1)-u(il,i)) &
                   -mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
         fv(il,i)=fv(il,i)+0.01*grav*(mp(il,i+1)*(vp(il,i+1)-v(il,i)) &
                   -mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv
        else  ! cvflag_grav
         fr(il,i)=fr(il,i)+0.5*sigd*(evap(il,i)+evap(il,i+1)) &
              +0.1*(mp(il,i+1)*(rp(il,i+1)-rr(il,i))-mp(il,i) &
                     *(rp(il,i)-rr(il,i-1)))*dpinv
         fu(il,i)=fu(il,i)+0.1*(mp(il,i+1)*(up(il,i+1)-u(il,i)) &
                   -mp(il,i)*(up(il,i)-u(il,i-1)))*dpinv
         fv(il,i)=fv(il,i)+0.1*(mp(il,i+1)*(vp(il,i+1)-v(il,i)) &
                   -mp(il,i)*(vp(il,i)-v(il,i-1)))*dpinv
        endif ! cvflag_grav

      endif ! i
1400  continue

! sb: interface with the cloud parameterization:          ! cld

      do k=i+1,nl
       do il=1,ncum
        if (k.le.inb(il) .and. i.le.inb(il)) then         ! cld
! (saturated downdrafts resulting from mixing)            ! cld
          qcond(il,i)=qcond(il,i)+elij(il,k,i)            ! cld
          nqcond(il,i)=nqcond(il,i)+1.                    ! cld
        endif                                             ! cld
       enddo                                              ! cld
      enddo                                               ! cld

! (particular case: no detraining level is found)         ! cld
      do il=1,ncum                                        ! cld
       if (i.le.inb(il) .and. nent(il,i).eq.0) then       ! cld
          qcond(il,i)=qcond(il,i)+(1.-ep(il,i))*clw(il,i) ! cld
          nqcond(il,i)=nqcond(il,i)+1.                    ! cld
       endif                                              ! cld
      enddo                                               ! cld

      do il=1,ncum                                        ! cld
       if (i.le.inb(il) .and. nqcond(il,i).ne.0.) then    ! cld
          qcond(il,i)=qcond(il,i)/nqcond(il,i)            ! cld
       endif                                              ! cld
      enddo

!      do j=1,ntra
!       do il=1,ncum
!        if (i.le.inb(il)) then
!         dpinv=1.0/(ph(il,i)-ph(il,i+1))
!         cpinv=1.0/cpn(il,i)

!         if (cvflag_grav) then
!          ftra(il,i,j)=ftra(il,i,j)+0.01*grav*dpinv
!     :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
!     :     -mp(il,i)*(trap(il,i,j)-tra(il,i-1,j)))
!         else
!          ftra(il,i,j)=ftra(il,i,j)+0.1*dpinv
!     :     *(mp(il,i+1)*(trap(il,i+1,j)-tra(il,i,j))
!     :     -mp(il,i)*(trap(il,i,j)-tra(il,i-1,j)))
!         endif
!        endif ! i
!       enddo
!      enddo

500   continue


!   ***   move the detrainment at level inb down to level inb-1   ***
!   ***        in such a way as to preserve the vertically        ***
!   ***          integrated enthalpy and water tendencies         ***
!
      do 503 il=1,ncum

      ax=0.1*ment(il,inb(il),inb(il))*(hp(il,inb(il))-h(il,inb(il)) &
       +t(il,inb(il))*(cpv-cpd) &
       *(rr(il,inb(il))-qent(il,inb(il),inb(il)))) &
        /(cpn(il,inb(il))*(ph(il,inb(il))-ph(il,inb(il)+1)))
      ft(il,inb(il))=ft(il,inb(il))-ax
      ft(il,inb(il)-1)=ft(il,inb(il)-1)+ax*cpn(il,inb(il)) &
          *(ph(il,inb(il))-ph(il,inb(il)+1))/(cpn(il,inb(il)-1) &
          *(ph(il,inb(il)-1)-ph(il,inb(il))))

      bx=0.1*ment(il,inb(il),inb(il))*(qent(il,inb(il),inb(il)) &
          -rr(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fr(il,inb(il))=fr(il,inb(il))-bx
      fr(il,inb(il)-1)=fr(il,inb(il)-1) &
         +bx*(ph(il,inb(il))-ph(il,inb(il)+1)) &
            /(ph(il,inb(il)-1)-ph(il,inb(il)))

      cx=0.1*ment(il,inb(il),inb(il))*(uent(il,inb(il),inb(il)) &
             -u(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fu(il,inb(il))=fu(il,inb(il))-cx
      fu(il,inb(il)-1)=fu(il,inb(il)-1) &
           +cx*(ph(il,inb(il))-ph(il,inb(il)+1)) &
              /(ph(il,inb(il)-1)-ph(il,inb(il)))

      dx=0.1*ment(il,inb(il),inb(il))*(vent(il,inb(il),inb(il)) &
            -v(il,inb(il)))/(ph(il,inb(il))-ph(il,inb(il)+1))
      fv(il,inb(il))=fv(il,inb(il))-dx
      fv(il,inb(il)-1)=fv(il,inb(il)-1) &
          +dx*(ph(il,inb(il))-ph(il,inb(il)+1)) &
             /(ph(il,inb(il)-1)-ph(il,inb(il)))

503   continue

!      do j=1,ntra
!       do il=1,ncum
!        ex=0.1*ment(il,inb(il),inb(il))
!     :      *(traent(il,inb(il),inb(il),j)-tra(il,inb(il),j))
!     :      /(ph(il,inb(il))-ph(il,inb(il)+1))
!        ftra(il,inb(il),j)=ftra(il,inb(il),j)-ex
!        ftra(il,inb(il)-1,j)=ftra(il,inb(il)-1,j)
!     :       +ex*(ph(il,inb(il))-ph(il,inb(il)+1))
!     :          /(ph(il,inb(il)-1)-ph(il,inb(il)))
!       enddo
!      enddo

!
!   ***    homoginize tendencies below cloud base    ***
!
!
      do il=1,ncum
       asum(il)=0.0
       bsum(il)=0.0
       csum(il)=0.0
       dsum(il)=0.0
      enddo

      do i=1,nl
       do il=1,ncum
        if (i.le.(icb(il)-1)) then
      asum(il)=asum(il)+ft(il,i)*(ph(il,i)-ph(il,i+1))
      bsum(il)=bsum(il)+fr(il,i)*(lv(il,i)+(cl-cpd)*(t(il,i)-t(il,1))) &
                        *(ph(il,i)-ph(il,i+1))
      csum(il)=csum(il)+(lv(il,i)+(cl-cpd)*(t(il,i)-t(il,1))) &
                            *(ph(il,i)-ph(il,i+1))
      dsum(il)=dsum(il)+t(il,i)*(ph(il,i)-ph(il,i+1))/th(il,i)
        endif
       enddo
      enddo

!!!!      do 700 i=1,icb(il)-1
      do i=1,nl
       do il=1,ncum
        if (i.le.(icb(il)-1)) then
         ft(il,i)=asum(il)*t(il,i)/(th(il,i)*dsum(il))
         fr(il,i)=bsum(il)/csum(il)
        endif
       enddo
      enddo

!
!   ***           reset counter and return           ***
!
      do il=1,ncum
       sig(il,nd)=2.0
      enddo


      do i=1,nd
       do il=1,ncum
        upwd(il,i)=0.0
        dnwd(il,i)=0.0
       enddo
      enddo

      do i=1,nl
       do il=1,ncum
        dnwd0(il,i)=-mp(il,i)
       enddo
      enddo
      do i=nl+1,nd
       do il=1,ncum
        dnwd0(il,i)=0.
       enddo
      enddo


      do i=1,nl
       do il=1,ncum
        if (i.ge.icb(il) .and. i.le.inb(il)) then
          upwd(il,i)=0.0
          dnwd(il,i)=0.0
        endif
       enddo
      enddo

      do i=1,nl
       do k=1,nl
        do il=1,ncum
          up1(il,k,i)=0.0
          dn1(il,k,i)=0.0
        enddo
       enddo
      enddo

      do i=1,nl
       do k=i,nl
        do n=1,i-1
         do il=1,ncum
          if (i.ge.icb(il).and.i.le.inb(il).and.k.le.inb(il)) then
             up1(il,k,i)=up1(il,k,i)+ment(il,n,k)
             dn1(il,k,i)=dn1(il,k,i)-ment(il,k,n)
          endif
         enddo
        enddo
       enddo
      enddo

      do i=2,nl
       do k=i,nl
        do il=1,ncum
!test         if (i.ge.icb(il).and.i.le.inb(il).and.k.le.inb(il)) then
         if (i.le.inb(il).and.k.le.inb(il)) then
            upwd(il,i)=upwd(il,i)+m(il,k)+up1(il,k,i)
            dnwd(il,i)=dnwd(il,i)+dn1(il,k,i)
         endif
        enddo
       enddo
      enddo


!!!!      DO il=1,ncum
!!!!      do i=icb(il),inb(il)
!!!!
!!!!      upwd(il,i)=0.0
!!!!      dnwd(il,i)=0.0
!!!!      do k=i,inb(il)
!!!!      up1=0.0
!!!!      dn1=0.0
!!!!      do n=1,i-1
!!!!      up1=up1+ment(il,n,k)
!!!!      dn1=dn1-ment(il,k,n)
!!!!      enddo
!!!!      upwd(il,i)=upwd(il,i)+m(il,k)+up1
!!!!      dnwd(il,i)=dnwd(il,i)+dn1
!!!!      enddo
!!!!      enddo
!!!!
!!!!      ENDDO

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        determination de la variation de flux ascendant entre
!        deux niveau non dilue mike
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nl
       do il=1,ncum
        mike(il,i)=m(il,i)
       enddo
      enddo

      do i=nl+1,nd
       do il=1,ncum
        mike(il,i)=0.
       enddo
      enddo

      do i=1,nd
       do il=1,ncum
        ma(il,i)=0
       enddo
      enddo

      do i=1,nl
       do j=i,nl
        do il=1,ncum
         ma(il,i)=ma(il,i)+m(il,j)
        enddo
       enddo
      enddo

      do i=nl+1,nd
       do il=1,ncum
        ma(il,i)=0.
       enddo
      enddo

      do i=1,nl
       do il=1,ncum
        if (i.le.(icb(il)-1)) then
         ma(il,i)=0
        endif
       enddo
      enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        icb represente de niveau ou se trouve la
!        base du nuage , et inb le top du nuage
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nd
       do il=1,ncum
        mke(il,i)=upwd(il,i)+dnwd(il,i)
       enddo
      enddo

      do i=1,nd
       DO 999 il=1,ncum
        rdcp=(rrd*(1.-rr(il,i))-rr(il,i)*rrv) &
              /(cpd*(1.-rr(il,i))+rr(il,i)*cpv)
        tls(il,i)=t(il,i)*(1000.0/p(il,i))**rdcp
        tps(il,i)=tp(il,i)
999    CONTINUE
      enddo

!
!   *** diagnose the in-cloud mixing ratio   ***            ! cld
!   ***           of condensed water         ***            ! cld
!                                                           ! cld

       do i=1,nd                                            ! cld
        do il=1,ncum                                        ! cld
         mac(il,i)=0.0                                      ! cld
         wa(il,i)=0.0                                       ! cld
         siga(il,i)=0.0                                     ! cld
         sax(il,i)=0.0                                      ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=minorig, nl                                     ! cld
        do k=i+1,nl+1                                       ! cld
         do il=1,ncum                                       ! cld
          if (i.le.inb(il) .and. k.le.(inb(il)+1)) then     ! cld
            mac(il,i)=mac(il,i)+m(il,k)                     ! cld
          endif                                             ! cld
         enddo                                              ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=1,nl                                            ! cld
        do j=1,i                                            ! cld
         do il=1,ncum                                       ! cld
          if (i.ge.icb(il) .and. i.le.(inb(il)-1) &
            .and. j.ge.icb(il) ) then                       ! cld
           sax(il,i)=sax(il,i)+rrd*(tvp(il,j)-tv(il,j)) &
              *(ph(il,j)-ph(il,j+1))/p(il,j)                ! cld
          endif                                             ! cld
         enddo                                              ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=1,nl                                            ! cld
        do il=1,ncum                                        ! cld
         if (i.ge.icb(il) .and. i.le.(inb(il)-1) &
             .and. sax(il,i).gt.0.0 ) then                  ! cld
           wa(il,i)=sqrt(2.*sax(il,i))                      ! cld
         endif                                              ! cld
        enddo                                               ! cld
       enddo                                                ! cld

       do i=1,nl                                            ! cld
        do il=1,ncum                                        ! cld
         if (wa(il,i).gt.0.0) &
           siga(il,i)=mac(il,i)/wa(il,i) &
               *rrd*tvp(il,i)/p(il,i)/100./delta            ! cld
          siga(il,i) = min(siga(il,i),1.0)                  ! cld
!IM cf. FH
         if (iflag_clw.eq.0) then
          qcondc(il,i)=siga(il,i)*clw(il,i)*(1.-ep(il,i)) &
                 + (1.-siga(il,i))*qcond(il,i)              ! cld
         else if (iflag_clw.eq.1) then
          qcondc(il,i)=qcond(il,i)              ! cld
         endif

        enddo                                               ! cld
       enddo                                                ! cld

        return
        end
