
      SUBROUTINE cv_yield(nloc,ncum,nd,nk,icb,inb,delt &
                   ,t,q,u,v,gz,p,ph,h,hp,lv,cpn &
                   ,ep,clw,frac,m,mp,qp,up,vp &
                   ,wt,water,evap &
                   ,ment,qent,uent,vent,nent,elij &
                   ,tv,tvp &
                   ,iflag,wd,qprime,tprime &
                   ,precip,cbmf,ft,fq,fu,fv,Ma,qcondc)
            use cvthermo
            use cvparam
      implicit none


! inputs
      integer ncum, nd, nloc
      integer nk(nloc), icb(nloc), inb(nloc)
      integer nent(nloc,nd)
      real, intent(in):: delt
      real t(nloc,nd), q(nloc,nd), u(nloc,nd), v(nloc,nd)
      real gz(nloc,nd)
      real p(nloc,nd), ph(nloc,nd+1), h(nloc,nd)
      real hp(nloc,nd), lv(nloc,nd)
      real cpn(nloc,nd), ep(nloc,nd), clw(nloc,nd), frac(nloc)
      real m(nloc,nd), mp(nloc,nd), qp(nloc,nd)
      real up(nloc,nd), vp(nloc,nd)
      real wt(nloc,nd), water(nloc,nd), evap(nloc,nd)
      real ment(nloc,nd,nd), qent(nloc,nd,nd), elij(nloc,nd,nd)
      real uent(nloc,nd,nd), vent(nloc,nd,nd)
      real tv(nloc,nd), tvp(nloc,nd)

! outputs
      integer iflag(nloc)  ! also an input
      real cbmf(nloc)      ! also an input
      real wd(nloc), tprime(nloc), qprime(nloc)
      real precip(nloc)
      real ft(nloc,nd), fq(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real Ma(nloc,nd)
      real qcondc(nloc,nd)

! local variables
      integer i,j,ij,k,num1
      real dpinv,cpinv,awat,fqold,ftold,fuold,fvold,delti
      real work(nloc), am(nloc),amp1(nloc),ad(nloc)
      real ents(nloc), uav(nloc),vav(nloc),lvcp(nloc,nd)
      real qcond(nloc,nd), nqcond(nloc,nd), wa(nloc,nd) ! cld
      real siga(nloc,nd), ax(nloc,nd), mac(nloc,nd)     ! cld


! -- initializations:

      delti = 1.0/delt

      do 160 i=1,ncum
      precip(i)=0.0
      wd(i)=0.0
      tprime(i)=0.0
      qprime(i)=0.0
       do 170 k=1,nl+1
        ft(i,k)=0.0
        fu(i,k)=0.0
        fv(i,k)=0.0
        fq(i,k)=0.0
        lvcp(i,k)=lv(i,k)/cpn(i,k)
        qcondc(i,k)=0.0              ! cld
        qcond(i,k)=0.0               ! cld
        nqcond(i,k)=0.0              ! cld
 170   continue
 160  continue

!
!   ***  Calculate surface precipitation in mm/day     ***
!
        do 1190 i=1,ncum
          if(iflag(i).le.1)then
            precip(i) = wt(i,1)*sigd*water(i,1)*86400/g
          endif
 1190   continue
!
!
!   ***  Calculate downdraft velocity scale and surface temperature and  ***
!   ***                    water vapor fluctuations                      ***
!
      do i=1,ncum
       wd(i)=betad*abs(mp(i,icb(i)))*0.01*rrd*t(i,icb(i)) &
                 /(sigd*p(i,icb(i)))
       qprime(i)=0.5*(qp(i,1)-q(i,1))
       tprime(i)=lv0*qprime(i)/cpd
      enddo
!
!   ***  Calculate tendencies of lowest level potential temperature  ***
!   ***                      and mixing ratio                        ***
!
        do 1200 i=1,ncum
          work(i)=0.01/(ph(i,1)-ph(i,2))
          am(i)=0.0
 1200   continue
        do 1220 k=2,nl
          do 1210 i=1,ncum
            if((nk(i).eq.1).and.(k.le.inb(i)).and.(nk(i).eq.1))then
              am(i)=am(i)+m(i,k)
            endif
 1210     continue
 1220   continue
        do 1240 i=1,ncum
          if((g*work(i)*am(i)).ge.delti)iflag(i)=1
          ft(i,1)=ft(i,1)+g*work(i)*am(i)*(t(i,2)-t(i,1) &
          +(gz(i,2)-gz(i,1))/cpn(i,1))
          ft(i,1)=ft(i,1)-lvcp(i,1)*sigd*evap(i,1)
          ft(i,1)=ft(i,1)+sigd*wt(i,2)*(cl-cpd)*water(i,2)*(t(i,2) &
           -t(i,1))*work(i)/cpn(i,1)
          fq(i,1)=fq(i,1)+g*mp(i,2)*(qp(i,2)-q(i,1))* &
          work(i)+sigd*evap(i,1)
          fq(i,1)=fq(i,1)+g*am(i)*(q(i,2)-q(i,1))*work(i)
          fu(i,1)=fu(i,1)+g*work(i)*(mp(i,2)*(up(i,2)-u(i,1)) &
          +am(i)*(u(i,2)-u(i,1)))
          fv(i,1)=fv(i,1)+g*work(i)*(mp(i,2)*(vp(i,2)-v(i,1)) &
          +am(i)*(v(i,2)-v(i,1)))
 1240   continue
        do 1260 j=2,nl
           do 1250 i=1,ncum
             if(j.le.inb(i))then
               fq(i,1)=fq(i,1) &
                       +g*work(i)*ment(i,j,1)*(qent(i,j,1)-q(i,1))
               fu(i,1)=fu(i,1) &
                       +g*work(i)*ment(i,j,1)*(uent(i,j,1)-u(i,1))
               fv(i,1)=fv(i,1) &
                       +g*work(i)*ment(i,j,1)*(vent(i,j,1)-v(i,1))
             endif
 1250      continue
 1260   continue
!
!   ***  Calculate tendencies of potential temperature and mixing ratio  ***
!   ***               at levels above the lowest level                   ***
!
!   ***  First find the net saturated updraft and downdraft mass fluxes  ***
!   ***                      through each level                          ***
!
        do 1500 i=2,nl+1
!
          num1=0
          do 1265 ij=1,ncum
            if(i.le.inb(ij))num1=num1+1
 1265     continue
          if(num1.le.0)go to 1500
!
          call zilch(amp1,ncum)
          call zilch(ad,ncum)
!
          do 1280 k=i+1,nl+1
            do 1270 ij=1,ncum
              if((i.ge.nk(ij)).and.(i.le.inb(ij)) &
                  .and.(k.le.(inb(ij)+1)))then
                amp1(ij)=amp1(ij)+m(ij,k)
              endif
 1270         continue
 1280     continue
!
          do 1310 k=1,i
            do 1300 j=i+1,nl+1
               do 1290 ij=1,ncum
                 if((j.le.(inb(ij)+1)).and.(i.le.inb(ij)))then
                   amp1(ij)=amp1(ij)+ment(ij,k,j)
                 endif
 1290          continue
 1300       continue
 1310     continue
          do 1340 k=1,i-1
            do 1330 j=i,nl+1
              do 1320 ij=1,ncum
                if((i.le.inb(ij)).and.(j.le.inb(ij)))then
                   ad(ij)=ad(ij)+ment(ij,j,k)
                endif
 1320         continue
 1330       continue
 1340     continue
!
          do 1350 ij=1,ncum
          if(i.le.inb(ij))then
            dpinv=0.01/(ph(ij,i)-ph(ij,i+1))
            cpinv=1.0/cpn(ij,i)
!
            ft(ij,i)=ft(ij,i) &
             +g*dpinv*(amp1(ij)*(t(ij,i+1)-t(ij,i) &
             +(gz(ij,i+1)-gz(ij,i))*cpinv) &
             -ad(ij)*(t(ij,i)-t(ij,i-1)+(gz(ij,i)-gz(ij,i-1))*cpinv)) &
             -sigd*lvcp(ij,i)*evap(ij,i)
            ft(ij,i)=ft(ij,i)+g*dpinv*ment(ij,i,i)*(hp(ij,i)-h(ij,i)+ &
              t(ij,i)*(cpv-cpd)*(q(ij,i)-qent(ij,i,i)))*cpinv
            ft(ij,i)=ft(ij,i)+sigd*wt(ij,i+1)*(cl-cpd)*water(ij,i+1)* &
              (t(ij,i+1)-t(ij,i))*dpinv*cpinv
            fq(ij,i)=fq(ij,i)+g*dpinv*(amp1(ij)*(q(ij,i+1)-q(ij,i))- &
              ad(ij)*(q(ij,i)-q(ij,i-1)))
            fu(ij,i)=fu(ij,i)+g*dpinv*(amp1(ij)*(u(ij,i+1)-u(ij,i))- &
              ad(ij)*(u(ij,i)-u(ij,i-1)))
            fv(ij,i)=fv(ij,i)+g*dpinv*(amp1(ij)*(v(ij,i+1)-v(ij,i))- &
              ad(ij)*(v(ij,i)-v(ij,i-1)))
         endif
 1350    continue
         do 1370 k=1,i-1
           do 1360 ij=1,ncum
             if(i.le.inb(ij))then
               awat=elij(ij,k,i)-(1.-ep(ij,i))*clw(ij,i)
               awat=max(awat,0.0)
               fq(ij,i)=fq(ij,i) &
               +g*dpinv*ment(ij,k,i)*(qent(ij,k,i)-awat-q(ij,i))
               fu(ij,i)=fu(ij,i) &
               +g*dpinv*ment(ij,k,i)*(uent(ij,k,i)-u(ij,i))
               fv(ij,i)=fv(ij,i) &
               +g*dpinv*ment(ij,k,i)*(vent(ij,k,i)-v(ij,i))
! (saturated updrafts resulting from mixing)               ! cld
               qcond(ij,i)=qcond(ij,i)+(elij(ij,k,i)-awat) ! cld
               nqcond(ij,i)=nqcond(ij,i)+1.                ! cld
             endif
 1360      continue
 1370    continue
         do 1390 k=i,nl+1
           do 1380 ij=1,ncum
             if((i.le.inb(ij)).and.(k.le.inb(ij)))then
               fq(ij,i)=fq(ij,i) &
                        +g*dpinv*ment(ij,k,i)*(qent(ij,k,i)-q(ij,i))
               fu(ij,i)=fu(ij,i) &
                        +g*dpinv*ment(ij,k,i)*(uent(ij,k,i)-u(ij,i))
               fv(ij,i)=fv(ij,i) &
                        +g*dpinv*ment(ij,k,i)*(vent(ij,k,i)-v(ij,i))
             endif
 1380      continue
 1390    continue
          do 1400 ij=1,ncum
           if(i.le.inb(ij))then
             fq(ij,i)=fq(ij,i) &
                      +sigd*evap(ij,i)+g*(mp(ij,i+1)* &
                      (qp(ij,i+1)-q(ij,i)) &
                      -mp(ij,i)*(qp(ij,i)-q(ij,i-1)))*dpinv
             fu(ij,i)=fu(ij,i) &
                      +g*(mp(ij,i+1)*(up(ij,i+1)-u(ij,i))-mp(ij,i)* &
                      (up(ij,i)-u(ij,i-1)))*dpinv
             fv(ij,i)=fv(ij,i) &
                     +g*(mp(ij,i+1)*(vp(ij,i+1)-v(ij,i))-mp(ij,i)* &
                     (vp(ij,i)-v(ij,i-1)))*dpinv
! (saturated downdrafts resulting from mixing)               ! cld
            do k=i+1,inb(ij)                                 ! cld
             qcond(ij,i)=qcond(ij,i)+elij(ij,k,i)            ! cld
             nqcond(ij,i)=nqcond(ij,i)+1.                    ! cld
            enddo                                            ! cld
! (particular case: no detraining level is found)            ! cld
            if (nent(ij,i).eq.0) then                        ! cld
             qcond(ij,i)=qcond(ij,i)+(1.-ep(ij,i))*clw(ij,i) ! cld
             nqcond(ij,i)=nqcond(ij,i)+1.                    ! cld
            endif                                            ! cld
            if (nqcond(ij,i).ne.0.) then                     ! cld
             qcond(ij,i)=qcond(ij,i)/nqcond(ij,i)            ! cld
            endif                                            ! cld
           endif
 1400     continue
 1500   continue
!
!   *** Adjust tendencies at top of convection layer to reflect  ***
!   ***       actual position of the level zero cape             ***
!
        do 503 ij=1,ncum
        fqold=fq(ij,inb(ij))
        fq(ij,inb(ij))=fq(ij,inb(ij))*(1.-frac(ij))
        fq(ij,inb(ij)-1)=fq(ij,inb(ij)-1) &
         +frac(ij)*fqold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/ &
         (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))*lv(ij,inb(ij)) &
         /lv(ij,inb(ij)-1)
        ftold=ft(ij,inb(ij))
        ft(ij,inb(ij))=ft(ij,inb(ij))*(1.-frac(ij))
        ft(ij,inb(ij)-1)=ft(ij,inb(ij)-1) &
         +frac(ij)*ftold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/ &
         (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))*cpn(ij,inb(ij)) &
         /cpn(ij,inb(ij)-1)
        fuold=fu(ij,inb(ij))
        fu(ij,inb(ij))=fu(ij,inb(ij))*(1.-frac(ij))
        fu(ij,inb(ij)-1)=fu(ij,inb(ij)-1) &
         +frac(ij)*fuold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/ &
         (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
        fvold=fv(ij,inb(ij))
        fv(ij,inb(ij))=fv(ij,inb(ij))*(1.-frac(ij))
        fv(ij,inb(ij)-1)=fv(ij,inb(ij)-1) &
        +frac(ij)*fvold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/ &
         (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
 503    continue
!
!   ***   Very slightly adjust tendencies to force exact   ***
!   ***     enthalpy, momentum and tracer conservation     ***
!
        do 682 ij=1,ncum
        ents(ij)=0.0
        uav(ij)=0.0
        vav(ij)=0.0
        do 681 i=1,inb(ij)
         ents(ij)=ents(ij) &
        +(cpn(ij,i)*ft(ij,i)+lv(ij,i)*fq(ij,i))*(ph(ij,i)-ph(ij,i+1))
         uav(ij)=uav(ij)+fu(ij,i)*(ph(ij,i)-ph(ij,i+1))
         vav(ij)=vav(ij)+fv(ij,i)*(ph(ij,i)-ph(ij,i+1))
  681    continue
  682   continue
        do 683 ij=1,ncum
        ents(ij)=ents(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
        uav(ij)=uav(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
        vav(ij)=vav(ij)/(ph(ij,1)-ph(ij,inb(ij)+1))
 683    continue
        do 642 ij=1,ncum
        do 641 i=1,inb(ij)
         ft(ij,i)=ft(ij,i)-ents(ij)/cpn(ij,i)
         fu(ij,i)=(1.-cu)*(fu(ij,i)-uav(ij))
         fv(ij,i)=(1.-cu)*(fv(ij,i)-vav(ij))
  641    continue
 642    continue
!
        do 1810 k=1,nl+1
          do 1800 i=1,ncum
            if((q(i,k)+delt*fq(i,k)).lt.0.0)iflag(i)=10
 1800     continue
 1810   continue
!
!
        do 1900 i=1,ncum
          if(iflag(i).gt.2)then
          precip(i)=0.0
          cbmf(i)=0.0
          endif
 1900   continue
        do 1920 k=1,nl
         do 1910 i=1,ncum
           if(iflag(i).gt.2)then
             ft(i,k)=0.0
             fq(i,k)=0.0
             fu(i,k)=0.0
             fv(i,k)=0.0
             qcondc(i,k)=0.0                               ! cld
           endif
 1910    continue
 1920   continue

        do k=1,nl+1
        do i=1,ncum
          Ma(i,k) = 0.
        enddo
        enddo
        do k=nl,1,-1
        do i=1,ncum
          Ma(i,k) = Ma(i,k+1)+m(i,k)
        enddo
        enddo

!
!   *** diagnose the in-cloud mixing ratio   ***
!   ***           of condensed water         ***
!
      DO ij=1,ncum
       do i=1,nd
        mac(ij,i)=0.0
        wa(ij,i)=0.0
        siga(ij,i)=0.0
       enddo
       do i=nk(ij),inb(ij)
       do k=i+1,inb(ij)+1
        mac(ij,i)=mac(ij,i)+m(ij,k)
       enddo
       enddo
       do i=icb(ij),inb(ij)-1
        ax(ij,i)=0.
        do j=icb(ij),i
         ax(ij,i)=ax(ij,i)+rrd*(tvp(ij,j)-tv(ij,j))             &
             *(ph(ij,j)-ph(ij,j+1))/p(ij,j)
        enddo
        if (ax(ij,i).gt.0.0) then
         wa(ij,i)=sqrt(2.*ax(ij,i))
        endif
       enddo
       do i=1,nl
        if (wa(ij,i).gt.0.0)                                 &
          siga(ij,i)=mac(ij,i)/wa(ij,i)                         &
              *rrd*tvp(ij,i)/p(ij,i)/100./delta
        siga(ij,i) = min(siga(ij,i),1.0)
        qcondc(ij,i)=siga(ij,i)*clw(ij,i)*(1.-ep(ij,i))         &
                + (1.-siga(ij,i))*qcond(ij,i)
       enddo
      ENDDO

        return
        end
