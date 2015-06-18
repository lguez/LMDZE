
      SUBROUTINE cv3_undilute2(nloc,ncum,nd,icb,icbs,nk &
                             ,tnk,qnk,gznk,t,qs,gz &
                             ,p,h,tv,lv,pbase,buoybase,plcl &
                             ,inb,tp,tvp,clw,hp,ep,sigp,buoy)
      use conema3_m
            use cv3_param_m
            use cvthermo
      implicit none

!---------------------------------------------------------------------
! Purpose:
!     FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
!     &
!     COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
!     FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
!     &
!     FIND THE LEVEL OF NEUTRAL BUOYANCY
!
! Main differences convect3/convect4:
!     - icbs (input) is the first level above LCL (may differ from icb)
!     - many minor differences in the iterations
!     - condensed water not removed from tvp in convect3
!   - vertical profile of buoyancy computed here (use of buoybase)
!   - the determination of inb is different
!   - no inb1, only inb in output
!---------------------------------------------------------------------


! inputs:
      integer, intent(in):: ncum, nd, nloc
      integer icb(nloc), icbs(nloc), nk(nloc)
      real t(nloc,nd), qs(nloc,nd), gz(nloc,nd)
      real p(nloc,nd)
      real tnk(nloc), qnk(nloc), gznk(nloc)
      real lv(nloc,nd), tv(nloc,nd), h(nloc,nd)
      real pbase(nloc), buoybase(nloc), plcl(nloc)

! outputs:
      integer inb(nloc)
      real tp(nloc,nd), tvp(nloc,nd), clw(nloc,nd)
      real ep(nloc,nd), sigp(nloc,nd), hp(nloc,nd)
      real buoy(nloc,nd)

! local variables:
      integer i, k
      real tg,qg,ahg,alv,s,tc,es,denom,rg
      real pden
      real ah0(nloc)

!=====================================================================
! --- SOME INITIALIZATIONS
!=====================================================================

      do 170 k=1,nl
      do 160 i=1,ncum
       ep(i,k)=0.0
       sigp(i,k)=spfac
 160  continue
 170  continue

!=====================================================================
! --- FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
!=====================================================================
!
! ---       The procedure is to solve the equation.
!              cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
!
!   ***  Calculate certain parcel quantities, including static energy   ***
!
!
      do 240 i=1,ncum
         ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) &
               +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
 240  continue
!
!
!    ***  Find lifted parcel quantities above cloud base    ***
!
!
      do 300 k=minorig+1,nl
        do 290 i=1,ncum
! ori        if(k.ge.(icb(i)+1))then
          if(k.ge.(icbs(i)+1))then ! convect3
            tg=t(i,k)
            qg=qs(i,k)
!debug         alv=lv0-clmcpv*(t(i,k)-t0)
            alv=lv0-clmcpv*(t(i,k)-273.15)
!
! First iteration.
!
! ori           s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
           s=cpd*(1.-qnk(i))+cl*qnk(i) &
            +alv*alv*qg/(rrv*t(i,k)*t(i,k)) ! convect3
             s=1./s
! ori           ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
           ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gz(i,k) ! convect3
             tg=tg+s*(ah0(i)-ahg)
! ori           tg=max(tg,35.0)
!debug          tc=tg-t0
             tc=tg-273.15
             denom=243.5+tc
           denom=MAX(denom,1.0) ! convect3
! ori           if(tc.ge.0.0)then
            es=6.112*exp(17.67*tc/denom)
! ori           else
! ori          es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori           endif
            qg=eps*es/(p(i,k)-es*(1.-eps))
!
! Second iteration.
!
! ori           s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
! ori           s=1./s
! ori           ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
           ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gz(i,k) ! convect3
             tg=tg+s*(ah0(i)-ahg)
! ori           tg=max(tg,35.0)
!debug          tc=tg-t0
             tc=tg-273.15
             denom=243.5+tc
           denom=MAX(denom,1.0) ! convect3
! ori           if(tc.ge.0.0)then
            es=6.112*exp(17.67*tc/denom)
! ori           else
! ori          es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori           endif
            qg=eps*es/(p(i,k)-es*(1.-eps))
!
!debug          alv=lv0-clmcpv*(t(i,k)-t0)
             alv=lv0-clmcpv*(t(i,k)-273.15)
!      print*,'cpd dans convect2 ',cpd
!      print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
!      print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd

! ori c approximation here:
! ori        tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd

! convect3: no approximation:
           tp(i,k)=(ah0(i)-gz(i,k)-alv*qg)/(cpd+(cl-cpd)*qnk(i))

               clw(i,k)=qnk(i)-qg
               clw(i,k)=max(0.0,clw(i,k))
               rg=qg/(1.-qnk(i))
! ori               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
! convect3: (qg utilise au lieu du vrai mixing ratio rg):
               tvp(i,k)=tp(i,k)*(1.+qg/eps-qnk(i)) ! whole thing
            endif
  290     continue
  300   continue
!
!=====================================================================
! --- SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF
! --- PRECIPITATION FALLING OUTSIDE OF CLOUD
! --- THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
!=====================================================================
!
! ori      do 320 k=minorig+1,nl
      do 320 k=1,nl ! convect3
        do 310 i=1,ncum
           pden=ptcrit-pbcrit
           ep(i,k)=(plcl(i)-p(i,k)-pbcrit)/pden*epmax
           ep(i,k)=amax1(ep(i,k),0.0)
           ep(i,k)=amin1(ep(i,k),epmax)
           sigp(i,k)=spfac
! ori          if(k.ge.(nk(i)+1))then
! ori            tca=tp(i,k)-t0
! ori            if(tca.ge.0.0)then
! ori              elacrit=elcrit
! ori            else
! ori              elacrit=elcrit*(1.0-tca/tlcrit)
! ori            endif
! ori            elacrit=max(elacrit,0.0)
! ori            ep(i,k)=1.0-elacrit/max(clw(i,k),1.0e-8)
! ori            ep(i,k)=max(ep(i,k),0.0 )
! ori            ep(i,k)=min(ep(i,k),1.0 )
! ori            sigp(i,k)=sigs
! ori          endif
 310    continue
 320  continue
!
!=====================================================================
! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
! --- VIRTUAL TEMPERATURE
!=====================================================================
!
! dans convect3, tvp est calcule en une seule fois, et sans retirer
! l'eau condensee (~> reversible CAPE)
!
! ori      do 340 k=minorig+1,nl
! ori        do 330 i=1,ncum
! ori        if(k.ge.(icb(i)+1))then
! ori          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
! oric         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
! oric         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
! ori        endif
! ori 330    continue
! ori 340  continue

! ori      do 350 i=1,ncum
! ori       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
! ori 350  continue

      do 350 i=1,ncum       ! convect3
       tp(i,nlp)=tp(i,nl)   ! convect3
 350  continue              ! convect3
!
!=====================================================================
!  --- EFFECTIVE VERTICAL PROFILE OF BUOYANCY (convect3 only):
!=====================================================================

!-- this is for convect3 only:

! first estimate of buoyancy:

      do 500 i=1,ncum
       do 501 k=1,nl
        buoy(i,k)=tvp(i,k)-tv(i,k)
 501   continue
 500  continue

! set buoyancy=buoybase for all levels below base
! for safety, set buoy(icb)=buoybase

      do 505 i=1,ncum
       do 506 k=1,nl
        if((k.ge.icb(i)).and.(k.le.nl).and.(p(i,k).ge.pbase(i)))then
         buoy(i,k)=buoybase(i)
        endif
 506   continue
       buoy(icb(i),k)=buoybase(i)
 505  continue

!-- end convect3

!=====================================================================
!  --- FIND THE FIRST MODEL LEVEL (INB) ABOVE THE PARCEL'S
!  --- LEVEL OF NEUTRAL BUOYANCY
!=====================================================================
!
!-- this is for convect3 only:

      do 510 i=1,ncum
       inb(i)=nl-1
 510  continue

      do 530 i=1,ncum
       do 535 k=1,nl-1
        if ((k.ge.icb(i)).and.(buoy(i,k).lt.dtovsh)) then
         inb(i)=MIN(inb(i),k)
        endif
 535   continue
 530  continue

!-- end convect3

! ori      do 510 i=1,ncum
! ori        cape(i)=0.0
! ori        capem(i)=0.0
! ori        inb(i)=icb(i)+1
! ori        inb1(i)=inb(i)
! ori 510  continue
!
! Originial Code
!
!     do 530 k=minorig+1,nl-1
!       do 520 i=1,ncum
!         if(k.ge.(icb(i)+1))then
!           by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
!           byp=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
!           cape(i)=cape(i)+by
!           if(by.ge.0.0)inb1(i)=k+1
!           if(cape(i).gt.0.0)then
!             inb(i)=k+1
!             capem(i)=cape(i)
!           endif
!         endif
!520    continue
!530  continue
!     do 540 i=1,ncum
!         byp=(tvp(i,nl)-tv(i,nl))*dph(i,nl)/p(i,nl)
!         cape(i)=capem(i)+byp
!         defrac=capem(i)-cape(i)
!         defrac=max(defrac,0.001)
!         frac(i)=-cape(i)/defrac
!         frac(i)=min(frac(i),1.0)
!         frac(i)=max(frac(i),0.0)
!540   continue
!
! K Emanuel fix
!
!     call zilch(byp,ncum)
!     do 530 k=minorig+1,nl-1
!       do 520 i=1,ncum
!         if(k.ge.(icb(i)+1))then
!           by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
!           cape(i)=cape(i)+by
!           if(by.ge.0.0)inb1(i)=k+1
!           if(cape(i).gt.0.0)then
!             inb(i)=k+1
!             capem(i)=cape(i)
!             byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
!           endif
!         endif
!520    continue
!530  continue
!     do 540 i=1,ncum
!         inb(i)=max(inb(i),inb1(i))
!         cape(i)=capem(i)+byp(i)
!         defrac=capem(i)-cape(i)
!         defrac=max(defrac,0.001)
!         frac(i)=-cape(i)/defrac
!         frac(i)=min(frac(i),1.0)
!         frac(i)=max(frac(i),0.0)
!540   continue
!
! J Teixeira fix
!
! ori      call zilch(byp,ncum)
! ori      do 515 i=1,ncum
! ori        lcape(i)=.true.
! ori 515  continue
! ori      do 530 k=minorig+1,nl-1
! ori        do 520 i=1,ncum
! ori          if(cape(i).lt.0.0)lcape(i)=.false.
! ori          if((k.ge.(icb(i)+1)).and.lcape(i))then
! ori            by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
! ori            byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
! ori            cape(i)=cape(i)+by
! ori            if(by.ge.0.0)inb1(i)=k+1
! ori            if(cape(i).gt.0.0)then
! ori              inb(i)=k+1
! ori              capem(i)=cape(i)
! ori            endif
! ori          endif
! ori 520    continue
! ori 530  continue
! ori      do 540 i=1,ncum
! ori          cape(i)=capem(i)+byp(i)
! ori          defrac=capem(i)-cape(i)
! ori          defrac=max(defrac,0.001)
! ori          frac(i)=-cape(i)/defrac
! ori          frac(i)=min(frac(i),1.0)
! ori          frac(i)=max(frac(i),0.0)
! ori 540  continue
!
!=====================================================================
! ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
!=====================================================================
!
!ym      do i=1,ncum*nlp
!ym       hp(i,1)=h(i,1)
!ym      enddo

      do k=1,nlp
        do i=1,ncum
        hp(i,k)=h(i,k)
      enddo
      enddo

      do 600 k=minorig+1,nl
        do 590 i=1,ncum
        if((k.ge.icb(i)).and.(k.le.inb(i)))then
          hp(i,k)=h(i,nk(i))+(lv(i,k)+(cpd-cpv)*t(i,k))*ep(i,k)*clw(i,k)
        endif
 590    continue
 600  continue

        return
        end
