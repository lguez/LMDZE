
      SUBROUTINE cv_undilute2(nloc,ncum,nd,icb,nk &
                             ,tnk,qnk,gznk,t,qs,gz &
                             ,p,dph,h,tv,lv &
                             ,inb,inb1,tp,tvp,clw,hp,ep,sigp,frac)
            use cvthermo
            use cv_param
      implicit none

!---------------------------------------------------------------------
! Purpose:
!     FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
!     &
!     COMPUTE THE PRECIPITATION EFFICIENCIES AND THE
!     FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
!     &
!     FIND THE LEVEL OF NEUTRAL BUOYANCY
!---------------------------------------------------------------------


! inputs:
      integer, intent(in):: ncum, nd, nloc
      integer icb(nloc), nk(nloc)
      real t(nloc,nd), qs(nloc,nd), gz(nloc,nd)
      real p(nloc,nd), dph(nloc,nd)
      real tnk(nloc), qnk(nloc), gznk(nloc)
      real lv(nloc,nd), tv(nloc,nd), h(nloc,nd)

! outputs:
      integer inb(nloc), inb1(nloc)
      real tp(nloc,nd), tvp(nloc,nd), clw(nloc,nd)
      real ep(nloc,nd), sigp(nloc,nd), hp(nloc,nd)
      real frac(nloc)

! local variables:
      integer i, k
      real tg,qg,ahg,alv,s,tc,es,denom,rg,tca,elacrit
      real by, defrac
      real ah0(nloc), cape(nloc), capem(nloc), byp(nloc)
      logical lcape(nloc)

!=====================================================================
! --- SOME INITIALIZATIONS
!=====================================================================

      do 170 k=1,nl
      do 160 i=1,ncum
       ep(i,k)=0.0
       sigp(i,k)=sigs
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
               +qnk(i)*(lv0-clmcpv*(tnk(i)-t0))+gznk(i)
 240  continue
!
!
!    ***  Find lifted parcel quantities above cloud base    ***
!
!
      do 300 k=minorig+1,nl
        do 290 i=1,ncum
          if(k.ge.(icb(i)+1))then
            tg=t(i,k)
            qg=qs(i,k)
            alv=lv0-clmcpv*(t(i,k)-t0)
!
! First iteration.
!
             s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
             s=1./s
             ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
             tg=tg+s*(ah0(i)-ahg)
             tg=max(tg,35.0)
             tc=tg-t0
             denom=243.5+tc
             if(tc.ge.0.0)then
            es=6.112*exp(17.67*tc/denom)
             else
            es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
             endif
            qg=eps*es/(p(i,k)-es*(1.-eps))
!
! Second iteration.
!
             s=cpd+alv*alv*qg/(rrv*t(i,k)*t(i,k))
             s=1./s
             ahg=cpd*tg+(cl-cpd)*qnk(i)*t(i,k)+alv*qg+gz(i,k)
             tg=tg+s*(ah0(i)-ahg)
             tg=max(tg,35.0)
             tc=tg-t0
             denom=243.5+tc
             if(tc.ge.0.0)then
            es=6.112*exp(17.67*tc/denom)
             else
            es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
             endif
            qg=eps*es/(p(i,k)-es*(1.-eps))
!
             alv=lv0-clmcpv*(t(i,k)-t0)
!      print*,'cpd dans convect2 ',cpd
!      print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
!      print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd
        tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd
!              if (.not.cpd.gt.1000.) then
!                  print*,'CPD=',cpd
!                  stop
!              endif
               clw(i,k)=qnk(i)-qg
               clw(i,k)=max(0.0,clw(i,k))
               rg=qg/(1.-qnk(i))
               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
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
      do 320 k=minorig+1,nl
        do 310 i=1,ncum
          if(k.ge.(nk(i)+1))then
            tca=tp(i,k)-t0
            if(tca.ge.0.0)then
              elacrit=elcrit
            else
              elacrit=elcrit*(1.0-tca/tlcrit)
            endif
            elacrit=max(elacrit,0.0)
            ep(i,k)=1.0-elacrit/max(clw(i,k),1.0e-8)
            ep(i,k)=max(ep(i,k),0.0 )
            ep(i,k)=min(ep(i,k),1.0 )
            sigp(i,k)=sigs
          endif
 310    continue
 320  continue
!
!=====================================================================
! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
! --- VIRTUAL TEMPERATURE
!=====================================================================
!
      do 340 k=minorig+1,nl
        do 330 i=1,ncum
        if(k.ge.(icb(i)+1))then
          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
!         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
!         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
        endif
 330    continue
 340  continue
      do 350 i=1,ncum
       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
 350  continue
!
!=====================================================================
!  --- FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S
!  --- HIGHEST LEVEL OF NEUTRAL BUOYANCY
!  --- AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)
!=====================================================================
!
      do 510 i=1,ncum
        cape(i)=0.0
        capem(i)=0.0
        inb(i)=icb(i)+1
        inb1(i)=inb(i)
 510  continue
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
      call zilch(byp,ncum)
      do 515 i=1,ncum
        lcape(i)=.true.
 515  continue
      do 530 k=minorig+1,nl-1
        do 520 i=1,ncum
          if(cape(i).lt.0.0)lcape(i)=.false.
          if((k.ge.(icb(i)+1)).and.lcape(i))then
            by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
            byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
            cape(i)=cape(i)+by
            if(by.ge.0.0)inb1(i)=k+1
            if(cape(i).gt.0.0)then
              inb(i)=k+1
              capem(i)=cape(i)
            endif
          endif
 520    continue
 530  continue
      do 540 i=1,ncum
          cape(i)=capem(i)+byp(i)
          defrac=capem(i)-cape(i)
          defrac=max(defrac,0.001)
          frac(i)=-cape(i)/defrac
          frac(i)=min(frac(i),1.0)
          frac(i)=max(frac(i),0.0)
 540  continue
!
!=====================================================================
! ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
!=====================================================================
!
! initialization:
      do i=1,ncum*nlp
       hp(i,1)=h(i,1)
      enddo

      do 600 k=minorig+1,nl
        do 590 i=1,ncum
        if((k.ge.icb(i)).and.(k.le.inb(i)))then
          hp(i,k)=h(i,nk(i))+(lv(i,k)+(cpd-cpv)*t(i,k))*ep(i,k)*clw(i,k)
        endif
 590    continue
 600  continue
!
        return
        end
