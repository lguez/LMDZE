!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/cv_routines.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      SUBROUTINE cv_param(nd)
      implicit none

c------------------------------------------------------------
c Set parameters for convectL
c (includes microphysical parameters and parameters that 
c  control the rate of approach to quasi-equilibrium) 
c------------------------------------------------------------

C   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
C   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
C   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
C   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
C   ***               BETWEEN 0 C AND TLCRIT)                        ***
C   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
C   ***                       FORMULATION                            ***
C   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
C   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
C   ***                        OF CLOUD                              ***
C   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
C   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
C   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF RAIN                             ***
C   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF SNOW                             ***
C   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
C   ***                         TRANSPORT                            ***
C   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
C   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
C   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
C   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
C   ***                   (DAMP MUST BE LESS THAN 1)                 ***

      include "cvparam.h"
      integer nd

c noff: integer limit for convection (nd-noff)
c minorig: First level of convection

      noff = 2
      minorig = 2

      nl=nd-noff
      nlp=nl+1
      nlm=nl-1

      elcrit=0.0011
      tlcrit=-55.0
      entp=1.5
      sigs=0.12
      sigd=0.05
      omtrain=50.0
      omtsnow=5.5
      coeffr=1.0
      coeffs=0.8
      dtmax=0.9
c
      cu=0.70
c
      betad=10.0
c
      damp=0.1
      alpha=0.2
c
      delta=0.01  ! cld
c
      return
      end

      SUBROUTINE cv_prelim(len,nd,ndp1,t,q,p,ph
     :                    ,lv,cpn,tv,gz,h,hm)
            use cvthermo
      implicit none

!=====================================================================
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
!=====================================================================

c inputs:
      integer len, nd, ndp1
      real t(len,nd), q(len,nd), p(len,nd), ph(len,ndp1)

c outputs:
      real lv(len,nd), cpn(len,nd), tv(len,nd)
      real gz(len,nd), h(len,nd), hm(len,nd)

c local variables:
      integer k, i
      real cpx(len,nd)

      include "cvparam.h"


      do 110 k=1,nlp
        do 100 i=1,len
          lv(i,k)= lv0-clmcpv*(t(i,k)-t0)
          cpn(i,k)=cpd*(1.0-q(i,k))+cpv*q(i,k)
          cpx(i,k)=cpd*(1.0-q(i,k))+cl*q(i,k)
          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
 100    continue
 110  continue
c
c gz = phi at the full levels (same as p).
c
      do 120 i=1,len
        gz(i,1)=0.0
 120  continue
      do 140 k=2,nlp
        do 130 i=1,len
          gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k))
     &         *(p(i,k-1)-p(i,k))/ph(i,k)
 130    continue
 140  continue
c
c h  = phi + cpT (dry static energy).
c hm = phi + cp(T-Tbase)+Lq
c
      do 170 k=1,nlp
        do 160 i=1,len
          h(i,k)=gz(i,k)+cpn(i,k)*t(i,k)
          hm(i,k)=gz(i,k)+cpx(i,k)*(t(i,k)-t(i,1))+lv(i,k)*q(i,k)
 160    continue
 170  continue

      return
      end

      SUBROUTINE cv_feed(len,nd,t,q,qs,p,hm,gz
     :                  ,nk,icb,icbmax,iflag,tnk,qnk,gznk,plcl)
      implicit none

C================================================================
C Purpose: CONVECTIVE FEED
C================================================================

      include "cvparam.h"

c inputs:
	  integer len, nd
      real t(len,nd), q(len,nd), qs(len,nd), p(len,nd)
      real hm(len,nd), gz(len,nd)

c outputs:
	  integer iflag(len), nk(len), icb(len), icbmax
      real tnk(len), qnk(len), gznk(len), plcl(len)

c local variables:
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
         if((hm(i,k).lt.work(i)).and.
     &      (hm(i,k).lt.hm(i,k-1)))then
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
c
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
       if(((t(i,nk(i)).lt.250.0).or.
     &      (q(i,nk(i)).le.0.0).or.
     &      (p(i,ihmin(i)).lt.400.0)).and.
     &      (iflag(i).eq.0))iflag(i)=7
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
c
        rh(i)=qnk(i)/qsnk(i)
        rh(i)=min(1.0,rh(i))
        chi(i)=tnk(i)/(1669.0-122.0*rh(i)-tnk(i))
        plcl(i)=pnk(i)*(rh(i)**chi(i))
        if(((plcl(i).lt.200.0).or.(plcl(i).ge.2000.0))
     &   .and.(iflag(i).eq.0))iflag(i)=8
 260   continue
!-------------------------------------------------------------------
! --- Calculate first level above lcl (=icb)
!-------------------------------------------------------------------
      do 270 i=1,len
       icb(i)=nlm
 270  continue
c
      do 290 k=minorig,nl
        do 280 i=1,len
          if((k.ge.(nk(i)+1)).and.(p(i,k).lt.plcl(i)))
     &    icb(i)=min(icb(i),k)
 280    continue
 290  continue
c
      do 300 i=1,len
        if((icb(i).ge.nlm).and.(iflag(i).eq.0))iflag(i)=9
 300  continue
c
c Compute icbmax.
c
      icbmax=2
      do 310 i=1,len
        icbmax=max(icbmax,icb(i))
 310  continue

      return
      end

      SUBROUTINE cv_undilute1(len,nd,t,q,qs,gz,p,nk,icb,icbmax
     :                       ,tp,tvp,clw)
            use cvthermo
      implicit none

      include "cvparam.h"

c inputs:
      integer len, nd
      integer nk(len), icb(len), icbmax
      real t(len,nd), q(len,nd), qs(len,nd), gz(len,nd)
      real p(len,nd)

c outputs:
      real tp(len,nd), tvp(len,nd), clw(len,nd)

c local variables:
      integer i, k
      real tg, qg, alv, s, ahg, tc, denom, es, rg
      real ah0(len), cpp(len)
      real tnk(len), qnk(len), gznk(len), ticb(len), gzicb(len)

!-------------------------------------------------------------------
! --- Calculates the lifted parcel virtual temperature at nk,
! --- the actual temperature, and the adiabatic
! --- liquid water content. The procedure is to solve the equation.
!     cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
!-------------------------------------------------------------------

      do 320 i=1,len
        tnk(i)=t(i,nk(i))
        qnk(i)=q(i,nk(i))
        gznk(i)=gz(i,nk(i))
        ticb(i)=t(i,icb(i))
        gzicb(i)=gz(i,icb(i))
 320  continue
c
c   ***  Calculate certain parcel quantities, including static energy   ***
c
      do 330 i=1,len
        ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i)
     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
        cpp(i)=cpd*(1.-qnk(i))+qnk(i)*cpv
 330  continue
c
c   ***   Calculate lifted parcel quantities below cloud base   ***
c
        do 350 k=minorig,icbmax-1
          do 340 i=1,len
           tp(i,k)=tnk(i)-(gz(i,k)-gznk(i))/cpp(i)
           tvp(i,k)=tp(i,k)*(1.+qnk(i)*epsi)
  340     continue
  350   continue
c
c    ***  Find lifted parcel quantities above cloud base    ***
c
        do 360 i=1,len
         tg=ticb(i)
         qg=qs(i,icb(i))
         alv=lv0-clmcpv*(ticb(i)-t0)
c
c First iteration.
c
          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
          s=1./s
          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          tg=tg+s*(ah0(i)-ahg)
          tg=max(tg,35.0)
          tc=tg-t0
          denom=243.5+tc
          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
          else
           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
          endif
          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
c
c Second iteration.
c
          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
          s=1./s
          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          tg=tg+s*(ah0(i)-ahg)
          tg=max(tg,35.0)
          tc=tg-t0
          denom=243.5+tc
          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
          else
           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
          end if
          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
c
         alv=lv0-clmcpv*(ticb(i)-273.15)
         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
     &   -gz(i,icb(i))-alv*qg)/cpd
         clw(i,icb(i))=qnk(i)-qg
         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         rg=qg/(1.-qnk(i))
         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
  360   continue
c
      do 380 k=minorig,icbmax
       do 370 i=1,len
         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
 370   continue
 380  continue
c
      return
      end

      SUBROUTINE cv_trigger(len,nd,icb,cbmf,tv,tvp,iflag)
      implicit none

!-------------------------------------------------------------------
! --- Test for instability.
! --- If there was no convection at last time step and parcel
! --- is stable at icb, then set iflag to 4.
!-------------------------------------------------------------------
 
      include "cvparam.h"

c inputs:
       integer len, nd, icb(len)
       real cbmf(len), tv(len,nd), tvp(len,nd)

c outputs:
       integer iflag(len) ! also an input

c local variables:
       integer i


      do 390 i=1,len
        if((cbmf(i).eq.0.0) .and.(iflag(i).eq.0).and.
     &  (tvp(i,icb(i)).le.(tv(i,icb(i))-dtmax)))iflag(i)=4
 390  continue
 
      return
      end

      SUBROUTINE cv_compress( len,nloc,ncum,nd
     :   ,iflag1,nk1,icb1
     :   ,cbmf1,plcl1,tnk1,qnk1,gznk1
     :   ,t1,q1,qs1,u1,v1,gz1
     :   ,h1,lv1,cpn1,p1,ph1,tv1,tp1,tvp1,clw1
     o   ,iflag,nk,icb
     o   ,cbmf,plcl,tnk,qnk,gznk
     o   ,t,q,qs,u,v,gz,h,lv,cpn,p,ph,tv,tp,tvp,clw 
     o   ,dph          )
      implicit none

      include "cvparam.h"

c inputs:
      integer len,ncum,nd,nloc
      integer iflag1(len),nk1(len),icb1(len)
      real cbmf1(len),plcl1(len),tnk1(len),qnk1(len),gznk1(len)
      real t1(len,nd),q1(len,nd),qs1(len,nd),u1(len,nd),v1(len,nd)
      real gz1(len,nd),h1(len,nd),lv1(len,nd),cpn1(len,nd)
      real p1(len,nd),ph1(len,nd+1),tv1(len,nd),tp1(len,nd)
      real tvp1(len,nd),clw1(len,nd)

c outputs:
      integer iflag(nloc),nk(nloc),icb(nloc)
      real cbmf(nloc),plcl(nloc),tnk(nloc),qnk(nloc),gznk(nloc)
      real t(nloc,nd),q(nloc,nd),qs(nloc,nd),u(nloc,nd),v(nloc,nd)
      real gz(nloc,nd),h(nloc,nd),lv(nloc,nd),cpn(nloc,nd)
      real p(nloc,nd),ph(nloc,nd+1),tv(nloc,nd),tp(nloc,nd)
      real tvp(nloc,nd),clw(nloc,nd)
      real dph(nloc,nd)

c local variables:
      integer i,k,nn


      do 110 k=1,nl+1
       nn=0
      do 100 i=1,len
      if(iflag1(i).eq.0)then
        nn=nn+1
        t(nn,k)=t1(i,k)
        q(nn,k)=q1(i,k)
        qs(nn,k)=qs1(i,k)
        u(nn,k)=u1(i,k)
        v(nn,k)=v1(i,k)
        gz(nn,k)=gz1(i,k)
        h(nn,k)=h1(i,k)
        lv(nn,k)=lv1(i,k)
        cpn(nn,k)=cpn1(i,k)
        p(nn,k)=p1(i,k)
        ph(nn,k)=ph1(i,k)
        tv(nn,k)=tv1(i,k)
        tp(nn,k)=tp1(i,k)
        tvp(nn,k)=tvp1(i,k)
        clw(nn,k)=clw1(i,k)
      endif
 100    continue
 110  continue

      if (nn.ne.ncum) then
         print*,'strange! nn not equal to ncum: ',nn,ncum
         stop
      endif

      nn=0
      do 150 i=1,len
      if(iflag1(i).eq.0)then
      nn=nn+1
      cbmf(nn)=cbmf1(i)
      plcl(nn)=plcl1(i)
      tnk(nn)=tnk1(i)
      qnk(nn)=qnk1(i)
      gznk(nn)=gznk1(i)
      nk(nn)=nk1(i)
      icb(nn)=icb1(i)
      iflag(nn)=iflag1(i)
      endif
 150  continue

      do 170 k=1,nl
       do 160 i=1,ncum
        dph(i,k)=ph(i,k)-ph(i,k+1)
 160   continue
 170  continue

      return
      end

      SUBROUTINE cv_undilute2(nloc,ncum,nd,icb,nk
     :                       ,tnk,qnk,gznk,t,q,qs,gz
     :                       ,p,dph,h,tv,lv
     o                       ,inb,inb1,tp,tvp,clw,hp,ep,sigp,frac)
            use cvthermo
      implicit none

C---------------------------------------------------------------------
C Purpose:
C     FIND THE REST OF THE LIFTED PARCEL TEMPERATURES
C     &
C     COMPUTE THE PRECIPITATION EFFICIENCIES AND THE 
C     FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
C     &
C     FIND THE LEVEL OF NEUTRAL BUOYANCY
C---------------------------------------------------------------------

      include "cvparam.h"

c inputs:
      integer ncum, nd, nloc
      integer icb(nloc), nk(nloc)
      real t(nloc,nd), q(nloc,nd), qs(nloc,nd), gz(nloc,nd)
      real p(nloc,nd), dph(nloc,nd)
      real tnk(nloc), qnk(nloc), gznk(nloc)
      real lv(nloc,nd), tv(nloc,nd), h(nloc,nd)

c outputs:
      integer inb(nloc), inb1(nloc)
      real tp(nloc,nd), tvp(nloc,nd), clw(nloc,nd)
      real ep(nloc,nd), sigp(nloc,nd), hp(nloc,nd)
      real frac(nloc)

c local variables:
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
c
c ---       The procedure is to solve the equation.
c              cp*tp+L*qp+phi=cp*tnk+L*qnk+gznk.
c
c   ***  Calculate certain parcel quantities, including static energy   ***
c
c
      do 240 i=1,ncum
         ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i)
     &         +qnk(i)*(lv0-clmcpv*(tnk(i)-t0))+gznk(i)
 240  continue
c
c
c    ***  Find lifted parcel quantities above cloud base    ***
c
c
	do 300 k=minorig+1,nl
	  do 290 i=1,ncum
	    if(k.ge.(icb(i)+1))then
	      tg=t(i,k)
	      qg=qs(i,k)
	      alv=lv0-clmcpv*(t(i,k)-t0)
c
c First iteration.
c
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
c
c Second iteration.
c
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
c
	       alv=lv0-clmcpv*(t(i,k)-t0)
c      print*,'cpd dans convect2 ',cpd
c      print*,'tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd'
c      print*,tp(i,k),ah0(i),cl,cpd,qnk(i),t(i,k),gz(i,k),alv,qg,cpd
        tp(i,k)=(ah0(i)-(cl-cpd)*qnk(i)*t(i,k)-gz(i,k)-alv*qg)/cpd
c              if (.not.cpd.gt.1000.) then
c                  print*,'CPD=',cpd
c                  stop
c              endif
               clw(i,k)=qnk(i)-qg
               clw(i,k)=max(0.0,clw(i,k))
               rg=qg/(1.-qnk(i))
               tvp(i,k)=tp(i,k)*(1.+rg*epsi)
            endif
  290     continue
  300   continue
c
!=====================================================================
! --- SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF
! --- PRECIPITATION FALLING OUTSIDE OF CLOUD
! --- THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)
!=====================================================================
c
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
c
!=====================================================================
! --- CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL
! --- VIRTUAL TEMPERATURE
!=====================================================================
c
      do 340 k=minorig+1,nl
        do 330 i=1,ncum
        if(k.ge.(icb(i)+1))then
          tvp(i,k)=tvp(i,k)*(1.0-qnk(i)+ep(i,k)*clw(i,k))
c         print*,'i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)'
c         print*, i,k,tvp(i,k),qnk(i),ep(i,k),clw(i,k)
        endif
 330    continue
 340  continue
      do 350 i=1,ncum
       tvp(i,nlp)=tvp(i,nl)-(gz(i,nlp)-gz(i,nl))/cpd
 350  continue
c
c=====================================================================
c  --- FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S
c  --- HIGHEST LEVEL OF NEUTRAL BUOYANCY
c  --- AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)
c=====================================================================
c
      do 510 i=1,ncum
        cape(i)=0.0
        capem(i)=0.0
        inb(i)=icb(i)+1
        inb1(i)=inb(i)
 510  continue
c
c Originial Code
c
c     do 530 k=minorig+1,nl-1
c       do 520 i=1,ncum
c         if(k.ge.(icb(i)+1))then
c           by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
c           byp=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
c           cape(i)=cape(i)+by
c           if(by.ge.0.0)inb1(i)=k+1
c           if(cape(i).gt.0.0)then
c             inb(i)=k+1
c             capem(i)=cape(i)
c           endif
c         endif
c520    continue
c530  continue
c     do 540 i=1,ncum
c         byp=(tvp(i,nl)-tv(i,nl))*dph(i,nl)/p(i,nl)
c         cape(i)=capem(i)+byp
c         defrac=capem(i)-cape(i)
c         defrac=max(defrac,0.001)
c         frac(i)=-cape(i)/defrac
c         frac(i)=min(frac(i),1.0)
c         frac(i)=max(frac(i),0.0)
c540   continue
c
c K Emanuel fix
c
c     call zilch(byp,ncum)
c     do 530 k=minorig+1,nl-1
c       do 520 i=1,ncum
c         if(k.ge.(icb(i)+1))then
c           by=(tvp(i,k)-tv(i,k))*dph(i,k)/p(i,k)
c           cape(i)=cape(i)+by
c           if(by.ge.0.0)inb1(i)=k+1
c           if(cape(i).gt.0.0)then
c             inb(i)=k+1
c             capem(i)=cape(i)
c             byp(i)=(tvp(i,k+1)-tv(i,k+1))*dph(i,k+1)/p(i,k+1)
c           endif
c         endif
c520    continue
c530  continue
c     do 540 i=1,ncum
c         inb(i)=max(inb(i),inb1(i))
c         cape(i)=capem(i)+byp(i)
c         defrac=capem(i)-cape(i)
c         defrac=max(defrac,0.001)
c         frac(i)=-cape(i)/defrac
c         frac(i)=min(frac(i),1.0)
c         frac(i)=max(frac(i),0.0)
c540   continue
c
c J Teixeira fix
c
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
c
c=====================================================================
c ---   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL
c=====================================================================
c
c initialization:
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
c
        return
        end
c
      SUBROUTINE cv_closure(nloc,ncum,nd,nk,icb
     :                     ,tv,tvp,p,ph,dph,plcl,cpn
     :                     ,iflag,cbmf)
            use cvthermo
      implicit none

c inputs:
      integer ncum, nd, nloc
      integer nk(nloc), icb(nloc)
      real tv(nloc,nd), tvp(nloc,nd), p(nloc,nd), dph(nloc,nd)
      real ph(nloc,nd+1) ! caution nd instead ndp1 to be consistent...
      real plcl(nloc), cpn(nloc,nd)

c outputs:
      integer iflag(nloc)
      real cbmf(nloc) ! also an input

c local variables:
      integer i, k, icbmax
      real dtpbl(nloc), dtmin(nloc), tvpplcl(nloc), tvaplcl(nloc)
      real work(nloc)

      include "cvparam.h"

c-------------------------------------------------------------------
c Compute icbmax. 
c-------------------------------------------------------------------

      icbmax=2
      do 230 i=1,ncum
       icbmax=max(icbmax,icb(i))
 230  continue

c=====================================================================
c ---  CALCULATE CLOUD BASE MASS FLUX 
c=====================================================================
c
c tvpplcl = parcel temperature lifted adiabatically from level
c           icb-1 to the LCL.
c tvaplcl = virtual temperature at the LCL.
c
      do 610 i=1,ncum
        dtpbl(i)=0.0
        tvpplcl(i)=tvp(i,icb(i)-1)
     &  -rrd*tvp(i,icb(i)-1)*(p(i,icb(i)-1)-plcl(i))
     &  /(cpn(i,icb(i)-1)*p(i,icb(i)-1))
        tvaplcl(i)=tv(i,icb(i))
     &  +(tvp(i,icb(i))-tvp(i,icb(i)+1))*(plcl(i)-p(i,icb(i)))
     &  /(p(i,icb(i))-p(i,icb(i)+1))
 610  continue

c-------------------------------------------------------------------
c --- Interpolate difference between lifted parcel and
c --- environmental temperatures to lifted condensation level
c-------------------------------------------------------------------
c
c dtpbl = average of tvp-tv in the PBL (k=nk to icb-1).
c
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
c
c-------------------------------------------------------------------
c --- Adjust cloud base mass flux
c-------------------------------------------------------------------
c
      do 650 i=1,ncum
       work(i)=cbmf(i)
       cbmf(i)=max(0.0,(1.0-damp)*cbmf(i)+0.1*alpha*dtmin(i))
       if((work(i).eq.0.0).and.(cbmf(i).eq.0.0))then
         iflag(i)=3
       endif
 650  continue

       return
       end

      SUBROUTINE cv_mixing(nloc,ncum,nd,icb,nk,inb,inb1
     :                    ,ph,t,q,qs,u,v,h,lv,qnk
     :                    ,hp,tv,tvp,ep,clw,cbmf
     :                    ,m,ment,qent,uent,vent,nent,sij,elij)
            use cvthermo
      implicit none

      include "cvparam.h"

c inputs:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc), inb1(nloc), nk(nloc)
      real cbmf(nloc), qnk(nloc)
      real ph(nloc,nd+1)
      real t(nloc,nd), q(nloc,nd), qs(nloc,nd), lv(nloc,nd)
      real u(nloc,nd), v(nloc,nd), h(nloc,nd), hp(nloc,nd)
      real tv(nloc,nd), tvp(nloc,nd), ep(nloc,nd), clw(nloc,nd)

c outputs:
      integer nent(nloc,nd)
      real m(nloc,nd), ment(nloc,nd,nd), qent(nloc,nd,nd)
      real uent(nloc,nd,nd), vent(nloc,nd,nd)
      real sij(nloc,nd,nd), elij(nloc,nd,nd)

c local variables:
      integer i, j, k, ij
      integer num1, num2
      real dbo, qti, bf2, anum, denom, dei, altem, cwat, stemp
      real alt, qp1, smid, sjmin, sjmax, delp, delm
      real work(nloc), asij(nloc), smin(nloc), scrit(nloc)
      real bsum(nloc,nd)
      logical lwork(nloc)

c=====================================================================
c --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
c=====================================================================
c
        do 360 i=1,ncum*nlp
          nent(i,1)=0
          m(i,1)=0.0
 360    continue
c
      do 400 k=1,nlp
       do 390 j=1,nlp
          do 385 i=1,ncum
            qent(i,k,j)=q(i,j)
            uent(i,k,j)=u(i,j)
            vent(i,k,j)=v(i,j)
            elij(i,k,j)=0.0
            ment(i,k,j)=0.0
            sij(i,k,j)=0.0
 385      continue
 390    continue
 400  continue
c
c-------------------------------------------------------------------
c --- Calculate rates of mixing,  m(i)
c-------------------------------------------------------------------
c
      call zilch(work,ncum)
c
      do 670 j=minorig+1,nl
        do 660 i=1,ncum
          if((j.ge.(icb(i)+1)).and.(j.le.inb(i)))then
             k=min(j,inb1(i))
             dbo=abs(tv(i,k+1)-tvp(i,k+1)-tv(i,k-1)+tvp(i,k-1))
     &       +entp*0.04*(ph(i,k)-ph(i,k+1))
             work(i)=work(i)+dbo
             m(i,j)=cbmf(i)*dbo
          endif
 660    continue
 670  continue
      do 690 k=minorig+1,nl
        do 680 i=1,ncum
          if((k.ge.(icb(i)+1)).and.(k.le.inb(i)))then
            m(i,k)=m(i,k)/work(i)
          endif
 680    continue
 690  continue
c
c
c=====================================================================
c --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
c --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
c --- FRACTION (sij)
c=====================================================================
c
c
       do 750 i=minorig+1,nl
         do 710 j=minorig+1,nl
           do 700 ij=1,ncum
             if((i.ge.(icb(ij)+1)).and.(j.ge.icb(ij))
     &         .and.(i.le.inb(ij)).and.(j.le.inb(ij)))then
               qti=qnk(ij)-ep(ij,i)*clw(ij,i)
               bf2=1.+lv(ij,j)*lv(ij,j)*qs(ij,j)
     &         /(rrv*t(ij,j)*t(ij,j)*cpd)
               anum=h(ij,j)-hp(ij,i)+(cpv-cpd)*t(ij,j)*(qti-q(ij,j))
               denom=h(ij,i)-hp(ij,i)+(cpd-cpv)*(q(ij,i)-qti)*t(ij,j)
               dei=denom
               if(abs(dei).lt.0.01)dei=0.01
               sij(ij,i,j)=anum/dei
               sij(ij,i,i)=1.0
               altem=sij(ij,i,j)*q(ij,i)+(1.-sij(ij,i,j))*qti-qs(ij,j)
               altem=altem/bf2
               cwat=clw(ij,j)*(1.-ep(ij,j))
               stemp=sij(ij,i,j)
               if((stemp.lt.0.0.or.stemp.gt.1.0.or.
     1           altem.gt.cwat).and.j.gt.i)then
                 anum=anum-lv(ij,j)*(qti-qs(ij,j)-cwat*bf2)
                 denom=denom+lv(ij,j)*(q(ij,i)-qti)
                 if(abs(denom).lt.0.01)denom=0.01
                 sij(ij,i,j)=anum/denom
                 altem=sij(ij,i,j)*q(ij,i)+(1.-sij(ij,i,j))*qti-qs(ij,j)
                 altem=altem-(bf2-1.)*cwat
               endif
               if(sij(ij,i,j).gt.0.0.and.sij(ij,i,j).lt.0.9)then
                 qent(ij,i,j)=sij(ij,i,j)*q(ij,i)
     &                        +(1.-sij(ij,i,j))*qti
                 uent(ij,i,j)=sij(ij,i,j)*u(ij,i)
     &                        +(1.-sij(ij,i,j))*u(ij,nk(ij))
                 vent(ij,i,j)=sij(ij,i,j)*v(ij,i)
     &                        +(1.-sij(ij,i,j))*v(ij,nk(ij))
                 elij(ij,i,j)=altem
                 elij(ij,i,j)=max(0.0,elij(ij,i,j))
                 ment(ij,i,j)=m(ij,i)/(1.-sij(ij,i,j))
                 nent(ij,i)=nent(ij,i)+1
               endif
             sij(ij,i,j)=max(0.0,sij(ij,i,j))
             sij(ij,i,j)=min(1.0,sij(ij,i,j))
             endif
  700      continue
  710    continue
c
c   ***   If no air can entrain at level i assume that updraft detrains  ***
c   ***   at that level and calculate detrained air flux and properties  ***
c
           do 740 ij=1,ncum
             if((i.ge.(icb(ij)+1)).and.(i.le.inb(ij))
     &       .and.(nent(ij,i).eq.0))then
               ment(ij,i,i)=m(ij,i)
               qent(ij,i,i)=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
               uent(ij,i,i)=u(ij,nk(ij))
               vent(ij,i,i)=v(ij,nk(ij))
               elij(ij,i,i)=clw(ij,i)
               sij(ij,i,i)=1.0
             endif
 740       continue
 750   continue
c
      do 770 i=1,ncum
        sij(i,inb(i),inb(i))=1.0
 770  continue
c
c=====================================================================
c   ---  NORMALIZE ENTRAINED AIR MASS FLUXES
c   ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
c=====================================================================
c
       call zilch(bsum,ncum*nlp)
       do 780 ij=1,ncum
         lwork(ij)=.false.
 780   continue
       do 789 i=minorig+1,nl
c
         num1=0
         do 779 ij=1,ncum
           if((i.ge.icb(ij)+1).and.(i.le.inb(ij)))num1=num1+1
 779     continue
         if(num1.le.0)go to 789
c
           do 781 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij)))then
                lwork(ij)=(nent(ij,i).ne.0)
                qp1=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
                anum=h(ij,i)-hp(ij,i)-lv(ij,i)*(qp1-qs(ij,i))
                denom=h(ij,i)-hp(ij,i)+lv(ij,i)*(q(ij,i)-qp1)
                if(abs(denom).lt.0.01)denom=0.01
                scrit(ij)=anum/denom
                alt=qp1-qs(ij,i)+scrit(ij)*(q(ij,i)-qp1)
                if(scrit(ij).lt.0.0.or.alt.lt.0.0)scrit(ij)=1.0
                asij(ij)=0.0
                smin(ij)=1.0
             endif
 781       continue
         do 783 j=minorig,nl
c
         num2=0
         do 778 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &       .and.(j.ge.icb(ij)).and.(j.le.inb(ij))
     &       .and.lwork(ij))num2=num2+1
 778     continue
         if(num2.le.0)go to 783
c
           do 782 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &       .and.(j.ge.icb(ij)).and.(j.le.inb(ij)).and.lwork(ij))then
                  if(sij(ij,i,j).gt.0.0.and.sij(ij,i,j).lt.0.9)then
                    if(j.gt.i)then
                      smid=min(sij(ij,i,j),scrit(ij))
                      sjmax=smid
                      sjmin=smid
                        if(smid.lt.smin(ij)
     &                  .and.sij(ij,i,j+1).lt.smid)then
                          smin(ij)=smid
                          sjmax=min(sij(ij,i,j+1),sij(ij,i,j),scrit(ij))
                          sjmin=max(sij(ij,i,j-1),sij(ij,i,j))
                          sjmin=min(sjmin,scrit(ij))
                        endif
                    else
                      sjmax=max(sij(ij,i,j+1),scrit(ij))
                      smid=max(sij(ij,i,j),scrit(ij))
                      sjmin=0.0
                      if(j.gt.1)sjmin=sij(ij,i,j-1)
                      sjmin=max(sjmin,scrit(ij))
                    endif
                    delp=abs(sjmax-smid)
                    delm=abs(sjmin-smid)
                    asij(ij)=asij(ij)+(delp+delm)
     &                           *(ph(ij,j)-ph(ij,j+1))
                    ment(ij,i,j)=ment(ij,i,j)*(delp+delm)
     &                           *(ph(ij,j)-ph(ij,j+1))
                  endif
              endif
  782    continue
  783    continue
            do 784 ij=1,ncum
            if((i.ge.icb(ij)+1).and.(i.le.inb(ij)).and.lwork(ij))then
               asij(ij)=max(1.0e-21,asij(ij))
               asij(ij)=1.0/asij(ij)
               bsum(ij,i)=0.0
            endif
 784        continue
            do 786 j=minorig,nl+1
              do 785 ij=1,ncum
                if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &          .and.(j.ge.icb(ij)).and.(j.le.inb(ij))
     &          .and.lwork(ij))then
                   ment(ij,i,j)=ment(ij,i,j)*asij(ij)
                   bsum(ij,i)=bsum(ij,i)+ment(ij,i,j)
                endif
 785     continue
 786     continue
             do 787 ij=1,ncum
               if((i.ge.icb(ij)+1).and.(i.le.inb(ij))
     &         .and.(bsum(ij,i).lt.1.0e-18).and.lwork(ij))then
                 nent(ij,i)=0
                 ment(ij,i,i)=m(ij,i)
                 qent(ij,i,i)=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
                 uent(ij,i,i)=u(ij,nk(ij))
                 vent(ij,i,i)=v(ij,nk(ij))
                 elij(ij,i,i)=clw(ij,i)
                 sij(ij,i,i)=1.0
               endif
  787        continue
  789  continue
c
       return
       end

      SUBROUTINE cv_unsat(nloc,ncum,nd,inb,t,q,qs,gz,u,v,p,ph
     :                  ,h,lv,ep,sigp,clw,m,ment,elij
     :                  ,iflag,mp,qp,up,vp,wt,water,evap)
            use cvthermo
      implicit none


      include "cvparam.h"

c inputs:
      integer ncum, nd, nloc
      integer inb(nloc)
      real t(nloc,nd), q(nloc,nd), qs(nloc,nd)
      real gz(nloc,nd), u(nloc,nd), v(nloc,nd)
      real p(nloc,nd), ph(nloc,nd+1), h(nloc,nd)
      real lv(nloc,nd), ep(nloc,nd), sigp(nloc,nd), clw(nloc,nd)
      real m(nloc,nd), ment(nloc,nd,nd), elij(nloc,nd,nd)

c outputs:
      integer iflag(nloc) ! also an input
      real mp(nloc,nd), qp(nloc,nd), up(nloc,nd), vp(nloc,nd)
      real water(nloc,nd), evap(nloc,nd), wt(nloc,nd)

c local variables:
      integer i,j,k,ij,num1
      integer jtt(nloc)
      real awat, coeff, qsm, afac, sigt, b6, c6, revap
      real dhdp, fac, qstm, rat
      real wdtrain(nloc)
      logical lwork(nloc)

c=====================================================================
c --- PRECIPITATING DOWNDRAFT CALCULATION
c=====================================================================
c
c Initializations:
c
         do i = 1, ncum
         do k = 1, nl+1
          wt(i,k) = omtsnow
          mp(i,k) = 0.0
          evap(i,k) = 0.0
          water(i,k) = 0.0
         enddo
         enddo

         do 420 i=1,ncum
          qp(i,1)=q(i,1)
          up(i,1)=u(i,1)
          vp(i,1)=v(i,1)
 420     continue

         do 440 k=2,nl+1
         do 430 i=1,ncum
          qp(i,k)=q(i,k-1)
          up(i,k)=u(i,k-1)
          vp(i,k)=v(i,k-1)
 430     continue
 440     continue


c   ***  Check whether ep(inb)=0, if so, skip precipitating    ***
c   ***             downdraft calculation                      ***
c
c
c   ***  Integrate liquid water equation to find condensed water   ***
c   ***                and condensed water flux                    ***
c
c
      do 890 i=1,ncum
        jtt(i)=2
        if(ep(i,inb(i)).le.0.0001)iflag(i)=2
        if(iflag(i).eq.0)then
          lwork(i)=.true.
        else
          lwork(i)=.false.
        endif
 890  continue
c
c    ***                    Begin downdraft loop                    ***
c
c
        call zilch(wdtrain,ncum)
        do 899 i=nl+1,1,-1
c
          num1=0
          do 879 ij=1,ncum
            if((i.le.inb(ij)).and.lwork(ij))num1=num1+1
 879      continue
          if(num1.le.0)go to 899
c
c
c    ***        Calculate detrained precipitation             ***
c
          do 891 ij=1,ncum
            if((i.le.inb(ij)).and.(lwork(ij)))then
            wdtrain(ij)=g*ep(ij,i)*m(ij,i)*clw(ij,i)
            endif
 891      continue
c
          if(i.gt.1)then
            do 893 j=1,i-1
              do 892 ij=1,ncum
                if((i.le.inb(ij)).and.(lwork(ij)))then
                  awat=elij(ij,j,i)-(1.-ep(ij,i))*clw(ij,i)
                  awat=max(0.0,awat)
                  wdtrain(ij)=wdtrain(ij)+g*awat*ment(ij,j,i)
                endif
 892          continue
 893      continue
          endif
c
c    ***    Find rain water and evaporation using provisional   ***
c    ***              estimates of qp(i)and qp(i-1)             ***
c
c
c  ***  Value of terminal velocity and coeffecient of evaporation for snow   ***
c
          do 894 ij=1,ncum
            if((i.le.inb(ij)).and.(lwork(ij)))then
            coeff=coeffs
            wt(ij,i)=omtsnow
c
c  ***  Value of terminal velocity and coeffecient of evaporation for rain   ***
c
            if(t(ij,i).gt.273.0)then
              coeff=coeffr
              wt(ij,i)=omtrain
            endif
            qsm=0.5*(q(ij,i)+qp(ij,i+1))
            afac=coeff*ph(ij,i)*(qs(ij,i)-qsm)
     &       /(1.0e4+2.0e3*ph(ij,i)*qs(ij,i))
            afac=max(afac,0.0)
            sigt=sigp(ij,i)
            sigt=max(0.0,sigt)
            sigt=min(1.0,sigt)
            b6=100.*(ph(ij,i)-ph(ij,i+1))*sigt*afac/wt(ij,i)
            c6=(water(ij,i+1)*wt(ij,i+1)+wdtrain(ij)/sigd)/wt(ij,i)
            revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
            evap(ij,i)=sigt*afac*revap
            water(ij,i)=revap*revap
c
c    ***  Calculate precipitating downdraft mass flux under     ***
c    ***              hydrostatic approximation                 ***
c
            if(i.gt.1)then
              dhdp=(h(ij,i)-h(ij,i-1))/(p(ij,i-1)-p(ij,i))
              dhdp=max(dhdp,10.0)
              mp(ij,i)=100.*ginv*lv(ij,i)*sigd*evap(ij,i)/dhdp
              mp(ij,i)=max(mp(ij,i),0.0)
c
c   ***   Add small amount of inertia to downdraft              ***
c
              fac=20.0/(ph(ij,i-1)-ph(ij,i))
              mp(ij,i)=(fac*mp(ij,i+1)+mp(ij,i))/(1.+fac)
c
c    ***      Force mp to decrease linearly to zero                 ***
c    ***      between about 950 mb and the surface                  ***
c
              if(p(ij,i).gt.(0.949*p(ij,1)))then
                 jtt(ij)=max(jtt(ij),i)
                 mp(ij,i)=mp(ij,jtt(ij))*(p(ij,1)-p(ij,i))
     &           /(p(ij,1)-p(ij,jtt(ij)))
              endif
            endif
c
c    ***       Find mixing ratio of precipitating downdraft     ***
c
            if(i.ne.inb(ij))then
              if(i.eq.1)then
                qstm=qs(ij,1)
              else
                qstm=qs(ij,i-1)
              endif
              if(mp(ij,i).gt.mp(ij,i+1))then
                 rat=mp(ij,i+1)/mp(ij,i)
                 qp(ij,i)=qp(ij,i+1)*rat+q(ij,i)*(1.0-rat)+100.*ginv*
     &             sigd*(ph(ij,i)-ph(ij,i+1))*(evap(ij,i)/mp(ij,i))
                 up(ij,i)=up(ij,i+1)*rat+u(ij,i)*(1.-rat)
                 vp(ij,i)=vp(ij,i+1)*rat+v(ij,i)*(1.-rat)
               else
                 if(mp(ij,i+1).gt.0.0)then
                   qp(ij,i)=(gz(ij,i+1)-gz(ij,i)
     &               +qp(ij,i+1)*(lv(ij,i+1)+t(ij,i+1)
     &               *(cl-cpd))+cpd*(t(ij,i+1)-t(ij,i)))
     &               /(lv(ij,i)+t(ij,i)*(cl-cpd))
                   up(ij,i)=up(ij,i+1)
                   vp(ij,i)=vp(ij,i+1)
                 endif
              endif
              qp(ij,i)=min(qp(ij,i),qstm)
              qp(ij,i)=max(qp(ij,i),0.0)
            endif
            endif
 894      continue
 899    continue
c
        return
        end

      SUBROUTINE cv_yield(nloc,ncum,nd,nk,icb,inb,delt
     :             ,t,q,u,v,gz,p,ph,h,hp,lv,cpn
     :             ,ep,clw,frac,m,mp,qp,up,vp
     :             ,wt,water,evap
     :             ,ment,qent,uent,vent,nent,elij
     :             ,tv,tvp
     o             ,iflag,wd,qprime,tprime
     o             ,precip,cbmf,ft,fq,fu,fv,Ma,qcondc)
            use cvthermo
      implicit none

      include "cvparam.h"

c inputs
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

c outputs
      integer iflag(nloc)  ! also an input
      real cbmf(nloc)      ! also an input
      real wd(nloc), tprime(nloc), qprime(nloc)
      real precip(nloc)
      real ft(nloc,nd), fq(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real Ma(nloc,nd)
      real qcondc(nloc,nd)

c local variables
      integer i,j,ij,k,num1
      real dpinv,cpinv,awat,fqold,ftold,fuold,fvold,delti
      real work(nloc), am(nloc),amp1(nloc),ad(nloc)
      real ents(nloc), uav(nloc),vav(nloc),lvcp(nloc,nd)
      real qcond(nloc,nd), nqcond(nloc,nd), wa(nloc,nd) ! cld
      real siga(nloc,nd), ax(nloc,nd), mac(nloc,nd)     ! cld

 
c -- initializations:

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

c
c   ***  Calculate surface precipitation in mm/day     ***
c
        do 1190 i=1,ncum
          if(iflag(i).le.1)then
            precip(i) = wt(i,1)*sigd*water(i,1)*86400/g
          endif
 1190   continue
c
c
c   ***  Calculate downdraft velocity scale and surface temperature and  ***
c   ***                    water vapor fluctuations                      ***
c
      do i=1,ncum
       wd(i)=betad*abs(mp(i,icb(i)))*0.01*rrd*t(i,icb(i))
     :           /(sigd*p(i,icb(i)))
       qprime(i)=0.5*(qp(i,1)-q(i,1))
       tprime(i)=lv0*qprime(i)/cpd
      enddo
c
c   ***  Calculate tendencies of lowest level potential temperature  ***
c   ***                      and mixing ratio                        ***
c
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
          ft(i,1)=ft(i,1)+g*work(i)*am(i)*(t(i,2)-t(i,1)
     &    +(gz(i,2)-gz(i,1))/cpn(i,1))
          ft(i,1)=ft(i,1)-lvcp(i,1)*sigd*evap(i,1)
          ft(i,1)=ft(i,1)+sigd*wt(i,2)*(cl-cpd)*water(i,2)*(t(i,2)
     &     -t(i,1))*work(i)/cpn(i,1)
          fq(i,1)=fq(i,1)+g*mp(i,2)*(qp(i,2)-q(i,1))*
     &    work(i)+sigd*evap(i,1)
          fq(i,1)=fq(i,1)+g*am(i)*(q(i,2)-q(i,1))*work(i)
          fu(i,1)=fu(i,1)+g*work(i)*(mp(i,2)*(up(i,2)-u(i,1))
     &    +am(i)*(u(i,2)-u(i,1)))
          fv(i,1)=fv(i,1)+g*work(i)*(mp(i,2)*(vp(i,2)-v(i,1))
     &    +am(i)*(v(i,2)-v(i,1)))
 1240   continue
        do 1260 j=2,nl
           do 1250 i=1,ncum
             if(j.le.inb(i))then
               fq(i,1)=fq(i,1)
     &                 +g*work(i)*ment(i,j,1)*(qent(i,j,1)-q(i,1))
               fu(i,1)=fu(i,1)
     &                 +g*work(i)*ment(i,j,1)*(uent(i,j,1)-u(i,1))
               fv(i,1)=fv(i,1)
     &                 +g*work(i)*ment(i,j,1)*(vent(i,j,1)-v(i,1))
             endif
 1250      continue
 1260   continue
c
c   ***  Calculate tendencies of potential temperature and mixing ratio  ***
c   ***               at levels above the lowest level                   ***
c
c   ***  First find the net saturated updraft and downdraft mass fluxes  ***
c   ***                      through each level                          ***
c
        do 1500 i=2,nl+1
c
          num1=0
          do 1265 ij=1,ncum
            if(i.le.inb(ij))num1=num1+1
 1265     continue
          if(num1.le.0)go to 1500
c
          call zilch(amp1,ncum)
          call zilch(ad,ncum)
c
          do 1280 k=i+1,nl+1
            do 1270 ij=1,ncum
              if((i.ge.nk(ij)).and.(i.le.inb(ij))
     &            .and.(k.le.(inb(ij)+1)))then
                amp1(ij)=amp1(ij)+m(ij,k)
              endif
 1270         continue
 1280     continue
c
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
c
          do 1350 ij=1,ncum
          if(i.le.inb(ij))then
            dpinv=0.01/(ph(ij,i)-ph(ij,i+1))
            cpinv=1.0/cpn(ij,i)
c
            ft(ij,i)=ft(ij,i)
     &       +g*dpinv*(amp1(ij)*(t(ij,i+1)-t(ij,i)
     &       +(gz(ij,i+1)-gz(ij,i))*cpinv)
     &       -ad(ij)*(t(ij,i)-t(ij,i-1)+(gz(ij,i)-gz(ij,i-1))*cpinv))
     &       -sigd*lvcp(ij,i)*evap(ij,i)
            ft(ij,i)=ft(ij,i)+g*dpinv*ment(ij,i,i)*(hp(ij,i)-h(ij,i)+
     &        t(ij,i)*(cpv-cpd)*(q(ij,i)-qent(ij,i,i)))*cpinv
            ft(ij,i)=ft(ij,i)+sigd*wt(ij,i+1)*(cl-cpd)*water(ij,i+1)*
     &        (t(ij,i+1)-t(ij,i))*dpinv*cpinv
            fq(ij,i)=fq(ij,i)+g*dpinv*(amp1(ij)*(q(ij,i+1)-q(ij,i))-
     &        ad(ij)*(q(ij,i)-q(ij,i-1)))
            fu(ij,i)=fu(ij,i)+g*dpinv*(amp1(ij)*(u(ij,i+1)-u(ij,i))-
     &        ad(ij)*(u(ij,i)-u(ij,i-1)))
            fv(ij,i)=fv(ij,i)+g*dpinv*(amp1(ij)*(v(ij,i+1)-v(ij,i))-
     &        ad(ij)*(v(ij,i)-v(ij,i-1)))
         endif
 1350    continue
         do 1370 k=1,i-1
           do 1360 ij=1,ncum
             if(i.le.inb(ij))then
               awat=elij(ij,k,i)-(1.-ep(ij,i))*clw(ij,i)
               awat=max(awat,0.0)
               fq(ij,i)=fq(ij,i)
     &         +g*dpinv*ment(ij,k,i)*(qent(ij,k,i)-awat-q(ij,i))
               fu(ij,i)=fu(ij,i)
     &         +g*dpinv*ment(ij,k,i)*(uent(ij,k,i)-u(ij,i))
               fv(ij,i)=fv(ij,i)
     &         +g*dpinv*ment(ij,k,i)*(vent(ij,k,i)-v(ij,i))
c (saturated updrafts resulting from mixing)               ! cld
               qcond(ij,i)=qcond(ij,i)+(elij(ij,k,i)-awat) ! cld
               nqcond(ij,i)=nqcond(ij,i)+1.                ! cld
             endif
 1360      continue
 1370    continue
         do 1390 k=i,nl+1
           do 1380 ij=1,ncum
             if((i.le.inb(ij)).and.(k.le.inb(ij)))then
               fq(ij,i)=fq(ij,i)
     &                  +g*dpinv*ment(ij,k,i)*(qent(ij,k,i)-q(ij,i))
               fu(ij,i)=fu(ij,i)
     &                  +g*dpinv*ment(ij,k,i)*(uent(ij,k,i)-u(ij,i))
               fv(ij,i)=fv(ij,i)
     &                  +g*dpinv*ment(ij,k,i)*(vent(ij,k,i)-v(ij,i))
             endif
 1380      continue
 1390    continue
          do 1400 ij=1,ncum
           if(i.le.inb(ij))then
             fq(ij,i)=fq(ij,i)
     &                +sigd*evap(ij,i)+g*(mp(ij,i+1)*
     &                (qp(ij,i+1)-q(ij,i))
     &                -mp(ij,i)*(qp(ij,i)-q(ij,i-1)))*dpinv
             fu(ij,i)=fu(ij,i)
     &                +g*(mp(ij,i+1)*(up(ij,i+1)-u(ij,i))-mp(ij,i)*
     &                (up(ij,i)-u(ij,i-1)))*dpinv
             fv(ij,i)=fv(ij,i)
     &               +g*(mp(ij,i+1)*(vp(ij,i+1)-v(ij,i))-mp(ij,i)*
     &               (vp(ij,i)-v(ij,i-1)))*dpinv
C (saturated downdrafts resulting from mixing)               ! cld
            do k=i+1,inb(ij)                                 ! cld
             qcond(ij,i)=qcond(ij,i)+elij(ij,k,i)            ! cld
             nqcond(ij,i)=nqcond(ij,i)+1.                    ! cld
            enddo                                            ! cld
C (particular case: no detraining level is found)            ! cld
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
c
c   *** Adjust tendencies at top of convection layer to reflect  ***
c   ***       actual position of the level zero cape             ***
c
        do 503 ij=1,ncum
        fqold=fq(ij,inb(ij))
        fq(ij,inb(ij))=fq(ij,inb(ij))*(1.-frac(ij))
        fq(ij,inb(ij)-1)=fq(ij,inb(ij)-1)
     &   +frac(ij)*fqold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))*lv(ij,inb(ij))
     &   /lv(ij,inb(ij)-1)
        ftold=ft(ij,inb(ij))
        ft(ij,inb(ij))=ft(ij,inb(ij))*(1.-frac(ij))
        ft(ij,inb(ij)-1)=ft(ij,inb(ij)-1)
     &   +frac(ij)*ftold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))*cpn(ij,inb(ij))
     &   /cpn(ij,inb(ij)-1)
        fuold=fu(ij,inb(ij))
        fu(ij,inb(ij))=fu(ij,inb(ij))*(1.-frac(ij))
        fu(ij,inb(ij)-1)=fu(ij,inb(ij)-1)
     &   +frac(ij)*fuold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
        fvold=fv(ij,inb(ij))
        fv(ij,inb(ij))=fv(ij,inb(ij))*(1.-frac(ij))
        fv(ij,inb(ij)-1)=fv(ij,inb(ij)-1)
     &  +frac(ij)*fvold*((ph(ij,inb(ij))-ph(ij,inb(ij)+1))/
     1   (ph(ij,inb(ij)-1)-ph(ij,inb(ij))))
 503    continue
c
c   ***   Very slightly adjust tendencies to force exact   ***
c   ***     enthalpy, momentum and tracer conservation     ***
c
        do 682 ij=1,ncum
        ents(ij)=0.0
        uav(ij)=0.0
        vav(ij)=0.0
        do 681 i=1,inb(ij)
         ents(ij)=ents(ij)
     &  +(cpn(ij,i)*ft(ij,i)+lv(ij,i)*fq(ij,i))*(ph(ij,i)-ph(ij,i+1))	
         uav(ij)=uav(ij)+fu(ij,i)*(ph(ij,i)-ph(ij,i+1))
         vav(ij)=vav(ij)+fv(ij,i)*(ph(ij,i)-ph(ij,i+1))
  681	continue
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
  641	continue
 642    continue
c
        do 1810 k=1,nl+1
          do 1800 i=1,ncum
            if((q(i,k)+delt*fq(i,k)).lt.0.0)iflag(i)=10
 1800     continue
 1810   continue
c
c
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

c
c   *** diagnose the in-cloud mixing ratio   ***            ! cld
c   ***           of condensed water         ***            ! cld
c                                                           ! cld
      DO ij=1,ncum                                          ! cld   
       do i=1,nd                                            ! cld 
        mac(ij,i)=0.0                                       ! cld   
        wa(ij,i)=0.0                                        ! cld
        siga(ij,i)=0.0                                      ! cld
       enddo                                                ! cld
       do i=nk(ij),inb(ij)                                  ! cld
       do k=i+1,inb(ij)+1                                   ! cld
        mac(ij,i)=mac(ij,i)+m(ij,k)                         ! cld
       enddo                                                ! cld
       enddo                                                ! cld
       do i=icb(ij),inb(ij)-1                               ! cld
        ax(ij,i)=0.                                         ! cld
        do j=icb(ij),i                                      ! cld
         ax(ij,i)=ax(ij,i)+rrd*(tvp(ij,j)-tv(ij,j))         ! cld   
     :       *(ph(ij,j)-ph(ij,j+1))/p(ij,j)                 ! cld   
        enddo                                               ! cld
        if (ax(ij,i).gt.0.0) then                           ! cld   
         wa(ij,i)=sqrt(2.*ax(ij,i))                         ! cld
        endif                                               ! cld
       enddo                                                ! cld
       do i=1,nl                                            ! cld
        if (wa(ij,i).gt.0.0)                                ! cld
     :    siga(ij,i)=mac(ij,i)/wa(ij,i)                     ! cld   
     :        *rrd*tvp(ij,i)/p(ij,i)/100./delta             ! cld   
        siga(ij,i) = min(siga(ij,i),1.0)                    ! cld
        qcondc(ij,i)=siga(ij,i)*clw(ij,i)*(1.-ep(ij,i))     ! cld   
     :          + (1.-siga(ij,i))*qcond(ij,i)               ! cld   
       enddo                                                ! cld
      ENDDO                                                 ! cld   

        return
        end

      SUBROUTINE cv_uncompress(nloc,len,ncum,nd,idcum
     :         ,iflag
     :         ,precip,cbmf
     :         ,ft,fq,fu,fv
     :         ,Ma,qcondc            
     :         ,iflag1
     :         ,precip1,cbmf1
     :         ,ft1,fq1,fu1,fv1
     :         ,Ma1,qcondc1            
     :                               )
      implicit none

      include "cvparam.h"

c inputs:
      integer len, ncum, nd, nloc
      integer idcum(nloc)
      integer iflag(nloc)
      real precip(nloc), cbmf(nloc)
      real ft(nloc,nd), fq(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real Ma(nloc,nd)
      real qcondc(nloc,nd) !cld

c outputs:
      integer iflag1(len)
      real precip1(len), cbmf1(len)
      real ft1(len,nd), fq1(len,nd), fu1(len,nd), fv1(len,nd)
      real Ma1(len,nd)
      real qcondc1(len,nd) !cld

c local variables:
      integer i,k

        do 2000 i=1,ncum
         precip1(idcum(i))=precip(i)
         cbmf1(idcum(i))=cbmf(i)
         iflag1(idcum(i))=iflag(i)
 2000   continue

        do 2020 k=1,nl
          do 2010 i=1,ncum
            ft1(idcum(i),k)=ft(i,k)
            fq1(idcum(i),k)=fq(i,k)
            fu1(idcum(i),k)=fu(i,k)
            fv1(idcum(i),k)=fv(i,k)
            Ma1(idcum(i),k)=Ma(i,k)
            qcondc1(idcum(i),k)=qcondc(i,k)
 2010     continue
 2020   continue

        return
        end

