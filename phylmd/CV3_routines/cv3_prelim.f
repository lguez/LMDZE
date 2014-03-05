
      SUBROUTINE cv3_prelim(len,nd,ndp1,t,q,p,ph &
                          ,lv,cpn,tv,gz,h,hm,th)
            use cv3_param_m
            use cvthermo
      implicit none

!=====================================================================
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
! "ori": from convect4.3 (vectorized)
! "convect3": to be exactly consistent with convect3
!=====================================================================

! inputs:
      integer len, nd, ndp1
      real, intent(in):: t(len,nd)
      real, intent(in):: q(len,nd)
      real p(len,nd), ph(len,ndp1)

! outputs:
      real lv(len,nd), cpn(len,nd), tv(len,nd)
      real gz(len,nd), h(len,nd), hm(len,nd)
      real th(len,nd)

! local variables:
      integer k, i
      real rdcp
      real tvx,tvy ! convect3
      real cpx(len,nd)



! ori      do 110 k=1,nlp
      do 110 k=1,nl ! convect3
        do 100 i=1,len
!debug          lv(i,k)= lv0-clmcpv*(t(i,k)-t0)
          lv(i,k)= lv0-clmcpv*(t(i,k)-273.15)
          cpn(i,k)=cpd*(1.0-q(i,k))+cpv*q(i,k)
          cpx(i,k)=cpd*(1.0-q(i,k))+cl*q(i,k)
! ori          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
          tv(i,k)=t(i,k)*(1.0+q(i,k)/eps-q(i,k))
          rdcp=(rrd*(1.-q(i,k))+q(i,k)*rrv)/cpn(i,k)
          th(i,k)=t(i,k)*(1000.0/p(i,k))**rdcp
 100    continue
 110  continue
!
! gz = phi at the full levels (same as p).
!
      do 120 i=1,len
        gz(i,1)=0.0
 120  continue
! ori      do 140 k=2,nlp
      do 140 k=2,nl ! convect3
        do 130 i=1,len
        tvx=t(i,k)*(1.+q(i,k)/eps-q(i,k))       !convect3
        tvy=t(i,k-1)*(1.+q(i,k-1)/eps-q(i,k-1)) !convect3
        gz(i,k)=gz(i,k-1)+0.5*rrd*(tvx+tvy)     &
                *(p(i,k-1)-p(i,k))/ph(i,k)      !convect3

! ori         gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k))
! ori    &         *(p(i,k-1)-p(i,k))/ph(i,k)
 130    continue
 140  continue
!
! h  = phi + cpT (dry static energy).
! hm = phi + cp(T-Tbase)+Lq
!
! ori      do 170 k=1,nlp
      do 170 k=1,nl ! convect3
        do 160 i=1,len
          h(i,k)=gz(i,k)+cpn(i,k)*t(i,k)
          hm(i,k)=gz(i,k)+cpx(i,k)*(t(i,k)-t(i,1))+lv(i,k)*q(i,k)
 160    continue
 170  continue

      return
      end
