
      SUBROUTINE cv_prelim(len,nd,ndp1,t,q,p,ph &
                          ,lv,cpn,tv,gz,h,hm)
            use cvthermo
            use cvparam
      implicit none

!=====================================================================
! --- CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY & STATIC ENERGY
!=====================================================================

! inputs:
      integer, intent(in):: len, nd, ndp1
      real, intent(in):: t(len,nd)
      real q(len,nd), p(len,nd), ph(len,ndp1)

! outputs:
      real lv(len,nd), cpn(len,nd), tv(len,nd)
      real gz(len,nd), h(len,nd), hm(len,nd)

! local variables:
      integer k, i
      real cpx(len,nd)



      do 110 k=1,nlp
        do 100 i=1,len
          lv(i,k)= lv0-clmcpv*(t(i,k)-t0)
          cpn(i,k)=cpd*(1.0-q(i,k))+cpv*q(i,k)
          cpx(i,k)=cpd*(1.0-q(i,k))+cl*q(i,k)
          tv(i,k)=t(i,k)*(1.0+q(i,k)*epsim1)
 100    continue
 110  continue
!
! gz = phi at the full levels (same as p).
!
      do 120 i=1,len
        gz(i,1)=0.0
 120  continue
      do 140 k=2,nlp
        do 130 i=1,len
          gz(i,k)=gz(i,k-1)+hrd*(tv(i,k-1)+tv(i,k)) &
               *(p(i,k-1)-p(i,k))/ph(i,k)
 130    continue
 140  continue
!
! h  = phi + cpT (dry static energy).
! hm = phi + cp(T-Tbase)+Lq
!
      do 170 k=1,nlp
        do 160 i=1,len
          h(i,k)=gz(i,k)+cpn(i,k)*t(i,k)
          hm(i,k)=gz(i,k)+cpx(i,k)*(t(i,k)-t(i,1))+lv(i,k)*q(i,k)
 160    continue
 170  continue

      return
      end
