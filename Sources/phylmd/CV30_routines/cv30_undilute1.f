
      SUBROUTINE cv30_undilute1(len,nd,t,q,qs,gz,plcl,p,nk,icb &
                             ,tp,tvp,clw,icbs)
            use cv30_param_m
            use cvthermo
      implicit none

!----------------------------------------------------------------
! Equivalent de TLIFT entre NK et ICB+1 inclus
!
! Differences with convect4:
!        - specify plcl in input
!       - icbs is the first level above LCL (may differ from icb)
!       - in the iterations, used x(icbs) instead x(icb)
!       - many minor differences in the iterations
!        - tvp is computed in only one time
!        - icbs: first level above Plcl (IMIN de TLIFT) in output
!       - if icbs=icb, compute also tp(icb+1),tvp(icb+1) & clw(icb+1)
!----------------------------------------------------------------


! inputs:
      integer, intent(in):: len, nd
      integer nk(len), icb(len)
      real, intent(in):: t(len,nd)
      real, intent(in):: q(len,nd), qs(len,nd), gz(len,nd)
      real, intent(in):: p(len,nd)
      real plcl(len) ! convect3

! outputs:
      real tp(len,nd), tvp(len,nd), clw(len,nd)

! local variables:
      integer i, k
      integer icb1(len), icbs(len), icbsmax2 ! convect3
      real tg, qg, alv, s, ahg, tc, denom, es
      real ah0(len), cpp(len)
      real tnk(len), qnk(len), gznk(len), ticb(len), gzicb(len)
      real qsicb(len) ! convect3
      real cpinv(len) ! convect3

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
! ori        ticb(i)=t(i,icb(i))
! ori        gzicb(i)=gz(i,icb(i))
 320  continue
!
!   ***  Calculate certain parcel quantities, including static energy   ***
!
      do 330 i=1,len
        ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) &
               +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
        cpp(i)=cpd*(1.-qnk(i))+qnk(i)*cpv
        cpinv(i)=1./cpp(i)
 330  continue
!
!   ***   Calculate lifted parcel quantities below cloud base   ***
!
        do i=1,len                      !convect3
         icb1(i)=MAX(icb(i),2)          !convect3
         icb1(i)=MIN(icb(i),nl)         !convect3
! if icb is below LCL, start loop at ICB+1:
! (icbs est le premier niveau au-dessus du LCL)
         icbs(i)=icb1(i)                !convect3
         if (plcl(i).lt.p(i,icb1(i))) then
             icbs(i)=MIN(icbs(i)+1,nl)  !convect3
         endif
        enddo                           !convect3

        do i=1,len                      !convect3
         ticb(i)=t(i,icbs(i))           !convect3
         gzicb(i)=gz(i,icbs(i))         !convect3
         qsicb(i)=qs(i,icbs(i))         !convect3
        enddo                           !convect3

!
! Re-compute icbsmax (icbsmax2):        !convect3
!                                       !convect3
      icbsmax2=2                        !convect3
      do 310 i=1,len                    !convect3
        icbsmax2=max(icbsmax2,icbs(i))  !convect3
 310  continue                          !convect3

! initialization outputs:

      do k=1,icbsmax2     ! convect3
       do i=1,len         ! convect3
        tp(i,k)  = 0.0    ! convect3
        tvp(i,k) = 0.0    ! convect3
        clw(i,k) = 0.0    ! convect3
       enddo              ! convect3
      enddo               ! convect3

! tp and tvp below cloud base:

        do 350 k=minorig,icbsmax2-1
          do 340 i=1,len
           tp(i,k)=tnk(i)-(gz(i,k)-gznk(i))*cpinv(i)
           tvp(i,k)=tp(i,k)*(1.+qnk(i)/eps-qnk(i)) !whole thing (convect3)
  340     continue
  350   continue
!
!    ***  Find lifted parcel quantities above cloud base    ***
!
        do 360 i=1,len
         tg=ticb(i)
! ori         qg=qs(i,icb(i))
         qg=qsicb(i) ! convect3
!debug         alv=lv0-clmcpv*(ticb(i)-t0)
         alv=lv0-clmcpv*(ticb(i)-273.15)
!
! First iteration.
!
! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
          s=cpd*(1.-qnk(i))+cl*qnk(i) &
            +alv*alv*qg/(rrv*ticb(i)*ticb(i)) ! convect3
          s=1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
!debug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
! ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          endif
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icbs(i))-es*(1.-eps))
!
! Second iteration.
!

! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
! ori          s=1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
!debug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
! ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          end if
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icbs(i))-es*(1.-eps))

         alv=lv0-clmcpv*(ticb(i)-273.15)

! ori c approximation here:
! ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
! ori     &   -gz(i,icb(i))-alv*qg)/cpd

! convect3: no approximation:
         tp(i,icbs(i))=(ah0(i)-gz(i,icbs(i))-alv*qg) &
                      /(cpd+(cl-cpd)*qnk(i))

! ori         clw(i,icb(i))=qnk(i)-qg
! ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         clw(i,icbs(i))=qnk(i)-qg
         clw(i,icbs(i))=max(0.0,clw(i,icbs(i)))

! convect3: (qg utilise au lieu du vrai mixing ratio rg)
         tvp(i,icbs(i))=tp(i,icbs(i))*(1.+qg/eps-qnk(i)) !whole thing

  360   continue
!
! ori      do 380 k=minorig,icbsmax2
! ori       do 370 i=1,len
! ori         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
! ori 370   continue
! ori 380  continue
!

! -- The following is only for convect3:
!
! * icbs is the first level above the LCL:
!    if plcl<p(icb), then icbs=icb+1
!    if plcl>p(icb), then icbs=icb
!
! * the routine above computes tvp from minorig to icbs (included).
!
! * to compute buoybase (in cv30_trigger.F), both tvp(icb) and tvp(icb+1)
!    must be known. This is the case if icbs=icb+1, but not if icbs=icb.
!
! * therefore, in the case icbs=icb, we compute tvp at level icb+1
!   (tvp at other levels will be computed in cv30_undilute2.F)
!

        do i=1,len
         ticb(i)=t(i,icb(i)+1)
         gzicb(i)=gz(i,icb(i)+1)
         qsicb(i)=qs(i,icb(i)+1)
        enddo

        do 460 i=1,len
         tg=ticb(i)
         qg=qsicb(i) ! convect3
!debug         alv=lv0-clmcpv*(ticb(i)-t0)
         alv=lv0-clmcpv*(ticb(i)-273.15)
!
! First iteration.
!
! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
          s=cpd*(1.-qnk(i))+cl*qnk(i) &
            +alv*alv*qg/(rrv*ticb(i)*ticb(i)) ! convect3
          s=1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
!debug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
! ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          endif
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icb(i)+1)-es*(1.-eps))
!
! Second iteration.
!

! ori          s=cpd+alv*alv*qg/(rrv*ticb(i)*ticb(i))
! ori          s=1./s
! ori          ahg=cpd*tg+(cl-cpd)*qnk(i)*ticb(i)+alv*qg+gzicb(i)
          ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i) ! convect3
          tg=tg+s*(ah0(i)-ahg)
! ori          tg=max(tg,35.0)
!debug          tc=tg-t0
          tc=tg-273.15
          denom=243.5+tc
          denom=MAX(denom,1.0) ! convect3
! ori          if(tc.ge.0.0)then
           es=6.112*exp(17.67*tc/denom)
! ori          else
! ori           es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
! ori          end if
! ori          qg=eps*es/(p(i,icb(i))-es*(1.-eps))
          qg=eps*es/(p(i,icb(i)+1)-es*(1.-eps))

         alv=lv0-clmcpv*(ticb(i)-273.15)

! ori c approximation here:
! ori         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i)
! ori     &   -gz(i,icb(i))-alv*qg)/cpd

! convect3: no approximation:
         tp(i,icb(i)+1)=(ah0(i)-gz(i,icb(i)+1)-alv*qg) &
                      /(cpd+(cl-cpd)*qnk(i))

! ori         clw(i,icb(i))=qnk(i)-qg
! ori         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         clw(i,icb(i)+1)=qnk(i)-qg
         clw(i,icb(i)+1)=max(0.0,clw(i,icb(i)+1))

! convect3: (qg utilise au lieu du vrai mixing ratio rg)
         tvp(i,icb(i)+1)=tp(i,icb(i)+1)*(1.+qg/eps-qnk(i)) !whole thing

  460   continue

      return
      end
