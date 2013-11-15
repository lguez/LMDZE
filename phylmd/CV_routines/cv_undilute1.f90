
      SUBROUTINE cv_undilute1(len,nd,t,q,qs,gz,p,nk,icb,icbmax &
                             ,tp,tvp,clw)
            use cvthermo
            use cvparam
      implicit none


! inputs:
      integer len, nd
      integer nk(len), icb(len), icbmax
      real, intent(in):: t(len,nd)
      real q(len,nd), qs(len,nd), gz(len,nd)
      real p(len,nd)

! outputs:
      real tp(len,nd), tvp(len,nd), clw(len,nd)

! local variables:
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
!
!   ***  Calculate certain parcel quantities, including static energy   ***
!
      do 330 i=1,len
        ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) &
               +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
        cpp(i)=cpd*(1.-qnk(i))+qnk(i)*cpv
 330  continue
!
!   ***   Calculate lifted parcel quantities below cloud base   ***
!
        do 350 k=minorig,icbmax-1
          do 340 i=1,len
           tp(i,k)=tnk(i)-(gz(i,k)-gznk(i))/cpp(i)
           tvp(i,k)=tp(i,k)*(1.+qnk(i)*epsi)
  340     continue
  350   continue
!
!    ***  Find lifted parcel quantities above cloud base    ***
!
        do 360 i=1,len
         tg=ticb(i)
         qg=qs(i,icb(i))
         alv=lv0-clmcpv*(ticb(i)-t0)
!
! First iteration.
!
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
!
! Second iteration.
!
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
!
         alv=lv0-clmcpv*(ticb(i)-273.15)
         tp(i,icb(i))=(ah0(i)-(cl-cpd)*qnk(i)*ticb(i) &
         -gz(i,icb(i))-alv*qg)/cpd
         clw(i,icb(i))=qnk(i)-qg
         clw(i,icb(i))=max(0.0,clw(i,icb(i)))
         rg=qg/(1.-qnk(i))
         tvp(i,icb(i))=tp(i,icb(i))*(1.+rg*epsi)
  360   continue
!
      do 380 k=minorig,icbmax
       do 370 i=1,len
         tvp(i,k)=tvp(i,k)-tp(i,k)*qnk(i)
 370   continue
 380  continue
!
      return
      end
