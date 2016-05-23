module cv30_undilute1_m

  implicit none

contains

  SUBROUTINE cv30_undilute1(t1, q1, qs1, gz1, plcl1, p1, nk1, icb1, tp1, &
       tvp1, clw1, icbs1)

    ! UNDILUTE (ADIABATIC) UPDRAFT / 1st part
    ! (up through ICB1 + 1)
    ! Calculates the lifted parcel virtual temperature at nk1, the
    ! actual temperature, and the adiabatic liquid water content.

    ! Equivalent de TLIFT entre NK1 et ICB1+1 inclus

    ! Differences with convect4:
    ! - icbs1 is the first level above LCL (may differ from icb1)
    ! - in the iterations, used x(icbs1) instead x(icb1)
    ! - tvp1 is computed in only one time
    ! - icbs1: first level above Plcl1 (IMIN de TLIFT) in output
    ! - if icbs1=icb1, compute also tp1(icb1+1), tvp1(icb1+1) & clw1(icb1+1)

    use cv30_param_m, only: minorig, nl
    use cv_thermo_m, only: cl, clmcpv, cpd, cpv, eps, lv0, rrv
    USE dimphy, ONLY: klev, klon

    ! inputs:
    integer, intent(in):: nk1(klon), icb1(klon)
    real, intent(in):: t1(klon, klev)
    real, intent(in):: q1(klon, klev), qs1(klon, klev), gz1(klon, klev)
    real, intent(in):: p1(klon, klev)
    real, intent(in):: plcl1(klon)

    ! outputs:
    real tp1(klon, klev), tvp1(klon, klev), clw1(klon, klev)

    ! local variables:
    integer i, k
    integer icbs1(klon), icbsmax2
    real tg, qg, alv, s, ahg, tc, denom, es
    real ah0(klon), cpp(klon)
    real tnk(klon), qnk(klon), gznk(klon), ticb(klon), gzicb(klon)
    real qsicb(klon)
    real cpinv(klon)

    !-------------------------------------------------------------------

    !  Calculates the lifted parcel virtual temperature at nk1,
    !  the actual temperature, and the adiabatic
    !  liquid water content. The procedure is to solve the equation.
    ! cp*tp1+L*qp+phi=cp*tnk+L*qnk+gznk.

    do i=1, klon
       tnk(i)=t1(i, nk1(i))
       qnk(i)=q1(i, nk1(i))
       gznk(i)=gz1(i, nk1(i))
    end do

    ! *** Calculate certain parcel quantities, including static energy ***

    do i=1, klon
       ah0(i)=(cpd*(1.-qnk(i))+cl*qnk(i))*tnk(i) &
            +qnk(i)*(lv0-clmcpv*(tnk(i)-273.15))+gznk(i)
       cpp(i)=cpd*(1.-qnk(i))+qnk(i)*cpv
       cpinv(i)=1./cpp(i)
    end do

    ! *** Calculate lifted parcel quantities below cloud base ***

    do i=1, klon
       ! if icb1 is below LCL, start loop at ICB1+1:
       ! (icbs1 est le premier niveau au-dessus du LCL)
       icbs1(i)=MIN(max(icb1(i), 2), nl)
       if (plcl1(i) < p1(i, icbs1(i))) then
          icbs1(i)=MIN(icbs1(i)+1, nl)
       endif
    enddo

    do i=1, klon
       ticb(i)=t1(i, icbs1(i))
       gzicb(i)=gz1(i, icbs1(i))
       qsicb(i)=qs1(i, icbs1(i))
    enddo

    ! Re-compute icbsmax (icbsmax2):
    icbsmax2=2
    do i=1, klon
       icbsmax2=max(icbsmax2, icbs1(i))
    end do

    ! initialization outputs:

    do k=1, icbsmax2
       do i=1, klon
          tp1(i, k) = 0.0
          tvp1(i, k) = 0.0
          clw1(i, k) = 0.0
       enddo
    enddo

    ! tp1 and tvp1 below cloud base:

    do k=minorig, icbsmax2-1
       do i=1, klon
          tp1(i, k)=tnk(i)-(gz1(i, k)-gznk(i))*cpinv(i)
          tvp1(i, k)=tp1(i, k)*(1.+qnk(i)/eps-qnk(i))
       end do
    end do

    ! *** Find lifted parcel quantities above cloud base ***

    do i=1, klon
       tg=ticb(i)
       qg=qsicb(i)
       !debug alv=lv0-clmcpv*(ticb(i)-t0)
       alv=lv0-clmcpv*(ticb(i)-273.15)

       ! First iteration.

       s=cpd*(1.-qnk(i))+cl*qnk(i) &
            +alv*alv*qg/(rrv*ticb(i)*ticb(i))
       s=1./s

       ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i)
       tg=tg+s*(ah0(i)-ahg)

       !debug tc=tg-t0
       tc=tg-273.15
       denom=243.5+tc
       denom=MAX(denom, 1.0)

       es=6.112*exp(17.67*tc/denom)
       qg=eps*es/(p1(i, icbs1(i))-es*(1.-eps))

       ! Second iteration.

       ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i)
       tg=tg+s*(ah0(i)-ahg)

       !debug tc=tg-t0
       tc=tg-273.15
       denom=243.5+tc
       denom=MAX(denom, 1.0)

       es=6.112*exp(17.67*tc/denom)

       qg=eps*es/(p1(i, icbs1(i))-es*(1.-eps))

       alv=lv0-clmcpv*(ticb(i)-273.15)

       ! no approximation:
       tp1(i, icbs1(i))=(ah0(i)-gz1(i, icbs1(i))-alv*qg) &
            /(cpd+(cl-cpd)*qnk(i))

       clw1(i, icbs1(i))=qnk(i)-qg
       clw1(i, icbs1(i))=max(0.0, clw1(i, icbs1(i)))

       ! (qg utilise au lieu du vrai mixing ratio rg)
       tvp1(i, icbs1(i))=tp1(i, icbs1(i))*(1.+qg/eps-qnk(i))

    end do

    ! * icbs1 is the first level above the LCL:
    ! if plcl1<p1(icb1), then icbs1=icb1+1
    ! if plcl1>p1(icb1), then icbs1=icb1

    ! * the routine above computes tvp1 from minorig to icbs1 (included).

    ! * to compute buoybase (in cv30_trigger.F), both tvp1(icb1) and
    ! tvp1(icb1+1) must be known. This is the case if icbs1=icb1+1,
    ! but not if icbs1=icb1.

    ! * therefore, in the case icbs1=icb1, we compute tvp1 at level icb1+1
    ! (tvp1 at other levels will be computed in cv30_undilute2.F)

    do i=1, klon
       ticb(i)=t1(i, icb1(i)+1)
       gzicb(i)=gz1(i, icb1(i)+1)
       qsicb(i)=qs1(i, icb1(i)+1)
    enddo

    do i=1, klon
       tg=ticb(i)
       qg=qsicb(i)
       !debug alv=lv0-clmcpv*(ticb(i)-t0)
       alv=lv0-clmcpv*(ticb(i)-273.15)

       ! First iteration.

       s=cpd*(1.-qnk(i))+cl*qnk(i) &
            +alv*alv*qg/(rrv*ticb(i)*ticb(i))
       s=1./s

       ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i)
       tg=tg+s*(ah0(i)-ahg)

       !debug tc=tg-t0
       tc=tg-273.15
       denom=243.5+tc
       denom=MAX(denom, 1.0)

       es=6.112*exp(17.67*tc/denom)

       qg=eps*es/(p1(i, icb1(i)+1)-es*(1.-eps))

       ! Second iteration.

       ahg=cpd*tg+(cl-cpd)*qnk(i)*tg+alv*qg+gzicb(i)
       tg=tg+s*(ah0(i)-ahg)

       !debug tc=tg-t0
       tc=tg-273.15
       denom=243.5+tc
       denom=MAX(denom, 1.0)

       es=6.112*exp(17.67*tc/denom)

       qg=eps*es/(p1(i, icb1(i)+1)-es*(1.-eps))

       alv=lv0-clmcpv*(ticb(i)-273.15)

       ! no approximation:
       tp1(i, icb1(i)+1)=(ah0(i)-gz1(i, icb1(i)+1)-alv*qg) &
            /(cpd+(cl-cpd)*qnk(i))

       clw1(i, icb1(i)+1)=qnk(i)-qg
       clw1(i, icb1(i)+1)=max(0.0, clw1(i, icb1(i)+1))

       ! (qg utilise au lieu du vrai mixing ratio rg)
       tvp1(i, icb1(i)+1)=tp1(i, icb1(i)+1)*(1.+qg/eps-qnk(i)) !whole thing
    end do

  end SUBROUTINE cv30_undilute1

end module cv30_undilute1_m
