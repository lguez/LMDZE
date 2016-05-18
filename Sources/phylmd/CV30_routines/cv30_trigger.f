module cv30_trigger_m

  implicit none

contains

  SUBROUTINE cv30_trigger(icb1, plcl1, p1, th1, tv1, tvp1, pbase1, buoybase1, &
       iflag1, sig1, w01)

    ! Triggering:
    ! - computes the cloud base
    ! - triggering (crude in this version)
    ! - relaxation of sig1 and w01 when no convection

    ! Caution 1: if no convection, we set iflag1 = 4

    ! Caution 2: at this stage, tvp1 (and thus buoy) are known up
    ! through icb1 only!  -> the buoyancy below cloud base not (yet)
    ! set to the cloud base buoyancy

    use cv30_param_m, only: alpha, beta, dtcrit, nl
    USE dimphy, ONLY: klev, klon

    integer, intent(in):: icb1(klon)
    real, intent(in):: plcl1(klon), p1(klon, klev)
    real, intent(in):: th1(klon, klev), tv1(klon, klev), tvp1(klon, klev)

    real, intent(out):: pbase1(klon), buoybase1(klon)

    integer, intent(inout):: iflag1(klon)
    real, intent(inout):: sig1(klon, klev), w01(klon, klev)

    ! Local:
    real, parameter:: dttrig = 5. ! (loose) condition for triggering
    real, parameter:: dpbase = - 40. ! definition cloud base (400 m above LCL)
    integer i, k
    real tvpbase, tvbase

    !---------------------------------------------------------------------

    ! Set cloud base buoyancy at plcl1 + dpbase level buoyancy:
    do i = 1, klon
       pbase1(i) = plcl1(i) + dpbase
       tvpbase = tvp1(i, icb1(i)) * (pbase1(i) - p1(i, icb1(i) + 1)) &
            /(p1(i, icb1(i)) - p1(i, icb1(i) + 1)) &
            + tvp1(i, icb1(i) + 1) * (p1(i, icb1(i)) - pbase1(i)) &
            /(p1(i, icb1(i)) - p1(i, icb1(i) + 1))
       tvbase = tv1(i, icb1(i)) * (pbase1(i) - p1(i, icb1(i) + 1)) &
            /(p1(i, icb1(i)) - p1(i, icb1(i) + 1)) &
            + tv1(i, icb1(i) + 1) * (p1(i, icb1(i)) - pbase1(i)) &
            /(p1(i, icb1(i)) - p1(i, icb1(i) + 1))
       buoybase1(i) = tvpbase - tvbase
    end do

    ! Make sure that column is dry adiabatic between the surface and
    ! cloud base, and that lifted air is positively buoyant at cloud
    ! base.  If not, return to calling program after resetting sig1(i)
    ! and w01(i).
    do k = 1, nl
       do i = 1, klon
          if (buoybase1(i) < dtcrit .or. th1(i, icb1(i) - 1) - dttrig &
               > th1(i, 1)) then
             sig1(i, k) = MAX(beta * sig1(i, k) - 2. * alpha &
                  * buoybase1(i)**2, 0.)
             w01(i, k) = beta * w01(i, k)
             iflag1(i) = 4
          endif
       end do
    end do

  end SUBROUTINE cv30_trigger

end module cv30_trigger_m
