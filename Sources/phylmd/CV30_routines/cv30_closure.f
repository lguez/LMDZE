module cv30_closure_m

  implicit none

contains

  SUBROUTINE cv30_closure(icb, inb, pbase, p, ph, tv, buoy, sig, w0, cape, m)

    ! CLOSURE
    ! Vectorization: S. Bony

    use cv30_param_m, only: alpha, beta, dtcrit, minorig, nl
    use cv_thermo_m, only: rrd
    USE dimphy, ONLY: klev, klon

    ! input:
    integer, intent(in):: icb(:), inb(:) ! (ncum)
    real pbase(klon)
    real p(:, :) ! (klon, klev)
    real, intent(in):: ph(:, :)  ! (ncum, klev + 1)
    real tv(klon, klev), buoy(klon, klev)

    ! input/output:
    real sig(klon, klev), w0(klon, klev)

    ! output:
    real cape(klon)
    real m(klon, klev)

    ! Local:
    integer ncum
    integer i, j, k, icbmax
    real deltap, fac, w, amu
    real dtmin(klon, klev), sigold(klon, klev)

    !-------------------------------------------------------

    ncum = size(icb)

    ! Initialization

    do k=1, nl
       do i=1, ncum
          m(i, k)=0.0
       enddo
    enddo

    ! Reset sig(i) and w0(i) for i>inb and i<icb

    ! update sig and w0 above LNB:

    do k=1, nl-1
       do i=1, ncum
          if ((inb(i) < (nl-1)).and.(k >= (inb(i) + 1)))then
             sig(i, k)=beta*sig(i, k) &
                  + 2.*alpha*buoy(i, inb(i))*ABS(buoy(i, inb(i)))
             sig(i, k)=AMAX1(sig(i, k), 0.0)
             w0(i, k)=beta*w0(i, k)
          endif
       end do
    end do

    ! compute icbmax:

    icbmax=2
    do i=1, ncum
       icbmax=MAX(icbmax, icb(i))
    end do

    ! update sig and w0 below cloud base:

    do k=1, icbmax
       do i=1, ncum
          if (k <= icb(i))then
             sig(i, k)=beta*sig(i, k)-2.*alpha*buoy(i, icb(i))*buoy(i, icb(i))
             sig(i, k)=amax1(sig(i, k), 0.0)
             w0(i, k)=beta*w0(i, k)
          endif
       end do
    end do

    ! Reset fractional areas of updrafts and w0 at initial time
    ! and after 10 time steps of no convection

    do k=1, nl-1
       do i=1, ncum
          if (sig(i, klev) < 1.5.or.sig(i, klev) > 12.0)then
             sig(i, k)=0.0
             w0(i, k)=0.0
          endif
       end do
    end do

    ! Calculate convective available potential energy (cape),
    ! vertical velocity (w), fractional area covered by
    ! undilute updraft (sig), and updraft mass flux (m)

    do i=1, ncum
       cape(i)=0.0
    end do

    ! compute dtmin (minimum buoyancy between ICB and given level k):

    do i=1, ncum
       do k=1, nl
          dtmin(i, k)=100.0
       enddo
    enddo

    do i=1, ncum
       do k=1, nl
          do j=minorig, nl
             if ((k >= (icb(i) + 1)).and.(k <= inb(i)).and. &
                  (j >= icb(i)).and.(j <= (k-1)))then
                dtmin(i, k)=AMIN1(dtmin(i, k), buoy(i, j))
             endif
          end do
       end do
    end do

    ! The interval on which cape is computed starts at pbase:

    do k=1, nl
       do i=1, ncum
          if ((k >= (icb(i) + 1)).and.(k <= inb(i))) then
             deltap = MIN(pbase(i), ph(i, k-1))-MIN(pbase(i), ph(i, k))
             cape(i)=cape(i) + rrd*buoy(i, k-1)*deltap/p(i, k-1)
             cape(i)=AMAX1(0.0, cape(i))
             sigold(i, k)=sig(i, k)

             sig(i, k)=beta*sig(i, k) + alpha*dtmin(i, k)*ABS(dtmin(i, k))
             sig(i, k)=amax1(sig(i, k), 0.0)
             sig(i, k)=amin1(sig(i, k), 0.01)
             fac=AMIN1(((dtcrit-dtmin(i, k))/dtcrit), 1.0)
             w=(1.-beta)*fac*SQRT(cape(i)) + beta*w0(i, k)
             amu=0.5*(sig(i, k) + sigold(i, k))*w
             m(i, k)=amu*0.007*p(i, k)*(ph(i, k)-ph(i, k + 1))/tv(i, k)
             w0(i, k)=w
          endif

       end do
    end do

    do i=1, ncum
       w0(i, icb(i))=0.5*w0(i, icb(i) + 1)
       m(i, icb(i))=0.5*m(i, icb(i) + 1) &
            *(ph(i, icb(i))-ph(i, icb(i) + 1)) &
            /(ph(i, icb(i) + 1)-ph(i, icb(i) + 2))
       sig(i, icb(i))=sig(i, icb(i) + 1)
       sig(i, icb(i)-1)=sig(i, icb(i))
    end do

  end SUBROUTINE cv30_closure

end module cv30_closure_m
