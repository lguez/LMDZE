module cv30_mixing_m

  implicit none

contains

  SUBROUTINE cv30_mixing(nloc, ncum, nd, na, icb, nk, inb, t, rr, rs, u, v, h, &
       lv, hp, ep, clw, m, sig, ment, qent, uent, vent, nent, sij, elij, &
       ments, qents)
    use cv30_param_m
    use cvthermo

    !---------------------------------------------------------------------
    ! a faire:
    ! - changer rr(il, 1) -> qnk(il)
    ! - vectorisation de la partie normalisation des flux (do 789)
    !---------------------------------------------------------------------

    ! inputs:
    integer, intent(in):: ncum, nd, na, nloc
    integer, intent(in):: icb(nloc), inb(nloc), nk(nloc)
    real sig(nloc, nd)
    real t(nloc, nd), rr(nloc, nd), rs(nloc, nd)
    real u(nloc, nd), v(nloc, nd)
    real lv(nloc, na), h(nloc, na), hp(nloc, na)
    real ep(nloc, na), clw(nloc, na)
    real m(nloc, na) ! input of convect3

    ! outputs:
    real ment(nloc, na, na), qent(nloc, na, na)
    real uent(nloc, na, na), vent(nloc, na, na)
    real sij(nloc, na, na), elij(nloc, na, na)
    real ments(nloc, nd, nd), qents(nloc, nd, nd)
    integer nent(nloc, nd)

    ! local variables:
    integer i, j, k, il, im, jm
    integer num1, num2
    real rti, bf2, anum, denom, dei, altem, cwat, stemp, qp
    real alt, smid, sjmin, sjmax, delp, delm
    real asij(nloc), smax(nloc), scrit(nloc)
    real asum(nloc, nd), bsum(nloc, nd), csum(nloc, nd)
    real wgh
    real zm(nloc, na)
    logical lwork(nloc)

    ! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS

    do j=1, nl
       do i=1, ncum
          nent(i, j)=0
          ! in convect3, m is computed in cv30_closure
          ! ori m(i, 1)=0.0
       end do
    end do

    do j=1, nl
       do k=1, nl
          do i=1, ncum
             qent(i, k, j)=rr(i, j)
             uent(i, k, j)=u(i, j)
             vent(i, k, j)=v(i, j)
             elij(i, k, j)=0.0
             !ym ment(i, k, j)=0.0
             !ym sij(i, k, j)=0.0
          end do
       end do
    end do

    !ym
    ment(1:ncum, 1:nd, 1:nd)=0.0
    sij(1:ncum, 1:nd, 1:nd)=0.0

    zm(:, :)=0.

    ! --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
    ! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
    ! --- FRACTION (sij)

    do i=minorig+1, nl

       do j=minorig, nl
          do il=1, ncum
             if((i >= icb(il)).and.(i <= inb(il)).and. &
                  (j >= (icb(il)-1)).and.(j <= inb(il)))then

                rti=rr(il, 1)-ep(il, i)*clw(il, i)
                bf2=1.+lv(il, j)*lv(il, j)*rs(il, j)/(rrv*t(il, j)*t(il, j)*cpd)
                anum=h(il, j)-hp(il, i)+(cpv-cpd)*t(il, j)*(rti-rr(il, j))
                denom=h(il, i)-hp(il, i)+(cpd-cpv)*(rr(il, i)-rti)*t(il, j)
                dei=denom
                if(abs(dei) < 0.01)dei=0.01
                sij(il, i, j)=anum/dei
                sij(il, i, i)=1.0
                altem=sij(il, i, j)*rr(il, i)+(1.-sij(il, i, j))*rti-rs(il, j)
                altem=altem/bf2
                cwat=clw(il, j)*(1.-ep(il, j))
                stemp=sij(il, i, j)
                if((stemp < 0.0.or.stemp > 1.0.or.altem > cwat) &
                     .and.j > i)then
                   anum=anum-lv(il, j)*(rti-rs(il, j)-cwat*bf2)
                   denom=denom+lv(il, j)*(rr(il, i)-rti)
                   if(abs(denom) < 0.01)denom=0.01
                   sij(il, i, j)=anum/denom
                   altem=sij(il, i, j)*rr(il, i)+(1.-sij(il, i, j))*rti-rs(il, j)
                   altem=altem-(bf2-1.)*cwat
                end if
                if(sij(il, i, j) > 0.0.and.sij(il, i, j) < 0.95)then
                   qent(il, i, j)=sij(il, i, j)*rr(il, i)+(1.-sij(il, i, j))*rti
                   uent(il, i, j)=sij(il, i, j)*u(il, i)+(1.-sij(il, i, j))*u(il, nk(il))
                   vent(il, i, j)=sij(il, i, j)*v(il, i)+(1.-sij(il, i, j))*v(il, nk(il))
                   elij(il, i, j)=altem
                   elij(il, i, j)=amax1(0.0, elij(il, i, j))
                   ment(il, i, j)=m(il, i)/(1.-sij(il, i, j))
                   nent(il, i)=nent(il, i)+1
                end if
                sij(il, i, j)=amax1(0.0, sij(il, i, j))
                sij(il, i, j)=amin1(1.0, sij(il, i, j))
             endif ! new
          end do
       end do

       ! *** if no air can entrain at level i assume that updraft detrains ***
       ! *** at that level and calculate detrained air flux and properties ***

       do il=1, ncum
          if ((i >= icb(il)).and.(i <= inb(il)).and.(nent(il, i) == 0)) then
             !@ if(nent(il, i) == 0)then
             ment(il, i, i)=m(il, i)
             qent(il, i, i)=rr(il, nk(il))-ep(il, i)*clw(il, i)
             uent(il, i, i)=u(il, nk(il))
             vent(il, i, i)=v(il, nk(il))
             elij(il, i, i)=clw(il, i)
             !MAF sij(il, i, i)=1.0
             sij(il, i, i)=0.0
          end if
       end do
    end do

    ! --- NORMALIZE ENTRAINED AIR MASS FLUXES
    ! --- TO REPRESENT EQUAL PROBABILITIES OF MIXING

    asum = 0.
    csum = 0.

    do il=1, ncum
       lwork(il) = .FALSE.
    enddo

    DO i=minorig+1, nl

       num1=0
       do il=1, ncum
          if (i >= icb(il) .and. i <= inb(il)) num1=num1+1
       enddo
       if (num1 <= 0) cycle

       do il=1, ncum
          if (i >= icb(il) .and. i <= inb(il)) then
             lwork(il)=(nent(il, i) /= 0)
             qp=rr(il, 1)-ep(il, i)*clw(il, i)
             anum=h(il, i)-hp(il, i)-lv(il, i)*(qp-rs(il, i)) &
                  +(cpv-cpd)*t(il, i)*(qp-rr(il, i))
             denom=h(il, i)-hp(il, i)+lv(il, i)*(rr(il, i)-qp) &
                  +(cpd-cpv)*t(il, i)*(rr(il, i)-qp)
             if(abs(denom) < 0.01)denom=0.01
             scrit(il)=anum/denom
             alt=qp-rs(il, i)+scrit(il)*(rr(il, i)-qp)
             if(scrit(il) <= 0.0.or.alt <= 0.0)scrit(il)=1.0
             smax(il)=0.0
             asij(il)=0.0
          endif
       end do

       do j=nl, minorig, -1

          num2=0
          do il=1, ncum
             if (i >= icb(il) .and. i <= inb(il) .and. &
                  j >= (icb(il)-1) .and. j <= inb(il) &
                  .and. lwork(il)) num2=num2+1
          enddo
          if (num2 <= 0) cycle

          do il=1, ncum
             if (i >= icb(il) .and. i <= inb(il) .and. &
                  j >= (icb(il)-1) .and. j <= inb(il) &
                  .and. lwork(il)) then

                if(sij(il, i, j) > 1.0e-16.and.sij(il, i, j) < 0.95)then
                   wgh=1.0
                   if(j > i)then
                      sjmax=amax1(sij(il, i, j+1), smax(il))
                      sjmax=amin1(sjmax, scrit(il))
                      smax(il)=amax1(sij(il, i, j), smax(il))
                      sjmin=amax1(sij(il, i, j-1), smax(il))
                      sjmin=amin1(sjmin, scrit(il))
                      if(sij(il, i, j) < (smax(il)-1.0e-16))wgh=0.0
                      smid=amin1(sij(il, i, j), scrit(il))
                   else
                      sjmax=amax1(sij(il, i, j+1), scrit(il))
                      smid=amax1(sij(il, i, j), scrit(il))
                      sjmin=0.0
                      if(j > 1)sjmin=sij(il, i, j-1)
                      sjmin=amax1(sjmin, scrit(il))
                   endif
                   delp=abs(sjmax-smid)
                   delm=abs(sjmin-smid)
                   asij(il)=asij(il)+wgh*(delp+delm)
                   ment(il, i, j)=ment(il, i, j)*(delp+delm)*wgh
                endif
             endif
          end do

       end do

       do il=1, ncum
          if (i >= icb(il).and.i <= inb(il).and.lwork(il)) then
             asij(il)=amax1(1.0e-16, asij(il))
             asij(il)=1.0/asij(il)
             asum(il, i)=0.0
             bsum(il, i)=0.0
             csum(il, i)=0.0
          endif
       enddo

       do j=minorig, nl
          do il=1, ncum
             if (i >= icb(il) .and. i <= inb(il) .and. lwork(il) &
                  .and. j >= (icb(il)-1) .and. j <= inb(il)) then
                ment(il, i, j)=ment(il, i, j)*asij(il)
             endif
          enddo
       end do

       do j=minorig, nl
          do il=1, ncum
             if (i >= icb(il) .and. i <= inb(il) .and. lwork(il) &
                  .and. j >= (icb(il)-1) .and. j <= inb(il)) then
                asum(il, i)=asum(il, i)+ment(il, i, j)
                ment(il, i, j)=ment(il, i, j)*sig(il, j)
                bsum(il, i)=bsum(il, i)+ment(il, i, j)
             endif
          enddo
       end do

       do il=1, ncum
          if (i >= icb(il).and.i <= inb(il).and.lwork(il)) then
             bsum(il, i)=amax1(bsum(il, i), 1.0e-16)
             bsum(il, i)=1.0/bsum(il, i)
          endif
       enddo

       do j=minorig, nl
          do il=1, ncum
             if (i >= icb(il) .and. i <= inb(il) .and. lwork(il) &
                  .and. j >= (icb(il)-1) .and. j <= inb(il)) then
                ment(il, i, j)=ment(il, i, j)*asum(il, i)*bsum(il, i)
             endif
          enddo
       end do

       do j=minorig, nl
          do il=1, ncum
             if (i >= icb(il) .and. i <= inb(il) .and. lwork(il) &
                  .and. j >= (icb(il)-1) .and. j <= inb(il)) then
                csum(il, i)=csum(il, i)+ment(il, i, j)
             endif
          enddo
       end do

       do il=1, ncum
          if (i >= icb(il) .and. i <= inb(il) .and. lwork(il) &
               .and. csum(il, i) < m(il, i)) then
             nent(il, i)=0
             ment(il, i, i)=m(il, i)
             qent(il, i, i)=rr(il, 1)-ep(il, i)*clw(il, i)
             uent(il, i, i)=u(il, nk(il))
             vent(il, i, i)=v(il, nk(il))
             elij(il, i, i)=clw(il, i)
             !MAF sij(il, i, i)=1.0
             sij(il, i, i)=0.0
          endif
       enddo ! il

    end DO

    ! MAF: renormalisation de MENT
    do jm=1, nd
       do im=1, nd
          do il=1, ncum
             zm(il, im)=zm(il, im)+(1.-sij(il, im, jm))*ment(il, im, jm)
          end do
       end do
    end do

    do jm=1, nd
       do im=1, nd
          do il=1, ncum
             if(zm(il, im) /= 0.) then
                ment(il, im, jm)=ment(il, im, jm)*m(il, im)/zm(il, im)
             endif
          end do
       end do
    end do

    do jm=1, nd
       do im=1, nd
          do il=1, ncum
             qents(il, im, jm)=qent(il, im, jm)
             ments(il, im, jm)=ment(il, im, jm)
          end do
       enddo
    enddo

  end SUBROUTINE cv30_mixing

end module cv30_mixing_m
