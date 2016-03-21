module cv30_unsat_m

  implicit none

contains

  SUBROUTINE cv30_unsat(nloc, ncum, nd, na, icb, inb, t, rr, rs, gz, u, v, p, &
       ph, th, tv, lv, cpn, ep, sigp, clw, m, ment, elij, delt, plcl, mp, rp, &
       up, vp, wt, water, evap, b)

    use cv30_param_m, only: nl, sigd
    use cvflag, only: cvflag_grav
    use cvthermo, only: cpd, ginv, grav

    ! inputs:
    integer, intent(in):: nloc, ncum, nd, na
    integer, intent(in):: icb(:), inb(:) ! (ncum)
    real t(nloc, nd), rr(nloc, nd), rs(nloc, nd)
    real gz(nloc, na)
    real u(nloc, nd), v(nloc, nd)
    real p(nloc, nd), ph(nloc, nd + 1)
    real th(nloc, na)
    real tv(nloc, na)
    real lv(nloc, na)
    real cpn(nloc, na)
    real ep(nloc, na), sigp(nloc, na), clw(nloc, na)
    real m(nloc, na), ment(nloc, na, na), elij(nloc, na, na)
    real, intent(in):: delt
    real plcl(nloc)

    ! outputs:
    real mp(nloc, na), rp(nloc, na), up(nloc, na), vp(nloc, na)
    real wt(nloc, na), water(nloc, na), evap(nloc, na)
    real b(:, :) ! (nloc, na)

    ! Local:
    integer i, j, il, num1
    real tinv, delti
    real awat, afac, afac1, afac2, bfac
    real pr1, pr2, sigt, b6, c6, revap, tevap, delth
    real amfac, amp2, xf, tf, fac2, ur, sru, fac, d, af, bf
    real ampmax
    real lvcp(nloc, na)
    real wdtrain(nloc)
    logical lwork(nloc)

    !------------------------------------------------------

    delti = 1. / delt
    tinv = 1. / 3.
    mp = 0.

    do i = 1, nl
       do il = 1, ncum
          mp(il, i) = 0.
          rp(il, i) = rr(il, i)
          up(il, i) = u(il, i)
          vp(il, i) = v(il, i)
          wt(il, i) = 0.001
          water(il, i) = 0.
          evap(il, i) = 0.
          b(il, i) = 0.
          lvcp(il, i) = lv(il, i) / cpn(il, i)
       enddo
    enddo

    ! check whether ep(inb) = 0, if so, skip precipitating 
    ! downdraft calculation 

    do il = 1, ncum
       lwork(il) = .TRUE.
       if (ep(il, inb(il)) < 0.0001) lwork(il) = .FALSE.
    enddo

    wdtrain(:ncum) = 0.

    downdraft_loop: DO i = nl + 1, 1, - 1
       num1 = 0

       do il = 1, ncum
          if (i <= inb(il) .and. lwork(il)) num1 = num1 + 1
       enddo

       if (num1 > 0) then
          ! integrate liquid water equation to find condensed water 
          ! and condensed water flux 

          ! calculate detrained precipitation 

          do il = 1, ncum
             if (i <= inb(il) .and. lwork(il)) then
                if (cvflag_grav) then
                   wdtrain(il) = grav * ep(il, i) * m(il, i) * clw(il, i)
                else
                   wdtrain(il) = 10. * ep(il, i) * m(il, i) * clw(il, i)
                endif
             endif
          enddo

          if (i > 1) then
             do j = 1, i - 1
                do il = 1, ncum
                   if (i <= inb(il) .and. lwork(il)) then
                      awat = elij(il, j, i) - (1. - ep(il, i)) * clw(il, i)
                      awat = amax1(awat, 0.)
                      if (cvflag_grav) then
                         wdtrain(il) = wdtrain(il) + grav * awat &
                              * ment(il, j, i)
                      else
                         wdtrain(il) = wdtrain(il) + 10. * awat * ment(il, j, i)
                      endif
                   endif
                enddo
             end do
          endif

          ! find rain water and evaporation using provisional 
          ! estimates of rp(i)and rp(i - 1) 

          do il = 1, ncum
             if (i <= inb(il) .and. lwork(il)) then
                wt(il, i) = 45.

                if (i < inb(il)) then
                   rp(il, i) = rp(il, i + 1) + (cpd * (t(il, i + 1) &
                        - t(il, i)) + gz(il, i + 1) - gz(il, i)) / lv(il, i)
                   rp(il, i) = 0.5 * (rp(il, i) + rr(il, i))
                endif
                rp(il, i) = amax1(rp(il, i), 0.)
                rp(il, i) = amin1(rp(il, i), rs(il, i))
                rp(il, inb(il)) = rr(il, inb(il))

                if (i == 1) then
                   afac = p(il, 1) * (rs(il, 1) - rp(il, 1)) &
                        / (1e4 + 2000. * p(il, 1) * rs(il, 1))
                else
                   rp(il, i - 1) = rp(il, i) + (cpd * (t(il, i) &
                        - t(il, i - 1)) + gz(il, i) - gz(il, i - 1)) / lv(il, i)
                   rp(il, i - 1) = 0.5 * (rp(il, i - 1) + rr(il, i - 1))
                   rp(il, i - 1) = amin1(rp(il, i - 1), rs(il, i - 1))
                   rp(il, i - 1) = amax1(rp(il, i - 1), 0.)
                   afac1 = p(il, i) * (rs(il, i) - rp(il, i)) &
                        / (1e4 + 2000. * p(il, i) * rs(il, i))
                   afac2 = p(il, i - 1) * (rs(il, i - 1) - rp(il, i - 1)) &
                        / (1e4 + 2000. * p(il, i - 1) * rs(il, i - 1))
                   afac = 0.5 * (afac1 + afac2)
                endif
                if (i == inb(il))afac = 0.
                afac = amax1(afac, 0.)
                bfac = 1. / (sigd * wt(il, i))

                ! prise en compte de la variation progressive de sigt dans
                ! les couches icb et icb - 1:
                ! pour plcl < ph(i + 1), pr1 = 0 & pr2 = 1
                ! pour plcl > ph(i), pr1 = 1 & pr2 = 0
                ! pour ph(i + 1) < plcl < ph(i), pr1 est la proportion a cheval
                ! sur le nuage, et pr2 est la proportion sous la base du
                ! nuage.
                pr1 = (plcl(il) - ph(il, i + 1)) / (ph(il, i) - ph(il, i + 1))
                pr1 = max(0., min(1., pr1))
                pr2 = (ph(il, i) - plcl(il)) / (ph(il, i) - ph(il, i + 1))
                pr2 = max(0., min(1., pr2))
                sigt = sigp(il, i) * pr1 + pr2

                b6 = bfac * 50. * sigd * (ph(il, i) - ph(il, i + 1)) * sigt &
                     * afac
                c6 = water(il, i + 1) + bfac * wdtrain(il) - 50. * sigd * bfac &
                     * (ph(il, i) - ph(il, i + 1)) * evap(il, i + 1)
                if (c6 > 0.) then
                   revap = 0.5 * (- b6 + sqrt(b6 * b6 + 4. * c6))
                   evap(il, i) = sigt * afac * revap
                   water(il, i) = revap * revap
                else
                   evap(il, i) = - evap(il, i + 1) + 0.02 * (wdtrain(il) &
                        + sigd * wt(il, i) * water(il, i + 1)) &
                        / (sigd * (ph(il, i) - ph(il, i + 1)))
                end if

                ! calculate precipitating downdraft mass flux under 
                ! hydrostatic approximation 

                if (i /= 1) then
                   tevap = amax1(0., evap(il, i))
                   delth = amax1(0.001, (th(il, i) - th(il, i - 1)))
                   if (cvflag_grav) then
                      mp(il, i) = 100. * ginv * lvcp(il, i) * sigd * tevap &
                           * (p(il, i - 1) - p(il, i)) / delth
                   else
                      mp(il, i) = 10. * lvcp(il, i) * sigd * tevap &
                           * (p(il, i - 1) - p(il, i)) / delth
                   endif

                   ! if hydrostatic assumption fails, 
                   ! solve cubic difference equation for downdraft theta 
                   ! and mass flux from two simultaneous differential eqns 

                   amfac = sigd * sigd * 70. * ph(il, i) &
                        * (p(il, i - 1) - p(il, i)) &
                        * (th(il, i) - th(il, i - 1)) / (tv(il, i) * th(il, i))
                   amp2 = abs(mp(il, i + 1) * mp(il, i + 1) - mp(il, i) &
                        * mp(il, i))
                   if (amp2 > (0.1 * amfac)) then
                      xf = 100. * sigd * sigd * sigd * (ph(il, i) &
                           - ph(il, i + 1))
                      tf = b(il, i) - 5. * (th(il, i) - th(il, i - 1)) &
                           * t(il, i) / (lvcp(il, i) * sigd * th(il, i))
                      af = xf * tf + mp(il, i + 1) * mp(il, i + 1) * tinv
                      bf = 2. * (tinv * mp(il, i + 1))**3 + tinv &
                           * mp(il, i + 1) * xf * tf + 50. * (p(il, i - 1) &
                           - p(il, i)) * xf * tevap
                      fac2 = 1.
                      if (bf < 0.)fac2 = - 1.
                      bf = abs(bf)
                      ur = 0.25 * bf * bf - af * af * af * tinv * tinv * tinv
                      if (ur >= 0.) then
                         sru = sqrt(ur)
                         fac = 1.
                         if ((0.5 * bf - sru) < 0.)fac = - 1.
                         mp(il, i) = mp(il, i + 1) * tinv &
                              + (0.5 * bf + sru)**tinv &
                              + fac * (abs(0.5 * bf - sru))**tinv
                      else
                         d = atan(2. * sqrt(- ur) / (bf + 1e-28))
                         if (fac2 < 0.)d = 3.14159 - d
                         mp(il, i) = mp(il, i + 1) * tinv + 2. &
                              * sqrt(af * tinv) * cos(d * tinv)
                      endif
                      mp(il, i) = amax1(0., mp(il, i))

                      if (cvflag_grav) then
                         ! Il y a vraisemblablement une erreur dans la
                         ! ligne 2 suivante: il faut diviser par (mp(il,
                         ! i) * sigd * grav) et non par (mp(il, i) + sigd
                         ! * 0.1).  Et il faut bien revoir les facteurs
                         ! 100.
                         b(il, i - 1) = b(il, i) + 100. * (p(il, i - 1) &
                              - p(il, i)) * tevap / (mp(il, i) + sigd * 0.1) &
                              - 10. * (th(il, i) - th(il, i - 1)) * t(il, i) &
                              / (lvcp(il, i) * sigd * th(il, i))
                      else
                         b(il, i - 1) = b(il, i) + 100. * (p(il, i - 1) &
                              - p(il, i)) * tevap / (mp(il, i) + sigd * 0.1) &
                              - 10. * (th(il, i) - th(il, i - 1)) * t(il, i) &
                              / (lvcp(il, i) * sigd * th(il, i))
                      endif
                      b(il, i - 1) = amax1(b(il, i - 1), 0.)
                   endif

                   ! limit magnitude of mp(i) to meet cfl condition 

                   ampmax = 2. * (ph(il, i) - ph(il, i + 1)) * delti
                   amp2 = 2. * (ph(il, i - 1) - ph(il, i)) * delti
                   ampmax = amin1(ampmax, amp2)
                   mp(il, i) = amin1(mp(il, i), ampmax)

                   ! force mp to decrease linearly to zero 
                   ! between cloud base and the surface 

                   if (p(il, i) > p(il, icb(il))) then
                      mp(il, i) = mp(il, icb(il)) * (p(il, 1) - p(il, i)) &
                           / (p(il, 1) - p(il, icb(il)))
                   endif
                endif ! i == 1

                ! find mixing ratio of precipitating downdraft 

                if (i /= inb(il)) then
                   rp(il, i) = rr(il, i)

                   if (mp(il, i) > mp(il, i + 1)) then
                      if (cvflag_grav) then
                         rp(il, i) = rp(il, i + 1) * mp(il, i + 1) + rr(il, i) &
                              * (mp(il, i) - mp(il, i + 1)) + 100. * ginv &
                              * 0.5 * sigd * (ph(il, i) - ph(il, i + 1)) &
                              * (evap(il, i + 1) + evap(il, i))
                      else
                         rp(il, i) = rp(il, i + 1) * mp(il, i + 1) + rr(il, i) &
                              * (mp(il, i) - mp(il, i + 1)) + 5. * sigd &
                              * (ph(il, i) - ph(il, i + 1)) &
                              * (evap(il, i + 1) + evap(il, i))
                      endif
                      rp(il, i) = rp(il, i) / mp(il, i)
                      up(il, i) = up(il, i + 1) * mp(il, i + 1) + u(il, i) &
                           * (mp(il, i) - mp(il, i + 1))
                      up(il, i) = up(il, i) / mp(il, i)
                      vp(il, i) = vp(il, i + 1) * mp(il, i + 1) + v(il, i) &
                           * (mp(il, i) - mp(il, i + 1))
                      vp(il, i) = vp(il, i) / mp(il, i)
                   else
                      if (mp(il, i + 1) > 1e-16) then
                         if (cvflag_grav) then
                            rp(il, i) = rp(il, i + 1) &
                                 + 100. * ginv * 0.5 * sigd * (ph(il, i) &
                                 - ph(il, i + 1)) &
                                 * (evap(il, i + 1) + evap(il, i)) &
                                 / mp(il, i + 1)
                         else
                            rp(il, i) = rp(il, i + 1) &
                                 + 5. * sigd * (ph(il, i) - ph(il, i + 1)) &
                                 * (evap(il, i + 1) + evap(il, i)) &
                                 / mp(il, i + 1)
                         endif
                         up(il, i) = up(il, i + 1)
                         vp(il, i) = vp(il, i + 1)
                      endif
                   endif
                   rp(il, i) = amin1(rp(il, i), rs(il, i))
                   rp(il, i) = amax1(rp(il, i), 0.)
                endif
             endif
          end do
       end if
    end DO downdraft_loop

  end SUBROUTINE cv30_unsat

end module cv30_unsat_m
