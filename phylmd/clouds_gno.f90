module CLOUDS_GNO_m

  IMPLICIT NONE

contains

  SUBROUTINE CLOUDS_GNO(Q_SERI, QSAT, QCONDC, PTCONV, RATQSC, RNEBCON0)

    ! From LMDZ4/libf/phylmd/clouds_gno.F, version 1.2, 2004/11/09 16:55:40

    ! Libraries:
    use numer_rec_95, only: nr_erf
    use jumble, only: pi, SQRT2

    use dimphy, only: klon, klev

    REAL, intent(in):: q_seri(:, :) ! (klon, llm)
    ! domain-averaged mixing ratio of total water

    REAL, intent(in):: QSAT(:, :) ! (klon, llm)
    ! mean saturation humidity mixing ratio within the gridbox

    REAL, intent(in):: QCONDC(:, :) ! (klon, llm)
    ! mixing ratio of condensed water within clouds associated
    ! with SUBGRID-SCALE condensation processes (here, it is
    ! predicted by the convection scheme)

    LOGICAL, intent(out):: PTCONV(:, :) ! (klon, llm) point convectif

    REAL, intent(out):: RATQSC(:, :) ! (klon, llm)
    ! largeur normalisee de la distribution

    REAL, intent(out):: RNEBCON0(:, :) ! (klon, llm) fraction nuageuse

    ! Local:

    ! Parameters controlling the iteration:

    INTEGER, PARAMETER:: nmax = 10
    ! maximum nb of iterations (hopefully never reached)

    REAL, parameter:: my_eps = 0.02 ! accuracy of the numerical resolution
    REAL, parameter:: vmax0 = 2.

    real vmax(klon)
    ! v-value above which we use an asymptotic expression for ERF(v)

    REAL, parameter:: min_mu = 1.e-12, min_Q = 1.e-12
    INTEGER i, K, n
    REAL mu(klon), delta(klon), beta(klon)
    real zu2(klon), zv2(klon)
    REAL xx(klon), aux(klon), coeff(klon), my_block(klon)
    REAL dist(klon), fprime(klon), det(klon)
    REAL u(klon), v(klon), erfcu(klon), erfcv(klon)
    REAL xx1(klon), xx2(klon)
    real sqrtpi, zx1, zx2, exdel
    LOGICAL lconv(klon) ! le calcul a converge (entre autres si qcondc < min_q)

    !--------------------------------------------------------------

    rnebcon0 = 0.
    sqrtpi = sqrt(pi)
    ptconv = .false.
    ratqsc = 0.

    loop_vertical: DO K = 1, klev
       do i = 1, klon
          mu(i) = MAX(Q_SERI(i, K), min_mu)
          delta(i) = log(mu(i)/MAX(QSAT(i, K), min_mu))
       enddo

       ! There is no subgrid-scale condensation; the scheme becomes
       ! equivalent to an "all-or-nothing" large-scale condensation
       ! scheme.

       ! Some condensation is produced at the subgrid-scale
       !
       ! PDF = generalized log-normal distribution (GNO)
       ! (k<0 because a lower bound is considered for the PDF)
       !
       ! -> Determine x (the parameter k of the GNO PDF) such that the
       ! contribution of subgrid-scale processes to the in-cloud water
       ! content is equal to QCONDC(K) (equations (13), (14), (15) +
       ! Appendix B of the paper)
       !
       ! Here, an iterative method is used for this purpose (other
       ! numerical methods might be more efficient)
       !
       ! NB: the "error function" is called ERF (ERF in double
       ! precision)

       ! On commence par eliminer les cas pour lesquels on n'a pas
       ! suffisamment d'eau nuageuse.

       do i = 1, klon
          IF ( QCONDC(i, K) .lt. min_Q ) THEN
             ptconv(i, k) = .false.
             ratqsc(i, k) = 0.
             lconv(i) = .true.
          ELSE
             lconv(i) = .FALSE.
             vmax(i) = vmax0

             beta(i) = QCONDC(i, K)/mu(i) + EXP( -MIN(0., delta(i)) )

             ! roots of equation v > vmax:

             det(i) = delta(i) + vmax(i)**2.
             if (det(i).LE.0.) vmax(i) = vmax0 + 1.
             det(i) = delta(i) + vmax(i)**2.

             if (det(i).LE.0.) then
                xx(i) = -0.0001
             else
                zx1 = -sqrt2*vmax(i)
                zx2 = SQRT(1.+delta(i)/(vmax(i)**2.))
                xx1(i) = zx1*(1.-zx2)
                xx2(i) = zx1*(1.+zx2)
                xx(i) = 1.01 * xx1(i)
                if ( xx1(i) .GE. 0. ) xx(i) = 0.5*xx2(i)
             endif
             if (delta(i).LT.0.) xx(i) = -0.5*SQRT(log(2.))
          ENDIF
       enddo

       ! It\'erations pour trouver la solution :
       loop_n: DO n = 1, nmax
          loop_horizontal: do i = 1, klon
             test_lconv: if (.not.lconv(i)) then
                u(i) = delta(i)/(xx(i)*sqrt2) + xx(i)/(2.*sqrt2)
                v(i) = delta(i)/(xx(i)*sqrt2) - xx(i)/(2.*sqrt2)

                IF ( v(i) .GT. vmax(i) ) THEN
                   IF ( ABS(u(i)) .GT. vmax(i) .AND. delta(i) .LT. 0. ) THEN
                      ! use asymptotic expression of erf for u and v large:
                      ! ( -> analytic solution for xx )
                      exdel = beta(i)*EXP(delta(i))
                      aux(i) = 2.*delta(i)*(1.-exdel) /(1.+exdel)
                      if (aux(i).lt.0.) then
                         aux(i) = 0.
                      endif
                      xx(i) = -SQRT(aux(i))
                      my_block(i) = EXP(-v(i)*v(i)) / v(i) / sqrtpi
                      dist(i) = 0.
                      fprime(i) = 1.
                   ELSE
                      ! erfv -> 1., use an asymptotic expression of
                      ! erfv for v large:

                      erfcu(i) = 1.-NR_ERF(u(i))
                      ! Attention : ajout d'un seuil pour l'exponentielle
                      aux(i) = sqrtpi*erfcu(i)*EXP(min(v(i)*v(i), 80.))
                      coeff(i) = 1. - 1./2./(v(i)**2.) + 3./4./(v(i)**4.)
                      my_block(i) = coeff(i) * EXP(-v(i)*v(i)) / v(i) / sqrtpi
                      dist(i) = v(i) * aux(i) / coeff(i) - beta(i)
                      fprime(i) = 2. / xx(i) * (v(i)**2.) &
                           * ( coeff(i)*EXP(-delta(i)) - u(i) * aux(i) ) &
                           / coeff(i) / coeff(i)
                   ENDIF
                ELSE
                   ! general case:

                   erfcu(i) = 1.-NR_ERF(u(i))
                   erfcv(i) = 1.-NR_ERF(v(i))
                   my_block(i) = erfcv(i)
                   dist(i) = erfcu(i) / erfcv(i) - beta(i)
                   zu2(i) = u(i)*u(i)
                   zv2(i) = v(i)*v(i)
                   if(zu2(i).gt.20..or. zv2(i).gt.20.) then
                      zu2(i) = 20.
                      zv2(i) = 20.
                      fprime(i) = 0.
                   else
                      fprime(i) = 2. /sqrtpi /xx(i) /erfcv(i)**2. &
                           * ( erfcv(i)*v(i)*EXP(-zu2(i)) &
                           - erfcu(i)*u(i)*EXP(-zv2(i)) )
                   endif
                ENDIF

                ! test numerical convergence:
                if ( ABS(dist(i)/beta(i)) .LT. my_eps ) then
                   ptconv(i, K) = .TRUE.
                   lconv(i) = .true.
                   ! borne pour l'exponentielle
                   ratqsc(i, k) = min(2.*(v(i)-u(i))**2, 20.)
                   ratqsc(i, k) = sqrt(exp(ratqsc(i, k))-1.)
                   RNEBCON0(i, K) = 0.5 * my_block(i)
                else
                   xx(i) = xx(i) - dist(i)/fprime(i)
                endif
             endif test_lconv
          enddo loop_horizontal
       ENDDO loop_n
    end DO loop_vertical

  END SUBROUTINE CLOUDS_GNO

end module CLOUDS_GNO_m
