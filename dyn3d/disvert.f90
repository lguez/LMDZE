module disvert_m

  implicit none

  private hybrid, funcd, compute_ab

  real, save, allocatable, protected:: ap(:) ! (llm+1) in Pa
  real, save, allocatable, protected:: bp(:) ! (llm+1)

  REAL, allocatable, protected:: s(:) ! (llm+1)
  ! "s(l)" is the atmospheric hybrid sigma-pressure coordinate at
  ! half-level, between layers "l" and "l-1"

  real, save, allocatable, protected:: presnivs(:) ! (llm)
  ! approximate full level pressure for a reference surface pressure, in Pa

  real, parameter:: preff = 101325. ! in Pa
  real, private:: y, ya ! for the hybrid function
  real, parameter, private:: pah = 5e4 ! in Pa

contains

  SUBROUTINE disvert

    ! From dyn3d/disvert.F, version 1.1.1.1, 2004/05/19 12:53:05
    ! Author: P. Le Van

    ! This procedure sets the vertical grid. It defines the host
    ! variables "ap", "bp", "presnivs".

    ! Libraries:
    use jumble, only: read_column, new_unit, pi, assert

    use dimensions, only: llm
    use unit_nml_m, only: unit_nml

    ! Local:

    real ds(llm)
    ! "ds(l)" : \'epaisseur de la couche "l" dans la coordonn\'ee "s"

    INTEGER l, unit
    REAL x(llm)
    real:: dsigmin = 1.
    real zz(llm) ! in km

    character(len=20):: vert_sampling = "tropo"
    ! Allowed values: "tropo", "strato_custom", "strato",
    ! "read_hybrid", "read_pressure".

    ! These variables are used only in the case vert_sampling ==
    ! "strato_custom", and all are in km:
    real:: vert_scale = 7. ! scale height
    real:: vert_dzmin = 0.017 ! width of first layer
    real:: vert_dzlow = 1. ! dz in the low atmosphere
    real:: vert_z0low = 8.7 ! height at which resolution reaches dzlow
    real:: vert_dzmid = 2. ! dz in the mid atmosphere
    real:: vert_z0mid = 70. ! height at which resolution reaches dzmid
    real:: vert_h_mid = 20. ! width of the transition
    real:: vert_dzhig = 11. ! dz in the high atmosphere
    real:: vert_z0hig = 75. ! height at which resolution reaches dz
    real:: vert_h_hig = 20. ! width of the transition

    real, allocatable:: p(:) ! (2:llm or llm + 1) pressure (in hPa)

    namelist /disvert_nml/vert_sampling, vert_scale, vert_dzmin, vert_dzlow, &
         vert_z0low, vert_dzmid, vert_z0mid, vert_h_mid, vert_dzhig, &
         vert_z0hig, vert_h_hig, dsigmin

    !-----------------------------------------------------------------------

    print *, "Call sequence information: disvert"

    allocate(ap(llm+1))
    allocate(bp(llm+1))
    allocate(s(llm+1))
    allocate(presnivs(llm))

    write(unit=*, nml=disvert_nml)
    print *, "Enter namelist 'disvert_nml'."
    read(unit=*, nml=disvert_nml)
    write(unit_nml, nml=disvert_nml)

    s(1) = 1.
    s(llm+1) = 0.

    select case (vert_sampling)
    case ("tropo")
       ! with llm = 19 and dsigmin = 1 for CMIP 3

       forall (l = 1: llm) ds(l) &
            = dsigmin + 7. * SIN(pi * (REAL(l) - 0.5) / real(llm + 1))**2
       ds = ds / sum(ds)

       DO l = llm, 2, -1
          s(l) = s(l+1) + ds(l)
       ENDDO

       call compute_ab
    case ("strato")
       ! with llm = 39 and dsigmin = 0.3 for CMIP5

       forall (l = 1: llm) x(l) = pi * (l - 0.5) / (llm + 1)

       ds = (dsigmin + 7. * SIN(x)**2) * (1. - tanh(2 * x / pi - 1.))**2 / 4.
       ds = ds / sum(ds)

       DO l = llm, 2, -1
          s(l) = s(l+1) + ds(l)
       ENDDO

       call compute_ab
    case ("strato_custom")
       ! with llm = 79 for CMIP 6

       zz(1) = 0.
       DO l = 1, llm - 1
          zz(l + 1) = zz(l) + vert_dzmin + vert_dzlow &
               * TANH(zz(l) / vert_z0low) + (vert_dzmid - vert_dzlow) &
               * (TANH((zz(l) - vert_z0mid) / vert_h_mid) &
               - TANH(- vert_z0mid / vert_h_mid)) &
               + (vert_dzhig - vert_dzmid - vert_dzlow) &
               * (TANH((zz(l) - vert_z0hig) / vert_h_hig) &
               - TANH(- vert_z0hig / vert_h_hig))
       ENDDO

       allocate(p(2: llm))
       p = preff * EXP(- zz(2:) / vert_scale)
       ya = pah / preff
       s(2: llm) = hybrid(p)

       call compute_ab
    case("read_hybrid")
       ! Read "ap" and "bp". First line is skipped (title line). "ap"
       ! should be in Pa. First couple of values should correspond to
       ! the surface, that is : "bp" should be in descending order.
       call new_unit(unit)
       open(unit, file="hybrid.csv", status="old", action="read", &
            position="rewind")
       read(unit, fmt=*) ! skip title line
       do l = 1, llm + 1
          read(unit, fmt=*) ap(l), bp(l)
       end do
       close(unit)
       ! Quick check:
       call assert(ap(1) == 0., ap(llm + 1) == 0., bp(1) == 1., &
            bp(llm + 1) == 0., "disvert: bad ap or bp values")
       s(2: llm) = ap(2: llm) / pah + bp(2: llm)
    case("read_pressure")
       ! Read pressure values, in Pa, in descending order, from preff
       ! to 0. First line is skipped (title line).
       call read_column(p, "pressure.txt", first = 2)
       call assert(size(p) == llm + 1, "disvert: bad number of pressure values")
       ! Quick check:
       call assert(p(1) == preff, p(llm + 1) == 0., &
            "disvert: bad pressure values")
       ya = pah / preff
       s(2: llm) = hybrid(p(2: llm))
       call compute_ab
    case default
       print *, 'Wrong value for "vert_sampling"'
       stop 1
    END select

    forall (l = 1: llm) presnivs(l) = 0.5 &
         * (ap(l) + bp(l) * preff + ap(l+1) + bp(l+1) * preff)

  END SUBROUTINE disvert

  !**********************************************************

  subroutine compute_ab

    ! Calcul de "ap" et "bp".

    where (s >= 1. / sqrt(1. - log(tiny(0.))))
       bp = exp(1. - 1. / s**2)
    elsewhere
       bp = 0.
    end where

    ap = pah * (s - bp)

  end subroutine compute_ab

  !**********************************************************

  function hybrid(p)

    ! This procedure computes the hybrid sigma-pressure coordinate
    ! from pressure values, assuming some reference surface
    ! pressure. The procedure assumes, and does not check, that
    ! pressure values are given in descending order.

    use numer_rec_95, only: rtsafe

    real, intent(in):: p(:) ! pressure (in hPa)
    real hybrid(size(p)) ! hybrid sigma-pressure coordinate

    ! Local:
    integer l

    !-------------------------------------------------------

    y = p(1) / preff
    hybrid(1) = rtsafe(funcd, x1 = 0., x2 = 1., xacc = 1e-4)

    do l = 2, size(p)
       y = p(l) / preff
       ! Assuming descending order in pressure:
       hybrid(l) = rtsafe(funcd, x1 = 0., x2 = hybrid(l - 1), &
            xacc = hybrid(l - 1) * 1e-4)
    end do

  end function hybrid

  !**********************************************************

  SUBROUTINE funcd(s, fval, fderiv)

    REAL, INTENT(IN):: s
    REAL, INTENT(OUT):: fval, fderiv

    ! Local:
    real b

    !------------------------------------

    if (s**3 > 1. / huge(0.)) then
       b = exp(1. - 1. / s**2)
       fval = ya * s + (1. - ya) * b - y
       fderiv = ya + 2 * (1. - ya) * b / s**3
    else
       fval = ya * s - y
       fderiv = ya
    end if

  END SUBROUTINE funcd

end module disvert_m
