module disvert_m

  use dimens_m, only: llm

  implicit none

  private llm

  real, save:: ap(llm+1), pa ! in Pa
  real, save:: bp(llm+1)

  real, save:: presnivs(llm)
  ! pressions approximatives des milieux de couches, en Pa

  real, parameter:: preff = 101325. ! in Pa

contains

  SUBROUTINE disvert

    ! From dyn3d/disvert.F, v 1.1.1.1 2004/05/19 12:53:05
    ! Author: P. Le Van

    ! This procedure sets the vertical grid. It defines the host
    ! variables "ap", "bp", "presnivs". "pa" should be defined before
    ! this procedure is called.

    use jumble, only: new_unit
    use nr_util, only: pi, assert
    use unit_nml_m, only: unit_nml

    ! Local:

    REAL s(llm+1)
    ! "s(l)" is the atmospheric hybrid sigma-pressure coordinate at
    ! the interface between layers "l" and "l-1"

    real ds(llm)
    ! "ds(l)" : épaisseur de la couche "l" dans la coordonnée "s"

    INTEGER l, unit
    REAL alpha, x(llm)

    character(len=7):: vert_sampling = "tropo"
    ! other allowed values are "param", "strato" and "read"

    real:: h = 7. ! scale height, in km
    ! used only if vert_sampling == "param"

    ! These variables are used only in the case vert_sampling == "param":
    real:: deltaz = 0.04 ! épaisseur de la première couche
    real:: beta = 1.3 ! facteur d'accroissement en haut
    real:: k0 = 20. ! nombre de couches dans la transition surface
    real:: k1 = 1.2 ! nombre de couches dans la transition haute

    namelist /disvert_nml/h, deltaz, beta, k0, k1, vert_sampling

    !-----------------------------------------------------------------------

    print *, "Call sequence information: disvert"

    print *, "Enter namelist 'disvert_nml'."
    read(unit=*, nml=disvert_nml)
    write(unit_nml, nml=disvert_nml)

    select case (vert_sampling)
    case ("param")
       s(1) = 1.
       s(llm+1) = 0.
       alpha = deltaz / tanh(1./k0) * 2.
       forall (l = 2: llm) s(l) &
            = cosh((l - 1) / k0) **(- alpha * k0 / h) &
            * exp(- alpha / h * tanh((llm - k1) / k0) &
            * beta **(l - 1 - (llm - k1)) / log(beta))
       call compute_ab
    case ("tropo")
       s(1) = 1.
       s(llm+1) = 0.
       forall (l = 1: llm) ds(l) &
            = 1. + 7. * SIN(pi * (REAL(l) - 0.5) / real(llm + 1))**2
       ds = ds / sum(ds)

       DO l = llm, 2, -1
          s(l) = s(l+1) + ds(l)
       ENDDO

       call compute_ab
    case ("strato")
       ! Recommended by F. Lott for a domain including the stratosphere
       s(1) = 1.
       s(llm+1) = 0.
       forall (l = 1: llm) x(l) = pi * (l - 0.5) / (llm + 1)

       ds = (0.3 + 7. * SIN(x)**2) * (1. - tanh(2 * x / pi - 1.))**2 / 4.
       ds = ds / sum(ds)

       DO l = llm, 2, -1
          s(l) = s(l+1) + ds(l)
       ENDDO

       call compute_ab
    case("read")
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
    case default
       print *, 'Wrong value for "vert_sampling"'
       stop 1
    END select

    forall (l = 1: llm) presnivs(l) = 0.5 &
         * (ap(l) + bp(l) * preff + ap(l+1) + bp(l+1) * preff)

  contains

    subroutine compute_ab

      ! Calcul de "ap" et "bp" :
      bp(:llm) = EXP(1. - 1. / s(:llm)**2)
      bp(llm + 1) = 0.
      ap = pa * (s - bp)

    end subroutine compute_ab

  END SUBROUTINE disvert

end module disvert_m
