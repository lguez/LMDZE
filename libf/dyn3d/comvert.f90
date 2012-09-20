module comvert

  use dimens_m, only: llm

  implicit none

  private llm

  real ap(llm+1), pa ! in Pa
  real bp(llm+1)
  real presnivs(llm) ! pressions approximatives des milieux de couches, en Pa
  real, parameter:: preff = 101325. ! in Pa
  real nivsigs(llm), nivsig(llm+1)

  save

contains

  SUBROUTINE disvert

    ! From dyn3d/disvert.F, v 1.1.1.1 2004/05/19 12:53:05
    ! Auteur : P. Le Van

    ! This procedure sets the vertical grid.
    ! It defines the host variables "ap", "bp", "presnivs", "nivsigs"
    ! and "nivsig".
    ! "pa" should be defined before this procedure is called.

    use nr_util, only: pi
    use jumble, only: new_unit
    use unit_nml_m, only: unit_nml

    REAL s(llm+1)
    ! "s(l)" is the atmospheric hybrid sigma-pressure coordinate at
    ! the interface between layers "l" and "l-1"

    real ds(llm)
    ! "ds(l)" : épaisseur de la couche "l" dans la coordonnée "s"

    INTEGER l, unit
    REAL alpha, x(llm)

    character(len=7):: s_sampling = "tropo"
    ! (other allowed values are "param", "strato1", "strato" and "read")

    real:: h = 7. ! scale height, in km
    ! (used only if "s_sampling" == "param" or "strato1")

    ! These variables are used only in the case 's_sampling == "param"':
    real:: deltaz = 0.04 ! épaisseur de la première couche
    real:: beta = 1.3 ! facteur d'accroissement en haut
    real:: k0 = 20. ! nombre de couches dans la transition surface
    real:: k1 = 1.2 ! nombre de couches dans la transition haute

    REAL ZZ(llm + 1), DZ(llm) ! in km

    namelist /disvert_nml/h, deltaz, beta, k0, k1, s_sampling

    !-----------------------------------------------------------------------

    print *, "Call sequence information: disvert"

    forall (l = 1: llm) nivsigs(l) = REAL(l)
    forall (l = 1: llm + 1) nivsig(l) = REAL(l)

    ! Compute "s":

    print *, "Enter namelist 'disvert_nml'."
    read(unit=*, nml=disvert_nml)
    write(unit_nml, nml=disvert_nml)

    select case (s_sampling)
    case ("param")
       s(1) = 1.
       s(llm+1) = 0.
       alpha = deltaz / tanh(1./k0) * 2.
       forall (l = 2: llm) s(l) &
            = cosh((l - 1) / k0) **(- alpha * k0 / h) &
            * exp(- alpha / h * tanh((llm - k1) / k0) &
            * beta **(l - 1 - (llm - k1)) / log(beta))
    case ("tropo")
       s(1) = 1.
       s(llm+1) = 0.
       forall (l = 1: llm) ds(l) &
            = 1. + 7. * SIN(pi * (REAL(l)-0.5) / real(llm+1))**2
       ds = ds / sum(ds)

       DO l = llm, 2, -1
          s(l) = s(l+1) + ds(l)
       ENDDO
    case ("strato1")
       ! F. Lott 70 niveaux et plus
       s(1) = 1.
       s(llm+1) = 0.
       forall (l = 1: llm) dz(l) = 1.56 + TANH(REAL(l - 12) / 5.) &
               + TANH(REAL(l - llm) / 10.) / 2.

       zz(1) = 0.
       DO l = 2, llm + 1
          zz(l) = zz(l - 1) + dz(l - 1)
       end DO

       s(2:llm) = (exp(- zz(2:llm) / h) - exp(- zz(llm + 1) / h)) &
            / (1. - exp(- zz(llm + 1) / h))
    case ("strato")
       ! Recommended by F. Lott for a domain including the stratosphere
       s(1) = 1.
       s(llm+1) = 0.
       forall (l = 1: llm) x(l) = pi * (l - 0.5) / (llm + 1)

       ds = (1. + 7. * SIN(x)**2) * (1. - tanh(2 * x / pi - 1.))**2 / 4.
       ds = ds / sum(ds)

       DO l = llm, 2, -1
          s(l) = s(l+1) + ds(l)
       ENDDO
    case("read")
       call new_unit(unit)
       open(unit, file="hybrid.csv", status="old", action="read", &
            position="rewind")
       read(unit, fmt=*) ! skip title line
       do l = 1, llm + 1
          read(unit, fmt=*) s(l)
       end do
       close(unit)
       ! Quick check:
       if (s(1) /= 1. .or. s(llm + 1) /= 0. .or. s(1) <= s(2)) then
          print *, '"s" should be in descending order, from 1 to 0.'
          stop 1
       end if
    case default
       print *, 'Wrong value for "s_sampling"'
       stop 1
    END select

    ! Calcul de "ap" et "bp" :
    bp(:llm) = EXP(1. - 1. / s(:llm)**2)
    bp(llm + 1) = 0.
    ap = pa * (s - bp)

    forall (l = 1: llm) presnivs(l) = 0.5 &
         * (ap(l) + bp(l) * preff + ap(l+1) + bp(l+1) * preff)

  END SUBROUTINE disvert

end module comvert
