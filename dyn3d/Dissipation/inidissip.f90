module inidissip_m

  use dimensions, only: llm

  IMPLICIT NONE

  private llm

  REAL, protected:: dtdiss ! in s
  integer, protected:: idissip ! période de la dissipation (en pas de temps)
  real, protected:: tetaudiv(llm), tetaurot(llm), tetah(llm) ! in s-1
  real, protected:: cdivu, crot, cdivh

contains

  SUBROUTINE inidissip

    ! From dyn3d/inidissip.F, version 1.1.1.1, 2004/05/19 12:53:06

    ! Initialisation de la dissipation horizontale. Calcul des valeurs
    ! propres des opérateurs par méthode itérative.

    USE comconst, ONLY: dtvr
    use comdissnew, only: nitergdiv, nitergrot, niterh, tetagdiv, tetagrot, &
         tetatemp
    USE conf_gcm_m, ONLY: iperiod
    USE dimensions, ONLY: iim, jjm
    USE disvert_m, ONLY: preff, presnivs
    use divgrad2_m, only: divgrad2
    use filtreg_scal_m, only: filtreg_scal
    use filtreg_v_m, only: filtreg_v
    use gradiv2_m, only: gradiv2
    use jumble, only: new_unit
    use nxgraro2_m, only: nxgraro2

    ! Local:
    REAL zvert(llm), max_zvert ! no dimension
    REAL, dimension(iim + 1, jjm + 1, 1):: zh, zu, gx, divgra, deltap
    real zv(iim + 1, jjm, 1), gy(iim + 1, jjm, 1)
    REAL zllm
    INTEGER l, seed_size, ii, unit
    REAL tetamin ! in s

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: inidissip'
    call random_seed(size=seed_size)
    call random_seed(put=[(1, ii = 1, seed_size)])

    PRINT *, 'Calcul des valeurs propres de divgrad'
    deltap = 1.
    call random_number(zh)
    zh = zh - 0.5
    CALL filtreg_scal(zh, direct = .true., intensive = .true.)

    DO l = 1, 50
       CALL divgrad2(1, zh, deltap, niterh, divgra, -1.)
       zllm = abs(maxval(divgra))
       zh = divgra / zllm
    END DO

    cdivh = 1. / zllm
    PRINT *, 'cdivh = ', cdivh

    PRINT *, 'Calcul des valeurs propres de gradiv'
    call random_number(zu)
    zu = zu - 0.5
    CALL filtreg_scal(zu, direct = .true., intensive = .true.)
    call random_number(zv)
    zv = zv - 0.5
    CALL filtreg_v(zv, intensive = .true.)

    DO l = 1, 50
       CALL gradiv2(zu, zv, nitergdiv, gx, gy, -1.)
       zllm = max(abs(maxval(gx)), abs(maxval(gy)))
       zu = gx / zllm
       zv = gy / zllm
    end DO

    cdivu = 1. / zllm
    PRINT *, 'cdivu = ', cdivu

    PRINT *, 'Calcul des valeurs propres de nxgraro2'
    call random_number(zu)
    zu = zu - 0.5
    CALL filtreg_scal(zu, direct = .true., intensive = .true.)
    call random_number(zv)
    zv = zv - 0.5
    CALL filtreg_v(zv, intensive = .true.)

    DO l = 1, 50
       CALL nxgraro2(zu, zv, nitergrot, gx, gy, crot = - 1.)
       zllm = max(abs(maxval(gx)), abs(maxval(gy)))
       zu = gx / zllm
       zv = gy / zllm
    end DO

    crot = 1. / zllm
    PRINT *, 'crot = ', crot

    ! Variation verticale du coefficient de dissipation :
    zvert = 2. - 1. / (1. + (preff / presnivs - 1.)**2)
    ! (between 1 and 2)

    tetaudiv = zvert / tetagdiv
    tetaurot = zvert / tetagrot
    tetah = zvert / tetatemp

    max_zvert = maxval(zvert)
    tetamin = min(1e6, tetagdiv / max_zvert, tetagrot / max_zvert, &
         tetatemp / max_zvert)
    PRINT *, 'tetamin = ', tetamin
    idissip = max(1, int(tetamin / (2 * dtvr * iperiod))) * iperiod
    PRINT *, 'idissip = ', idissip
    dtdiss = idissip * dtvr
    PRINT *, 'dtdiss = ', dtdiss, "s"

    call new_unit(unit)
    open(unit, file="inidissip.csv", status="replace", action="write")

    ! Title line:
    write(unit, fmt=*) '"presnivs (hPa)" "dtdiss * tetaudiv" ' &
         // '"dtdiss * tetaurot" "dtdiss * tetah"'

    do l = 1, llm
       write(unit, fmt=*) presnivs(l) / 100., dtdiss * tetaudiv(l), &
            dtdiss * tetaurot(l), dtdiss * tetah(l)
    end do
    close(unit)
    print *, 'Created file "inidissip.csv".'

  END SUBROUTINE inidissip

end module inidissip_m
