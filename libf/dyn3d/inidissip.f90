module inidissip_m

  use dimens_m, only: llm

  IMPLICIT NONE

  private llm

  REAL dtdiss
  integer idissip ! période de la dissipation (en pas de temps)
  real tetaudiv(llm), tetaurot(llm), tetah(llm) 
  real cdivu, crot, cdivh

contains

  SUBROUTINE inidissip

    ! From dyn3d/inidissip.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Initialisation de la dissipation horizontale. Calcul des valeurs
    ! propres des opérateurs par méthode itérative.

    USE comconst, ONLY : dtvr
    use comdissnew, only: lstardis, nitergdiv, nitergrot, niterh, tetagdiv, &
         tetagrot, tetatemp
    USE comvert, ONLY : preff, presnivs
    USE conf_gcm_m, ONLY : iperiod
    USE dimens_m, ONLY : iim, jjm, llm
    USE paramet_m, ONLY : jjp1
    use jumble, only: new_unit
    use filtreg_m, only: filtreg
    use gradiv2_m, only: gradiv2

    ! Variables local to the procedure:
    REAL zvert(llm), max_zvert
    REAL, dimension(iim + 1, jjm + 1):: zh, zu
    real zv(iim + 1, jjm), deltap(iim + 1, jjm + 1, llm)
    REAL zllm
    INTEGER l, seed_size, ii, unit
    REAL tetamin ! in s

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: inidissip'
    call random_seed(size=seed_size)
    call random_seed(put=(/(0, ii = 1, seed_size)/))

    PRINT *, 'Calcul des valeurs propres de divgrad'
    deltap = 1.
    call random_number(zh)
    zh = zh - 0.5
    CALL filtreg(zh, jjp1, 1, 2, 1, .TRUE., 1)

    DO l = 1, 50
       IF (lstardis) THEN
          CALL divgrad2(1, zh, deltap, niterh, zh, -1.)
       ELSE
          CALL divgrad(1, zh, niterh, zh, -1.)
       END IF

       zllm = abs(maxval(zh))
       zh = zh / zllm
    END DO

    IF (lstardis) THEN
       cdivh = 1. / zllm
    ELSE
       cdivh = zllm**(- 1. / niterh)
    END IF
    PRINT *, 'cdivh = ', cdivh

    PRINT *, 'Calcul des valeurs propres de gradiv'
    call random_number(zu)
    zu = zu - 0.5
    CALL filtreg(zu, jjp1, 1, 2, 1, .TRUE., 1)
    call random_number(zv)
    zv = zv - 0.5
    CALL filtreg(zv, jjm, 1, 2, 1, .FALSE., 1)

    DO l = 1, 50
       IF (lstardis) THEN
          CALL gradiv2(1, zu, zv, nitergdiv, zu, zv, -1.)
       ELSE
          CALL gradiv(1, zu, zv, nitergdiv, zu, zv, -1.)
       END IF

       zllm = max(abs(maxval(zu)), abs(maxval(zv)))
       zu = zu / zllm
       zv = zv / zllm
    end DO

    IF (lstardis) THEN
       cdivu = 1. / zllm
    ELSE
       cdivu = zllm**(- 1. / nitergdiv)
    END IF
    PRINT *, 'cdivu = ', cdivu

    PRINT *, 'Calcul des valeurs propres de nxgrarot'
    call random_number(zu)
    zu = zu - 0.5
    CALL filtreg(zu, jjp1, 1, 2, 1, .TRUE., 1)
    call random_number(zv)
    zv = zv - 0.5
    CALL filtreg(zv, jjm, 1, 2, 1, .FALSE., 1)

    DO l = 1, 50
       IF (lstardis) THEN
          CALL nxgraro2(1, zu, zv, nitergrot, zu, zv, -1.)
       ELSE
          CALL nxgrarot(1, zu, zv, nitergrot, zu, zv, -1.)
       END IF

       zllm = max(abs(maxval(zu)), abs(maxval(zv)))
       zu = zu / zllm
       zv = zv / zllm
    end DO

    IF (lstardis) THEN
       crot = 1. / zllm
    ELSE
       crot = zllm**(-1. / nitergrot)
    END IF
    PRINT *, 'crot = ', crot

    ! Variation verticale du coefficient de dissipation :
    zvert = 2. - 1. / (1. + (preff / presnivs - 1.)**2)
    ! (between 1 and 2)

    tetaudiv = zvert / tetagdiv
    tetaurot = zvert / tetagrot
    tetah = zvert / tetatemp

    call new_unit(unit)
    open(unit, file="inidissip.csv", status="replace", action="write")
    write(unit, fmt=*) "tetaudiv tetaurot tetah" ! title line
    do l = 1, llm
       write(unit, fmt=*) tetaudiv(l), tetaurot(l), tetah(l)
    end do
    close(unit)
    print *, 'Created file "inidissip.csv".'

    max_zvert = maxval(zvert)
    tetamin = min(1e6, tetagdiv / max_zvert, tetagrot / max_zvert, &
         tetatemp / max_zvert)
    PRINT *, 'tetamin = ', tetamin
    idissip = max(1, int(tetamin / (2 * dtvr * iperiod))) * iperiod
    PRINT *, 'idissip = ', idissip
    dtdiss = idissip * dtvr
    PRINT *, 'dtdiss = ', dtdiss

  END SUBROUTINE inidissip

end module inidissip_m
