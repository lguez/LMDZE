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
    ! Initialisation de la dissipation horizontale 

    USE comconst, ONLY : dtvr
    use comdissnew, only: lstardis, nitergdiv, nitergrot, niterh, tetagdiv, &
         tetagrot, tetatemp
    USE comvert, ONLY : preff, presnivs
    USE conf_gcm_m, ONLY : iperiod
    USE dimens_m, ONLY : jjm, llm
    USE paramet_m, ONLY : iip1, ip1jm, ip1jmp1, jjp1
    use jumble, only: new_unit
    use filtreg_m, only: filtreg

    ! Variables local to the procedure:
    REAL zvert(llm), max_zvert
    REAL zh(ip1jmp1), zu(ip1jmp1), zv(ip1jm), deltap(ip1jmp1, llm)
    REAL zhmin, zhmax
    REAL zllm
    INTEGER l, ij, idum, ii, unit
    REAL tetamin ! in s
    REAL ran1

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: inidissip'

    ! Calcul des valeurs propres des opérateurs par méthode itérative :

    crot = -1.
    cdivu = -1.
    cdivh = -1.

    ! Calcul de la valeur propre de divgrad :

    deltap = 1.
    idum = -1
    zh(1) = ran1(idum) - 0.5
    idum = 0
    DO ij = 2, ip1jmp1
       zh(ij) = ran1(idum) - 0.5
    END DO

    CALL filtreg(zh, jjp1, 1, 2, 1, .TRUE., 1)

    CALL minmax(iip1*jjp1, zh, zhmin, zhmax)
    IF (zhmin >= zhmax) THEN
       PRINT *, 'zhmin zhmax', zhmin, zhmax
       print *, 'Problème générateur aléatoire dans inidissip'
       STOP 1
    END IF

    DO l = 1, 50
       IF (lstardis) THEN
          CALL divgrad2(1, zh, deltap, niterh, zh)
       ELSE
          CALL divgrad(1, zh, niterh, zh)
       END IF

       zllm = abs(maxval(zh))
       zh = zh / zllm
    END DO

    IF (lstardis) THEN
       cdivh = 1. / zllm
    ELSE
       cdivh = zllm**(- 1. / niterh)
    END IF

    ! Calcul des valeurs propres de gradiv (ii = 1) et nxgrarot (ii = 2) 

    PRINT *, 'Calcul des valeurs propres'

    DO ii = 1, 2
       DO ij = 1, ip1jmp1
          zu(ij) = ran1(idum) - 0.5
       END DO
       CALL filtreg(zu, jjp1, 1, 2, 1, .TRUE., 1)
       DO ij = 1, ip1jm
          zv(ij) = ran1(idum) - 0.5
       END DO
       CALL filtreg(zv, jjm, 1, 2, 1, .FALSE., 1)

       DO l = 1, 50
          IF (ii==1) THEN
             IF (lstardis) THEN
                CALL gradiv2(1, zu, zv, nitergdiv, zu, zv)
             ELSE
                CALL gradiv(1, zu, zv, nitergdiv, zu, zv)
             END IF
          ELSE
             IF (lstardis) THEN
                CALL nxgraro2(1, zu, zv, nitergrot, zu, zv)
             ELSE
                CALL nxgrarot(1, zu, zv, nitergrot, zu, zv)
             END IF
          END IF

          zllm = max(abs(maxval(zu)), abs(maxval(zv)))
          zu = zu / zllm
          zv = zv / zllm
       end DO

       IF (ii==1) THEN
          IF (lstardis) THEN
             cdivu = 1. / zllm
          ELSE
             cdivu = zllm**(- 1. / nitergdiv)
          END IF
       ELSE
          IF (lstardis) THEN
             crot = 1./zllm
          ELSE
             crot = zllm**(-1. / nitergrot)
          END IF
       END IF
    END DO

    PRINT *, 'cdivu = ', cdivu
    PRINT *, 'crot = ', crot
    PRINT *, 'cdivh = ', cdivh

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
    tetamin = min(1E6, tetagdiv / max_zvert, tetagrot / max_zvert, &
         tetatemp / max_zvert)
    PRINT *, 'tetamin = ', tetamin
    idissip = max(1, int(tetamin / (2 * dtvr * iperiod))) * iperiod
    PRINT *, 'idissip = ', idissip
    dtdiss = idissip * dtvr
    PRINT *, 'dtdiss = ', dtdiss

  END SUBROUTINE inidissip

end module inidissip_m
