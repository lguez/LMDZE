module inidissip_m

  use dimens_m, only: llm

  IMPLICIT NONE

  private llm

  REAL dtdiss
  integer idissip ! periode de la dissipation (en pas)
  real tetaudiv(llm),tetaurot(llm),tetah(llm)   
  real cdivu, crot, cdivh

contains

  SUBROUTINE inidissip(lstardis, nitergdiv, nitergrot, niterh, tetagdiv, &
       tetagrot, tetatemp)

    ! From dyn3d/inidissip.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Initialisation de la dissipation horizontale                        

    USE comconst, ONLY : dtvr
    USE comvert, ONLY : preff, presnivs
    USE conf_gcm_m, ONLY : iperiod
    USE dimens_m, ONLY : jjm, llm
    USE paramet_m, ONLY : iip1, ip1jm, ip1jmp1, jjp1
    use new_unit_m, only: new_unit

    LOGICAL, intent(in):: lstardis
    INTEGER, intent(in):: nitergdiv, nitergrot, niterh
    REAL, intent(in):: tetagdiv, tetagrot, tetatemp

    ! Variables local to the procedure:
    REAL zvert(llm)
    REAL zh(ip1jmp1), zu(ip1jmp1), zv(ip1jm), deltap(ip1jmp1, llm)
    REAL ullm, vllm, umin, vmin, zhmin, zhmax
    REAL zllm, z1llm
    INTEGER l, ij, idum, ii, unit
    REAL tetamin
    REAL ran1

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: inidissip'

    !   calcul des valeurs propres des operateurs par methode iterrative:   

    crot = -1.
    cdivu = -1.
    cdivh = -1.

    !   calcul de la valeur propre de divgrad:                              

    idum = 0
    DO l = 1, llm
       DO ij = 1, ip1jmp1
          deltap(ij, l) = 1.
       END DO
    END DO

    idum = -1
    zh(1) = ran1(idum) - .5
    idum = 0
    DO ij = 2, ip1jmp1
       zh(ij) = ran1(idum) - .5
    END DO

    CALL filtreg(zh, jjp1, 1, 2, 1, .TRUE., 1)

    CALL minmax(iip1*jjp1, zh, zhmin, zhmax)

    IF (zhmin>=zhmax) THEN
       PRINT *, '  Inidissip  zh min max  ', zhmin, zhmax
       STOP 'probleme generateur alleatoire dans inidissip'
    END IF

    zllm = abs(zhmax)
    DO l = 1, 50
       IF (lstardis) THEN
          CALL divgrad2(1, zh, deltap, niterh, zh)
       ELSE
          CALL divgrad(1, zh, niterh, zh)
       END IF

       CALL minmax(iip1*jjp1, zh, zhmin, zhmax)

       zllm = abs(zhmax)
       z1llm = 1./zllm
       DO ij = 1, ip1jmp1
          zh(ij) = zh(ij)*z1llm
       END DO
    END DO

    IF (lstardis) THEN
       cdivh = 1./zllm
    ELSE
       cdivh = zllm**(-1./niterh)
    END IF

    !   calcul des valeurs propres de gradiv (ii =1) et  nxgrarot(ii=2)     

    PRINT *, 'calcul des valeurs propres'

    DO  ii = 1, 2

       DO ij = 1, ip1jmp1
          zu(ij) = ran1(idum) - .5
       END DO
       CALL filtreg(zu, jjp1, 1, 2, 1, .TRUE., 1)
       DO ij = 1, ip1jm
          zv(ij) = ran1(idum) - .5
       END DO
       CALL filtreg(zv, jjm, 1, 2, 1, .FALSE., 1)

       CALL minmax(iip1*jjp1, zu, umin, ullm)
       CALL minmax(iip1*jjm, zv, vmin, vllm)

       ullm = abs(ullm)
       vllm = abs(vllm)

       DO  l = 1, 50
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

          CALL minmax(iip1*jjp1, zu, umin, ullm)
          CALL minmax(iip1*jjm, zv, vmin, vllm)

          ullm = abs(ullm)
          vllm = abs(vllm)

          zllm = max(ullm, vllm)
          z1llm = 1./zllm
          DO ij = 1, ip1jmp1
             zu(ij) = zu(ij)*z1llm
          END DO
          DO ij = 1, ip1jm
             zv(ij) = zv(ij)*z1llm
          END DO
       end DO

       IF (ii==1) THEN
          IF (lstardis) THEN
             cdivu = 1./zllm
          ELSE
             cdivu = zllm**(-1./nitergdiv)
          END IF
       ELSE
          IF (lstardis) THEN
             crot = 1./zllm
          ELSE
             crot = zllm**(-1./nitergrot)
          END IF
       END IF

    END DO

    PRINT *, 'cdivu = ', cdivu
    PRINT *, 'crot = ', crot
    PRINT *, 'cdivh = ', cdivh

    ! Variation verticale du coefficient de dissipation :
    zvert = 2. - 1. / (1. + (1. - preff / presnivs)**2)

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

    tetamin = min(1E6, minval(1. / tetaudiv), minval(1. / tetaurot), &
         minval(1. / tetah))
    PRINT *, 'tetamin = ', tetamin
    idissip = max(iperiod, int(tetamin / (2 * dtvr * iperiod)) * iperiod)
    PRINT *, 'idissip = ', idissip
    dtdiss = idissip * dtvr
    PRINT *, 'dtdiss = ', dtdiss

  END SUBROUTINE inidissip

end module inidissip_m
