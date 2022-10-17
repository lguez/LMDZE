module diagcld1_m

  IMPLICIT none

contains

  SUBROUTINE diagcld1(paprs, pplay, rain, snow, kbot, ktop, diafra, dialiq)

    ! Laurent Li (LMD/CNRS), le 12 octobre 1998
    ! (adaptation du code ECMWF)

    ! Dans certains cas, le sch\'ema pronostique des nuages n'est pas
    ! suffisamment performant. On a donc besoin de diagnostiquer ces
    ! nuages. Je dois avouer que c'est une frustration.

    use dimphy, only: klon, klev

    REAL, intent(in):: paprs(klon, klev+1) ! pression (Pa) a inter-couche
    REAL, intent(in):: pplay(klon, klev) ! pression (Pa) au milieu de couche
    REAL, intent(in):: rain(klon) ! pluie convective (kg/m2/s)
    REAL, intent(in):: snow(klon) ! neige convective (kg/m2/s)
    INTEGER, intent(in):: kbot(klon) ! bas de la convection
    INTEGER, intent(in):: ktop(klon) ! sommet de la convection
    REAL, intent(out):: diafra(klon, klev) ! fraction nuageuse diagnostiquee
    REAL, intent(out):: dialiq(klon, klev) ! eau liquide nuageuse

    ! Local:

    ! Constantes ajustables:
    REAL, PARAMETER:: CANVA = 2., CANVB = 0.3, CANVH = 0.4
    REAL, PARAMETER:: CCA = 0.125, CCB = 1.5, CCC = 0.8
    REAL, PARAMETER:: CCFCT = 0.400
    REAL, PARAMETER:: CCSCAL = 1E+11
    REAL, PARAMETER:: CETAHB = 0.45
    REAL, PARAMETER:: CCLWMR = 1E-4
    REAL, PARAMETER:: ZEPSCR = 1E-10

    INTEGER i, k
    REAL zcc(klon)

    !---------------------------------------------------------------------

    ! Initialisation:
    diafra = 0.
    dialiq = 0.

    ! Calculer la fraction nuageuse:
    DO i = 1, klon
       zcc(i) = 0.
       IF((rain(i)+snow(i)) > 0.) THEN
          zcc(i) = CCA * LOG(MAX(ZEPSCR, (rain(i)+snow(i))*CCSCAL))-CCB
          zcc(i) = MIN(CCC, MAX(0., zcc(i)))
       ENDIF
    ENDDO

    ! Pour traiter les enclumes:
    DO i = 1, klon
       diafra(i, ktop(i)) = MAX(diafra(i, ktop(i)), zcc(i)*CCFCT)
       IF ((zcc(i) >= CANVH) .AND. &
            (pplay(i, ktop(i)) <= CETAHB*paprs(i, 1))) &
            diafra(i, ktop(i)) = MAX(diafra(i, ktop(i)), &
            MAX(zcc(i)*CCFCT, CANVA*(zcc(i)-CANVB)))
       dialiq(i, ktop(i)) = CCLWMR*diafra(i, ktop(i))
    ENDDO

    ! Nuages convectifs (sauf enclumes):
    DO k = 1, klev
       DO i = 1, klon
          IF (k < ktop(i) .AND. k >= kbot(i)) THEN
             diafra(i, k) = MAX(diafra(i, k), zcc(i)*CCFCT)
             dialiq(i, k) = CCLWMR*diafra(i, k)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE diagcld1

end module diagcld1_m
