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

    ! Arguments d'entree:
    REAL, intent(in):: paprs(klon, klev+1) ! pression (Pa) a inter-couche
    REAL, intent(in):: pplay(klon, klev) ! pression (Pa) au milieu de couche
    REAL, intent(in):: rain(klon) ! pluie convective (kg/m2/s)
    REAL, intent(in):: snow(klon) ! neige convective (kg/m2/s)
    INTEGER, intent(in):: kbot(klon) ! bas de la convection
    INTEGER, intent(in):: ktop(klon) ! sommet de la convection

    ! Arguments de sortie:
    REAL, intent(out):: diafra(klon, klev) ! fraction nuageuse diagnostiquee
    REAL, intent(out):: dialiq(klon, klev) ! eau liquide nuageuse

    ! Constantes ajustables:
    REAL CANVA, CANVB, CANVH
    PARAMETER (CANVA=2.0, CANVB=0.3, CANVH=0.4)
    REAL CCA, CCB, CCC
    PARAMETER (CCA=0.125, CCB=1.5, CCC=0.8)
    REAL CCFCT, CCSCAL
    PARAMETER (CCFCT=0.400)
    PARAMETER (CCSCAL=1.0E+11)
    REAL CETAHB
    PARAMETER (CETAHB=0.45)
    REAL CCLWMR
    PARAMETER (CCLWMR=1.E-04)
    REAL ZEPSCR
    PARAMETER (ZEPSCR=1.0E-10)

    ! Variables locales:
    INTEGER i, k
    REAL zcc(klon)

    ! Initialisation:

    DO k = 1, klev
       DO i = 1, klon
          diafra(i, k) = 0.0
          dialiq(i, k) = 0.0
       ENDDO
    ENDDO

    DO i = 1, klon ! Calculer la fraction nuageuse
       zcc(i) = 0.0
       IF((rain(i)+snow(i)) > 0.) THEN
          zcc(i)= CCA * LOG(MAX(ZEPSCR, (rain(i)+snow(i))*CCSCAL))-CCB
          zcc(i)= MIN(CCC, MAX(0.0, zcc(i)))
       ENDIF
    ENDDO

    DO i = 1, klon ! pour traiter les enclumes
       diafra(i, ktop(i)) = MAX(diafra(i, ktop(i)), zcc(i)*CCFCT)
       IF ((zcc(i) >= CANVH) .AND. &
            (pplay(i, ktop(i)) <= CETAHB*paprs(i, 1))) &
            diafra(i, ktop(i)) = MAX(diafra(i, ktop(i)), &
            MAX(zcc(i)*CCFCT, CANVA*(zcc(i)-CANVB)))
       dialiq(i, ktop(i))=CCLWMR*diafra(i, ktop(i))
    ENDDO

    DO k = 1, klev ! nuages convectifs (sauf enclumes)
       DO i = 1, klon
          IF (k < ktop(i) .AND. k >= kbot(i)) THEN
             diafra(i, k)=MAX(diafra(i, k), zcc(i)*CCFCT)
             dialiq(i, k)=CCLWMR*diafra(i, k)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE diagcld1

end module diagcld1_m
