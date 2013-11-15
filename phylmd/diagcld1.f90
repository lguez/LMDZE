      SUBROUTINE diagcld1(paprs,pplay,rain,snow,kbot,ktop, &
                         diafra,dialiq)
      use dimens_m
      use dimphy
      use SUPHEC_M
      IMPLICIT none
!
! Laurent Li (LMD/CNRS), le 12 octobre 1998
!                        (adaptation du code ECMWF)
!
! Dans certains cas, le schema pronostique des nuages n'est
! pas suffisament performant. On a donc besoin de diagnostiquer
! ces nuages. Je dois avouer que c'est une frustration.
!
!
! Arguments d'entree:
      REAL, intent(in):: paprs(klon,klev+1) ! pression (Pa) a inter-couche
      REAL, intent(in):: pplay(klon,klev) ! pression (Pa) au milieu de couche
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! humidite specifique (Kg/Kg)
      REAL rain(klon) ! pluie convective (kg/m2/s)
      REAL snow(klon) ! neige convective (kg/m2/s)
      INTEGER ktop(klon) ! sommet de la convection
      INTEGER kbot(klon) ! bas de la convection
!
! Arguments de sortie:
      REAL diafra(klon,klev) ! fraction nuageuse diagnostiquee
      REAL dialiq(klon,klev) ! eau liquide nuageuse
!
! Constantes ajustables:
      REAL CANVA, CANVB, CANVH
      PARAMETER (CANVA=2.0, CANVB=0.3, CANVH=0.4)
      REAL CCA, CCB, CCC
      PARAMETER (CCA=0.125, CCB=1.5, CCC=0.8)
      REAL CCFCT, CCSCAL
      PARAMETER (CCFCT=0.400)
      PARAMETER (CCSCAL=1.0E+11)
      REAL CETAHB, CETAMB
      PARAMETER (CETAHB=0.45, CETAMB=0.80)
      REAL CCLWMR
      PARAMETER (CCLWMR=1.E-04)
      REAL ZEPSCR
      PARAMETER (ZEPSCR=1.0E-10)
!
! Variables locales:
      INTEGER i, k
      REAL zcc(klon)
!
! Initialisation:
!
      DO k = 1, klev
      DO i = 1, klon
         diafra(i,k) = 0.0
         dialiq(i,k) = 0.0
      ENDDO
      ENDDO
!
      DO i = 1, klon ! Calculer la fraction nuageuse
      zcc(i) = 0.0
      IF((rain(i)+snow(i)).GT.0.) THEN
         zcc(i)= CCA * LOG(MAX(ZEPSCR,(rain(i)+snow(i))*CCSCAL))-CCB
         zcc(i)= MIN(CCC,MAX(0.0,zcc(i)))
      ENDIF
      ENDDO
!
      DO i = 1, klon ! pour traiter les enclumes
      diafra(i,ktop(i)) = MAX(diafra(i,ktop(i)),zcc(i)*CCFCT)
      IF ((zcc(i).GE.CANVH) .AND. &
          (pplay(i,ktop(i)).LE.CETAHB*paprs(i,1))) &
       diafra(i,ktop(i)) = MAX(diafra(i,ktop(i)), &
                               MAX(zcc(i)*CCFCT,CANVA*(zcc(i)-CANVB)))
      dialiq(i,ktop(i))=CCLWMR*diafra(i,ktop(i))
      ENDDO
!
      DO k = 1, klev ! nuages convectifs (sauf enclumes)
      DO i = 1, klon
      IF (k.LT.ktop(i) .AND. k.GE.kbot(i)) THEN
         diafra(i,k)=MAX(diafra(i,k),zcc(i)*CCFCT)
         dialiq(i,k)=CCLWMR*diafra(i,k)
      ENDIF
      ENDDO
      ENDDO
!
      RETURN
      END
