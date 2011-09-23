SUBROUTINE ajsec(paprs, pplay, t,q, d_t,d_q)

  ! From LMDZ4/libf/phylmd/ajsec.F,v 1.1.1.1 2004/05/19 12:53:08

  use dimens_m
  use dimphy
  use SUPHEC_M
  IMPLICIT none
  !======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: ajustement sec (adaptation du GCM du LMD)
  !======================================================================
  ! Arguments:
  ! t-------input-R- Temperature
  !
  ! d_t-----output-R-Incrementation de la temperature
  !======================================================================
  REAL, intent(in):: paprs(klon,klev+1)
  real, intent(in):: pplay(klon,klev)
  REAL, intent(in):: t(klon,klev)
  real q(klon,klev)
  REAL d_t(klon,klev), d_q(klon,klev)
  !
  INTEGER limbas, limhau ! les couches a ajuster
  !cc      PARAMETER (limbas=klev-3, limhau=klev)
  PARAMETER (limbas=1, limhau=klev)
  !
  LOGICAL mixq
  !cc      PARAMETER (mixq=.TRUE.)
  PARAMETER (mixq=.FALSE.)
  !
  REAL zh(klon,klev)
  REAL zq(klon,klev)
  REAL zpk(klon,klev)
  REAL zpkdp(klon,klev)
  REAL hm, sm, qm
  LOGICAL modif(klon), down
  INTEGER i, k, k1, k2
  !
  ! Initialisation:
  !
  DO k = 1, klev
     DO i = 1, klon
        d_t(i,k) = 0.0
        d_q(i,k) = 0.0
     ENDDO
  ENDDO
  !------------------------------------- detection des profils a modifier
  DO k = limbas, limhau
     DO i = 1, klon
        zpk(i,k) = pplay(i,k)**RKAPPA
        zh(i,k) = RCPD * t(i,k)/ zpk(i,k)
        zq(i,k) = q(i,k)
     ENDDO
  ENDDO
  !
  DO k = limbas, limhau
     DO i = 1, klon
        zpkdp(i,k) = zpk(i,k) * (paprs(i,k)-paprs(i,k+1))
     ENDDO
  ENDDO
  !
  DO i = 1, klon
     modif(i) = .FALSE.
  ENDDO
  DO k = limbas+1, limhau
     DO i = 1, klon
        IF (.NOT.modif(i)) THEN
           IF ( zh(i,k).LT.zh(i,k-1) ) modif(i) = .TRUE.
        ENDIF
     ENDDO
  ENDDO
  !------------------------------------- correction des profils instables
  DO  i = 1, klon
     IF (modif(i)) THEN
        k2 = limbas
        do
           k2 = k2 + 1
           IF (k2 .GT. limhau) exit
           IF (zh(i,k2) .LT. zh(i,k2-1)) THEN
              k1 = k2 - 1
              k = k1
              sm = zpkdp(i,k2)
              hm = zh(i,k2)
              qm = zq(i,k2)
              do
                 sm = sm +zpkdp(i,k)
                 hm = hm +zpkdp(i,k) * (zh(i,k)-hm) / sm
                 qm = qm +zpkdp(i,k) * (zq(i,k)-qm) / sm
                 down = .FALSE.
                 IF (k1 .ne. limbas) THEN
                    IF (hm .LT. zh(i,k1-1)) down = .TRUE.
                 ENDIF
                 IF (down) THEN
                    k1 = k1 - 1
                    k = k1
                 ELSE
                    IF ((k2 .EQ. limhau)) exit
                    IF ((zh(i,k2+1).GE.hm)) exit
                    k2 = k2 + 1
                    k = k2
                 ENDIF
              end do

              !------------ nouveau profil : constant (valeur moyenne)
              DO k = k1, k2
                 zh(i,k) = hm
                 zq(i,k) = qm
              ENDDO
              k2 = k2 + 1
           ENDIF
        end do
     ENDIF
  end DO
  !
  DO k = limbas, limhau
     DO i = 1, klon
        d_t(i,k) = zh(i,k)*zpk(i,k)/RCPD - t(i,k)
        d_q(i,k) = zq(i,k) - q(i,k)
     ENDDO
  ENDDO
  !
  IF (limbas.GT.1) THEN
     DO k = 1, limbas-1
        DO i = 1, klon
           d_t(i,k) = 0.0
           d_q(i,k) = 0.0
        ENDDO
     ENDDO
  ENDIF
  !
  IF (limhau.LT.klev) THEN
     DO k = limhau+1, klev
        DO i = 1, klon
           d_t(i,k) = 0.0
           d_q(i,k) = 0.0
        ENDDO
     ENDDO
  ENDIF
  !
  IF (.NOT.mixq) THEN
     DO k = 1, klev
        DO i = 1, klon
           d_q(i,k) = 0.0
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE ajsec
