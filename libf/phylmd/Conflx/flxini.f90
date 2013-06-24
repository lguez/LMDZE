module flxini_m

  IMPLICIT none

contains

  SUBROUTINE flxini(pten, pqen, pqsen, pgeo, paph, pgeoh, ptenh, pqenh, &
       pqsenh, ptu, pqu, ptd, pqd, pmfd, pmfds, pmfdq, pdmfdp, pmfu, pmfus, &
       pmfuq, pdmfup, pdpmel, plu, plude, klab,pen_u, pde_u, pen_d, pde_d)

    ! This routine interpolates large-scale fields of T,q etc. to half
    ! levels (i.e. grid for massflux scheme), and initializes values
    ! for updrafts.

    USE dimphy, ONLY: klev, klon
    use flxadjtq_m, only: flxadjtq
    USE suphec_m, ONLY: rcpd

    REAL, intent(in):: pten(klon,klev)   ! temperature (environnement)
    REAL, intent(in):: pqen(klon,klev)   ! humidite (environnement)
    REAL, intent(in):: pqsen(klon,klev)  ! humidite saturante (environnement)
    REAL, intent(in):: pgeo(klon,klev)   ! geopotentiel (g * metre)
    REAL pgeoh(klon,klev)  ! geopotentiel aux demi-niveaux
    REAL paph(klon,klev+1) ! pression aux demi-niveaux
    REAL ptenh(klon,klev)  ! temperature aux demi-niveaux
    REAL pqenh(klon,klev)  ! humidite aux demi-niveaux
    REAL pqsenh(klon,klev) ! humidite saturante aux demi-niveaux
    !
    REAL ptu(klon,klev)    ! temperature du panache ascendant (p-a)
    REAL pqu(klon,klev)    ! humidite du p-a
    REAL plu(klon,klev)    ! eau liquide du p-a
    REAL pmfu(klon,klev)   ! flux de masse du p-a
    REAL pmfus(klon,klev)  ! flux de l'energie seche dans le p-a
    REAL pmfuq(klon,klev)  ! flux de l'humidite dans le p-a
    REAL pdmfup(klon,klev) ! quantite de l'eau precipitee dans p-a
    REAL plude(klon,klev)  ! quantite de l'eau liquide jetee du
    !                              p-a a l'environnement
    REAL pdpmel(klon,klev) ! quantite de neige fondue
    !
    REAL ptd(klon,klev)    ! temperature du panache descendant (p-d)
    REAL pqd(klon,klev)    ! humidite du p-d
    REAL pmfd(klon,klev)   ! flux de masse du p-d
    REAL pmfds(klon,klev)  ! flux de l'energie seche dans le p-d
    REAL pmfdq(klon,klev)  ! flux de l'humidite dans le p-d
    REAL pdmfdp(klon,klev) ! quantite de precipitation dans p-d
    !
    REAL pen_u(klon,klev) ! quantite de masse entrainee pour p-a
    REAL pde_u(klon,klev) ! quantite de masse detrainee pour p-a
    REAL pen_d(klon,klev) ! quantite de masse entrainee pour p-d
    REAL pde_d(klon,klev) ! quantite de masse detrainee pour p-d
    !
    INTEGER  klab(klon,klev)
    LOGICAL  llflag(klon)
    INTEGER k, i, icall
    REAL zzs
    !----------------------------------------------------------------------
    ! SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
    ! ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
    !----------------------------------------------------------------------
    DO  k = 2, klev
       !
       DO i = 1, klon
          pgeoh(i,k)=pgeo(i,k)+(pgeo(i,k-1)-pgeo(i,k))*0.5
          ptenh(i,k)=(MAX(RCPD*pten(i,k-1)+pgeo(i,k-1), &
               RCPD*pten(i,k)+pgeo(i,k))-pgeoh(i,k))/RCPD
          pqsenh(i,k)=pqsen(i,k-1)
          llflag(i)=.TRUE.
       ENDDO
       !
       icall=0
       CALL flxadjtq(paph(1,k),ptenh(1,k),pqsenh(1,k),llflag,icall)
       !
       DO i = 1, klon
          pqenh(i,k)=MIN(pqen(i,k-1),pqsen(i,k-1)) &
               +(pqsenh(i,k)-pqsen(i,k-1))
          pqenh(i,k)=MAX(pqenh(i,k),0.)
       ENDDO
       !
    end DO
    !
    DO  i = 1, klon
       ptenh(i,klev)=(RCPD*pten(i,klev)+pgeo(i,klev)- &
            pgeoh(i,klev))/RCPD
       pqenh(i,klev)=pqen(i,klev)
       ptenh(i,1)=pten(i,1)
       pqenh(i,1)=pqen(i,1)
       pgeoh(i,1)=pgeo(i,1)
    end DO
    !
    DO  k = klev-1, 2, -1
       DO  i = 1, klon
          zzs = MAX(RCPD*ptenh(i,k)+pgeoh(i,k), &
               RCPD*ptenh(i,k+1)+pgeoh(i,k+1))
          ptenh(i,k) = (zzs-pgeoh(i,k))/RCPD
       end DO
    end DO
    !
    !-----------------------------------------------------------------------
    ! INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
    !-----------------------------------------------------------------------
    DO k = 1, klev
       DO i = 1, klon
          ptu(i,k) = ptenh(i,k)
          pqu(i,k) = pqenh(i,k)
          plu(i,k) = 0.
          pmfu(i,k) = 0.
          pmfus(i,k) = 0.
          pmfuq(i,k) = 0.
          pdmfup(i,k) = 0.
          pdpmel(i,k) = 0.
          plude(i,k) = 0.
          !
          klab(i,k) = 0
          !
          ptd(i,k) = ptenh(i,k)
          pqd(i,k) = pqenh(i,k)
          pmfd(i,k) = 0.0
          pmfds(i,k) = 0.0
          pmfdq(i,k) = 0.0
          pdmfdp(i,k) = 0.0
          !
          pen_u(i,k) = 0.0
          pde_u(i,k) = 0.0
          pen_d(i,k) = 0.0
          pde_d(i,k) = 0.0
       ENDDO
    ENDDO

  END SUBROUTINE flxini

end module flxini_m
