module flxini_m

  IMPLICIT none

contains

  SUBROUTINE flxini(ten, pqen, pqsen, pgeo, paph, pgeoh, ptenh, pqenh, &
       pqsenh, ptu, pqu, ptd, pqd, mfd, pmfds, pmfdq, pdmfdp, mfu, mfus, &
       mfuq, pdmfup, pdpmel, plu, plude, klab, pen_u, pde_u, pen_d, pde_d)

    ! This routine interpolates large-scale fields of T, q etc. to
    ! half levels (i. e. grid for massflux scheme), and initializes
    ! values for updrafts.

    USE dimphy, ONLY: klev, klon
    use flxadjtq_m, only: flxadjtq
    USE suphec_m, ONLY: rcpd

    REAL, intent(in):: ten(klon, klev) ! temperature (environnement)
    REAL, intent(in):: pqen(klon, klev) ! humidite (environnement)
    REAL, intent(in):: pqsen(klon, klev) ! humidite saturante (environnement)
    REAL, intent(in):: pgeo(klon, klev) ! geopotentiel (g * metre)
    REAL paph(klon, klev+1) ! pression aux demi-niveaux
    REAL pgeoh(klon, klev) ! geopotentiel aux demi-niveaux
    REAL ptenh(klon, klev) ! temperature aux demi-niveaux
    REAL pqenh(klon, klev) ! humidite aux demi-niveaux
    REAL pqsenh(klon, klev) ! humidite saturante aux demi-niveaux
    REAL ptu(klon, klev) ! temperature du panache ascendant
    REAL pqu(klon, klev) ! humidite du panache ascendant
    REAL ptd(klon, klev) ! temperature du panache descendant
    REAL pqd(klon, klev) ! humidite du panache descendant
    REAL, intent(out):: mfd(klon, klev) ! flux de masse du panache descendant
    REAL pmfds(klon, klev) ! flux de l'energie seche dans le panache descendant
    REAL pmfdq(klon, klev) ! flux de l'humidite dans le panache descendant
    REAL pdmfdp(klon, klev) ! quantite de precipitation dans panache descendant
    REAL, intent(out):: mfu(klon, klev) ! flux de masse du panache ascendant
    REAL mfus(klon, klev) ! flux de l'energie seche dans le panache ascendant
    REAL mfuq(klon, klev) ! flux de l'humidite dans le panache ascendant

    REAL pdmfup(klon, klev)
    ! quantite de l'eau precipitee dans panache ascendant

    REAL pdpmel(klon, klev) ! quantite de neige fondue
    REAL plu(klon, klev) ! eau liquide du panache ascendant

    REAL plude(klon, klev) 
    ! quantite de l'eau liquide jetee du panache ascendant a l'environnement

    INTEGER klab(klon, klev)
    REAL pen_u(klon, klev) ! quantite de masse entrainee pour panache ascendant
    REAL pde_u(klon, klev) ! quantite de masse detrainee pour panache ascendant
    REAL pen_d(klon, klev) ! quantite de masse entrainee pour panache descendant
    REAL pde_d(klon, klev) ! quantite de masse detrainee pour panache descendant

    ! Local:
    LOGICAL llflag(klon)
    INTEGER k, i, icall
    REAL zzs

    !----------------------------------------------------------------------

    ! Specify large scale parameters at half levels. Adjust
    ! temperature fields if statically unstable.

    DO k = 2, klev
       DO i = 1, klon
          pgeoh(i, k)=pgeo(i, k)+(pgeo(i, k-1)-pgeo(i, k))*0.5
          ptenh(i, k)=(MAX(RCPD*ten(i, k-1)+pgeo(i, k-1), &
               RCPD*ten(i, k)+pgeo(i, k))-pgeoh(i, k))/RCPD
          pqsenh(i, k)=pqsen(i, k-1)
          llflag(i)=.TRUE.
       ENDDO

       icall=0
       CALL flxadjtq(paph(1, k), ptenh(1, k), pqsenh(1, k), llflag, icall)

       DO i = 1, klon
          pqenh(i, k)=MIN(pqen(i, k-1), pqsen(i, k-1)) &
               +(pqsenh(i, k)-pqsen(i, k-1))
          pqenh(i, k)=MAX(pqenh(i, k), 0.)
       ENDDO
    end DO

    DO i = 1, klon
       ptenh(i, klev)=(RCPD*ten(i, klev)+pgeo(i, klev)- pgeoh(i, klev))/RCPD
       pqenh(i, klev)=pqen(i, klev)
       ptenh(i, 1)=ten(i, 1)
       pqenh(i, 1)=pqen(i, 1)
       pgeoh(i, 1)=pgeo(i, 1)
    end DO

    DO k = klev-1, 2, -1
       DO i = 1, klon
          zzs = MAX(RCPD*ptenh(i, k)+pgeoh(i, k), &
               RCPD*ptenh(i, k+1)+pgeoh(i, k+1))
          ptenh(i, k) = (zzs-pgeoh(i, k))/RCPD
       end DO
    end DO

    ! Initialize values for updrafts and downdrafts
    DO k = 1, klev
       DO i = 1, klon
          ptu(i, k) = ptenh(i, k)
          pqu(i, k) = pqenh(i, k)
          plu(i, k) = 0.
          mfu(i, k) = 0.
          mfus(i, k) = 0.
          mfuq(i, k) = 0.
          pdmfup(i, k) = 0.
          pdpmel(i, k) = 0.
          plude(i, k) = 0.
          klab(i, k) = 0
          ptd(i, k) = ptenh(i, k)
          pqd(i, k) = pqenh(i, k)
          mfd(i, k) = 0.0
          pmfds(i, k) = 0.0
          pmfdq(i, k) = 0.0
          pdmfdp(i, k) = 0.0
          pen_u(i, k) = 0.0
          pde_u(i, k) = 0.0
          pen_d(i, k) = 0.0
          pde_d(i, k) = 0.0
       ENDDO
    ENDDO

  END SUBROUTINE flxini

end module flxini_m
