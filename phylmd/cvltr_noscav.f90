module cvltr_noscav_m

  IMPLICIT NONE

contains

  SUBROUTINE cvltr_noscav(pdtime, da, phi, mp, paprs, x, upd, dnd, dx)

    ! From LMDZ4/libf/phylmd/cvltr.F, version 1.1 2005/04/15 12:36:17

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rg

    ! Objet : convection des traceurs / Kerry Emanuel
    ! Authors: M.-A. Filiberti and J.-Y. Grandpeix

    REAL, intent(in):: pdtime
    real, intent(in):: da(klon, klev), phi(klon, klev, klev), mp(klon, klev)
    REAL, intent(in):: paprs(klon, klev + 1) ! pression aux 1/2 couches
    REAL, intent(in):: x(klon, klev) ! q de traceur (bas en haut)
    REAL, intent(in):: upd(klon, klev) ! saturated updraft mass flux
    REAL, intent(in):: dnd(klon, klev) ! saturated downdraft mass flux
    REAL, intent(out):: dx(klon, klev) ! tendance de traceur (bas en haut)

    ! Local:
    real zed(klon, klev), zmd(klon, klev, klev)
    real za(klon, klev, klev)
    real zmfd(klon, klev), zmfa(klon, klev)
    real zmfp(klon, klev), zmfu(klon, klev)
    integer i, k, j

    !------------------------------------------------------------

    ! calcul des tendances liees au downdraft

    zed = 0.
    zmfd = 0.
    zmfa = 0.
    zmfu = 0.
    zmfp = 0.
    zmd = 0.
    za = 0.
    ! entrainement
    do k = 1, klev - 1
       do i = 1, klon
          zed(i, k) = max(0., mp(i, k) - mp(i, k + 1))
       end do
    end do

    ! calcul de la matrice d echange
    ! matrice de distribution de la masse entrainee en k

    do k = 1, klev
       do i = 1, klon
          zmd(i, k, k) = zed(i, k)
       end do
    end do
    do k = 2, klev
       do j = k - 1, 1, - 1
          do i = 1, klon
             if(mp(i, j + 1) /= 0) then
                zmd(i, j, k) = zmd(i, j + 1, k) * min(1., mp(i, j)/mp(i, j + 1))
             endif
          end do
       end do
    end do
    do k = 1, klev
       do j = 1, klev - 1
          do i = 1, klon
             za(i, j, k) = max(0., zmd(i, j + 1, k) - zmd(i, j, k))
          end do
       end do
    end do

    ! rajout du terme lie a l'ascendance induite

    do j = 2, klev
       do i = 1, klon
          za(i, j, j - 1) = za(i, j, j - 1) + mp(i, j)
       end do
    end do

    ! tendances

    do k = 1, klev
       do j = 1, klev
          do i = 1, klon
             zmfd(i, j) = zmfd(i, j) + za(i, j, k) * (x(i, k) - x(i, j))
          end do
       end do
    end do

    ! calcul des tendances liees aux flux satures

    do j = 1, klev
       do i = 1, klon
          zmfa(i, j) = da(i, j) * (x(i, 1) - x(i, j))
       end do
    end do
    do k = 1, klev
       do j = 1, klev
          do i = 1, klon
             zmfp(i, j) = zmfp(i, j) + phi(i, j, k) * (x(i, k) - x(i, j))
          end do
       end do
    end do
    do j = 1, klev - 1
       do i = 1, klon
          zmfu(i, j) = max(0., upd(i, j + 1) + dnd(i, j + 1)) &
               * (x(i, j + 1) - x(i, j))
       end do
    end do
    do j = 2, klev
       do i = 1, klon
          zmfu(i, j) = zmfu(i, j) &
               + min(0., upd(i, j) + dnd(i, j)) * (x(i, j) - x(i, j - 1))
       end do
    end do

    ! calcul final des tendances

    do k = 1, klev
       do i = 1, klon
          dx(i, k) = (zmfd(i, k) + zmfu(i, k) &
               + zmfa(i, k) + zmfp(i, k)) * pdtime &
               * RG/(paprs(i, k) - paprs(i, k + 1))
       enddo
    enddo

  end SUBROUTINE cvltr_noscav

end module cvltr_noscav_m
