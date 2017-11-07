module coefkzmin_m

  IMPLICIT NONE

contains

  SUBROUTINE coefkzmin(ngrid, ypaprs, ypplay, yu, yv, yt, yq, ycoefm, km, kn)

    ! From LMDZ4/libf/phylmd/coefkzmin.F, version 1.1.1.1 2004/05/19 12:53:08

    ! Entrées modifiées en attendant une version où les zlev et zlay
    ! soient disponibles.

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rd, rg, rkappa

    integer, intent(in):: ngrid
    REAL, intent(in):: ypaprs(klon, klev+1), ypplay(klon, klev)
    REAL, intent(in):: yu(klon, klev), yv(klon, klev) ! wind, in m s-1
    REAL, intent(in):: yt(klon, klev) ! temperature, in K
    REAL, intent(in):: yq(klon, klev)
    REAL, intent(in):: ycoefm(:) ! (ngrid) drag coefficient

    REAL, intent(inout):: km(:, 2:) ! (knon, 2:klev)
    ! coefficient de diffusion turbulente de quantité de mouvement (au
    ! bas de chaque couche) (en sortie : la valeur à la fin du pas de
    ! temps), m2 s-1

    REAL, intent(inout):: kn(:, 2:) ! (knon, 2:klev)
    ! coefficient de diffusion turbulente des scalaires (au bas de
    ! chaque couche) (en sortie : la valeur à la fin du pas de temps), m2 s-1

    ! Local:

    real ustar(ngrid) ! u*
    real zlay(ngrid, klev) ! in m
    integer i, k
    real pblhmin(ngrid)
    real, parameter:: coriol = 1e-4

    REAL zlev(ngrid, 2: klev)
    ! altitude at level (interface between layer with same index), in m

    REAL teta(ngrid, klev)
    ! température potentielle au centre de chaque couche (la valeur au
    ! debut du pas de temps)

    real, PARAMETER:: kap = 0.4

    !---------------------------------------------------------------------

    ! Debut de la partie qui doit etre incluse a terme dans clmain.

    do i = 1, ngrid
       zlay(i, 1) = RD * yt(i, 1) * 2 / (ypaprs(i, 1) + ypplay(i, 1)) &
            * (ypaprs(i, 1) - ypplay(i, 1)) / RG
    enddo

    do k = 2, klev
       do i = 1, ngrid
          zlay(i, k) = zlay(i, k-1) + RD * 0.5 * (yt(i, k - 1) + yt(i, k)) &
               / ypaprs(i, k) * (ypplay(i, k - 1) - ypplay(i, k)) / RG
       enddo
    enddo

    do k=1, klev
       do i = 1, ngrid
          ! Attention : on passe la temperature potentielle virtuelle
          ! pour le calcul de K.
          teta(i, k) = yt(i, k) * (ypaprs(i, 1) / ypplay(i, k))**rkappa &
               * (1. + 0.61 * yq(i, k))
       enddo
    enddo

    forall (k = 2: klev) zlev(:, k) = 0.5 * (zlay(:, k) + zlay(:, k-1))
    ustar = SQRT(ycoefm * (yu(:ngrid, 1)**2 + yv(:ngrid, 1)**2))

    ! Fin de la partie qui doit être incluse à terme dans clmain

    ! Cette routine est ecrite pour avoir en entree ustar, teta et zlev
    ! Ici, on a inclus le calcul de ces trois variables dans la routine
    ! coefkzmin en attendant une nouvelle version de la couche limite
    ! ou ces variables seront disponibles.

    ! Debut de la routine coefkzmin proprement dite

    pblhmin = 0.07 * ustar / coriol

    do k = 2, klev
       do i = 1, ngrid
          if (teta(i, 2) > teta(i, 1)) then
             kn(i, k) = kap * zlev(i, k) * ustar(i) &
                  * (max(1. - zlev(i, k) / pblhmin(i), 0.))**2
          else
             kn(i, k) = 0. ! min n'est utilisé que pour les SL stables
          endif
          km(i, k) = kn(i, k)
       enddo
    enddo

  end SUBROUTINE coefkzmin

end module coefkzmin_m
