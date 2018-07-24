module coefkzmin_m

  IMPLICIT NONE

contains

  SUBROUTINE coefkzmin(ypaprs, ypplay, yu, yv, yt, yq, ycoefm, kn)

    ! From LMDZ4/libf/phylmd/coefkzmin.F, version 1.1.1.1, 2004/05/19 12:53:08

    ! Entrées modifiées en attendant une version où les zlev et zlay
    ! soient disponibles.

    USE dimphy, ONLY: klev
    USE suphec_m, ONLY: rd, rg, rkappa

    REAL, intent(in):: ypaprs(:, :) ! (knon, klev+1)
    REAL, intent(in):: ypplay(:, :) ! (knon, klev)
    REAL, intent(in):: yu(:, :), yv(:, :) ! (knon, klev) wind, in m s-1
    REAL, intent(in):: yt(:, :) ! (knon, klev) temperature, in K
    REAL, intent(in):: yq(:, :) ! (knon, klev)
    REAL, intent(in):: ycoefm(:) ! (knon) drag coefficient

    REAL, intent(out):: kn(:, 2:) ! (knon, 2:klev) coefficient de
    ! diffusion turbulente de la quantité de mouvement et des
    ! scalaires (au bas de chaque couche) (en sortie : la valeur à la
    ! fin du pas de temps), m2 s-1

    ! Local:

    integer knon
    real ustar(size(ypaprs, 1)) ! (knon) u*
    real zlay(size(ypaprs, 1), klev) ! (knon, klev) in m
    integer i, k
    real pblhmin(size(ypaprs, 1)) ! (knon)
    real, parameter:: coriol = 1e-4

    REAL zlev(size(ypaprs, 1), 2: klev) ! (knon, 2: klev)
    ! altitude at level (interface between layer with same index), in m

    REAL teta(size(ypaprs, 1), klev) ! (knon, klev)
    ! température potentielle au centre de chaque couche (la valeur au
    ! debut du pas de temps)

    real, PARAMETER:: kap = 0.4

    !---------------------------------------------------------------------

    knon = size(ypaprs, 1)
    
    ! Debut de la partie qui doit etre incluse a terme dans pbl_surface.

    do i = 1, knon
       zlay(i, 1) = RD * yt(i, 1) * 2 / (ypaprs(i, 1) + ypplay(i, 1)) &
            * (ypaprs(i, 1) - ypplay(i, 1)) / RG
    enddo

    do k = 2, klev
       do i = 1, knon
          zlay(i, k) = zlay(i, k-1) + RD * 0.5 * (yt(i, k - 1) + yt(i, k)) &
               / ypaprs(i, k) * (ypplay(i, k - 1) - ypplay(i, k)) / RG
       enddo
    enddo

    do k=1, klev
       do i = 1, knon
          ! Attention : on passe la temperature potentielle virtuelle
          ! pour le calcul de K.
          teta(i, k) = yt(i, k) * (ypaprs(i, 1) / ypplay(i, k))**rkappa &
               * (1. + 0.61 * yq(i, k))
       enddo
    enddo

    forall (k = 2: klev) zlev(:, k) = 0.5 * (zlay(:, k) + zlay(:, k-1))
    ustar = SQRT(ycoefm * (yu(:, 1)**2 + yv(:, 1)**2))

    ! Fin de la partie qui doit être incluse à terme dans pbl_surface

    ! Cette routine est ecrite pour avoir en entree ustar, teta et zlev
    ! Ici, on a inclus le calcul de ces trois variables dans la routine
    ! coefkzmin en attendant une nouvelle version de la couche limite
    ! ou ces variables seront disponibles.

    ! Debut de la routine coefkzmin proprement dite

    pblhmin = 0.07 * ustar / coriol

    do k = 2, klev
       do i = 1, knon
          if (teta(i, 2) > teta(i, 1)) then
             kn(i, k) = kap * zlev(i, k) * ustar(i) &
                  * (max(1. - zlev(i, k) / pblhmin(i), 0.))**2
          else
             kn(i, k) = 0. ! min n'est utilisé que pour les SL stables
          endif
       enddo
    enddo

  end SUBROUTINE coefkzmin

end module coefkzmin_m
