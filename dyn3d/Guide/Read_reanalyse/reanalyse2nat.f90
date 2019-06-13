module reanalyse2nat_m

  implicit none

contains

  subroutine reanalyse2nat(invert_y, psi, unc, vnc, tnc, qnc, pl, u, v, t, q, &
       pk)

    ! Inversion nord-sud de la grille et interpolation verticale sur
    ! les niveaux du modèle.

    USE comconst, ONLY: cpp, kappa
    USE comgeom, ONLY: aireu_2d, airev_2d, aire_2d
    USE dimensions, ONLY: jjm, llm
    USE disvert_m, ONLY: ap, bp, preff
    USE exner_hyb_m, ONLY: exner_hyb
    use massbar_m, only: massbar
    USE paramet_m, ONLY: iip1, jjp1, llmp1
    use pres2lev_m, only: pres2lev

    logical, intent(in):: invert_y
    real, intent(in):: psi(:, :) ! (iip1, jjp1)

    real, intent(in):: unc(:, :, :) ! (iip1, jjp1, :)
    real, intent(in):: vnc(:, :, :) ! (iip1, jjm, :)
    real, intent(in):: tnc(:, :, :) ! (iip1, jjp1, :)
    real, intent(in):: qnc(:, :, :) ! (iip1, jjp1, :)
    real, intent(in):: pl(:)

    real, intent(out):: u(:, :, :) ! (iip1, jjp1, llm)
    real, intent(out):: v(:, :, :) ! (iip1, jjm, llm)
    real, intent(out):: t(:, :, :), q(:, :, :) ! (iip1, jjp1, llm)
    real, intent(out):: pk(:, :, :) ! (iip1, jjp1, llm)

    ! Local:

    real zu(iip1, jjp1, llm), zv(iip1, jjm, llm)
    real zt(iip1, jjp1, llm), zq(iip1, jjp1, llm)

    real pext(iip1, jjp1, llm)
    real pbarx(iip1, jjp1, llm), pbary(iip1, jjm, llm)
    real plunc(iip1, jjp1, llm), plvnc(iip1, jjm, llm)
    real plsnc(iip1, jjp1, llm)

    real p(iip1, jjp1, llmp1)
    real pks(iip1, jjp1)
    real pls(iip1, jjp1, llm)
    real unskap

    integer i, j, l

    ! -----------------------------------------------------------------

    ! calcul de la pression au milieu des couches
    forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * psi
    CALL exner_hyb(psi, p, pks, pk)

    ! Calcul de pls, pression au milieu des couches, en Pascals
    unskap=1./kappa
    DO l = 1, llm
       DO j=1, jjp1
          DO i =1, iip1
             pls(i, j, l) = preff * ( pk(i, j, l)/cpp) ** unskap
          ENDDO
       ENDDO
    ENDDO

    ! calcul des pressions pour les grilles u et v

    do l=1, llm
       do j=1, jjp1
          do i=1, iip1
             pext(i, j, l)=pls(i, j, l)*aire_2d(i, j)
          enddo
       enddo
    enddo
    call massbar(pext, pbarx, pbary )
    do l=1, llm
       do j=1, jjp1
          do i=1, iip1
             plunc(i, jjp1+1-j, l)=pbarx(i, j, l)/aireu_2d(i, j)
             plsnc(i, jjp1+1-j, l)=pls(i, j, l)
          enddo
       enddo
    enddo
    do l=1, llm
       do j=1, jjm
          do i=1, iip1
             plvnc(i, jjm+1-j, l)=pbary(i, j, l)/airev_2d(i, j)
          enddo
       enddo
    enddo

    call pres2lev(unc, zu, pl, plunc)
    call pres2lev(vnc, zv, pl, plvnc )
    call pres2lev(tnc, zt, pl, plsnc)
    call pres2lev(qnc, zq, pl, plsnc)

    if (invert_y) then
       ! Inversion Nord/Sud
       u=zu(:, jjp1:1:-1, :)
       v=zv(:, jjm:1:-1, :)
       t=zt(:, jjp1:1:-1, :)
       q=zq(:, jjp1:1:-1, :)
    else
       u = zu
       v = zv
       t = zt
       q = zq
    end if

  end subroutine reanalyse2nat

end module reanalyse2nat_m
