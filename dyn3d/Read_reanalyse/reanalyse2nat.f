module reanalyse2nat_m

  implicit none

contains

  subroutine reanalyse2nat(nlevnc, psi, unc, vnc, tnc, qnc, pl, u, v, t, q, &
       masse, pk)

    ! Inversion nord-sud de la grille et interpolation sur les niveaux
    ! verticaux du modèle.

    USE dimens_m, ONLY: iim, jjm, llm
    USE paramet_m, ONLY: iip1, jjp1, llmp1
    USE comconst, ONLY: cpp, kappa
    USE disvert_m, ONLY: ap, bp, preff
    USE comgeom, ONLY: aireu_2d, airev_2d, aire_2d
    USE exner_hyb_m, ONLY: exner_hyb
    use massdair_m, only: massdair

    integer nlevnc
    real, intent(in):: psi(iip1, jjp1)
    real unc(iip1, jjp1, nlevnc), vnc(iip1, jjm, nlevnc)
    real tnc(iip1, jjp1, nlevnc)
    real qnc(iip1, jjp1, nlevnc)
    real pl(nlevnc)
    real u(iip1, jjp1, llm), v(iip1, jjm, llm)
    real t(iip1, jjp1, llm), q(iip1, jjp1, llm)
    real masse(iip1, jjp1, llm)
    real pk(iip1, jjp1, llm)

    ! Local:

    real zu(iip1, jjp1, llm), zv(iip1, jjm, llm)
    real zt(iip1, jjp1, llm), zq(iip1, jjp1, llm)

    real pext(iip1, jjp1, llm)
    real pbarx(iip1, jjp1, llm), pbary(iip1, jjm, llm)
    real plunc(iip1, jjp1, llm), plvnc(iip1, jjm, llm)
    real plsnc(iip1, jjp1, llm)

    real p(iip1, jjp1, llmp1)
    real pks(iip1, jjp1)
    real pkf(iip1, jjp1, llm)
    real pls(iip1, jjp1, llm)
    real prefkap, unskap

    integer i, j, l

    ! -----------------------------------------------------------------

    ! calcul de la pression au milieu des couches
    forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * psi
    call massdair(p, masse)
    CALL exner_hyb(psi, p, pks, pk, pkf)

    ! Calcul de pls, pression au milieu des couches, en Pascals
    unskap=1./kappa
    prefkap = preff ** kappa
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

    call pres2lev(unc, zu, nlevnc, llm, pl, plunc, iip1, jjp1)
    call pres2lev(vnc, zv, nlevnc, llm, pl, plvnc, iip1, jjm )
    call pres2lev(tnc, zt, nlevnc, llm, pl, plsnc, iip1, jjp1)
    call pres2lev(qnc, zq, nlevnc, llm, pl, plsnc, iip1, jjp1)

    ! Inversion Nord/Sud
    do l=1, llm
       do j=1, jjp1
          do i=1, iim
             u(i, j, l)=zu(i, jjp1+1-j, l)
             t(i, j, l)=zt(i, jjp1+1-j, l)
             q(i, j, l)=zq(i, jjp1+1-j, l)
          enddo
          u(iip1, j, l)=u(1, j, l)
          t(iip1, j, l)=t(1, j, l)
          q(iip1, j, l)=q(1, j, l)
       enddo
    enddo

    do l=1, llm
       do j=1, jjm
          do i=1, iim
             v(i, j, l)=zv(i, jjm+1-j, l)
          enddo
          v(iip1, j, l)=v(1, j, l)
       enddo
    enddo

  end subroutine reanalyse2nat

end module reanalyse2nat_m
