module radiornpb_m

    IMPLICIT none

contains

  function radiornpb(tr_seri, pdtphys, tautr)

    ! From phylmd/radiornpb.F, v 1.2 2005/05/25 13:10:09

    ! Auteurs: AA + CG (LGGE/CNRS) Date 24-06-94
    ! Objet: Decroissance radioactive d'un traceur dans l'atmosphere
    !G 24 06 94 : Pour un traceur, le radon
    !G 16 12 94 : Plus un 2eme traceur, le 210Pb. Le radon decroit en plomb.

    ! Le pas de temps "pdtphys" est supposé beaucoup plus petit que la
    ! constante de temps de décroissance.

    use dimens_m, only: llm
    use dimphy, only: klon, nbtr
    use nr_util, only: assert

    REAL, intent(in):: tr_seri(:, :, :), pdtphys, tautr(:)
    real radiornpb(klon, llm, 2)

    ! Variable local to the procedure:
    INTEGER it

    !-----------------------------------------------

    call assert(shape(tr_seri) == (/klon, llm, nbtr/), "radiornpb tr_seri")
    call assert(size(tautr) == nbtr, "radiornpb tautr")

    DO it = 1, 2
       IF (tautr(it) > 0.) THEN
          radiornpb(:, :, it) = - tr_seri(:, :, it) * pdtphys / tautr(it)
       ELSE
          radiornpb(:, :, it) = 0.
       END IF
    END DO

    !G161294 : Cas particulier radon 1 => plomb 2
    radiornpb(:, :, 2) = radiornpb(:, :, 2) - radiornpb(:, :, 1)

  END function radiornpb

end module radiornpb_m
