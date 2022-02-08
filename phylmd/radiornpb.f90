module radiornpb_m

    IMPLICIT none

contains

  function radiornpb(tr_seri, pdtphys, tautr)

    ! From phylmd/radiornpb.F, version 1.2, 2005/05/25 13:10:09

    ! Auteurs : AA + Christophe Genthon (LGGE/CNRS)
    ! Date: June 1994
    ! Objet: d\'ecroissance radioactive d'un traceur dans l'atmosph\`ere

    ! Pour un traceur, le radon. Plus un deuxi\`eme traceur, le 210
    ! Pb. Le radon d\'ecro\^it en plomb. Le pas de temps "pdtphys" est
    ! suppos\'e beaucoup plus petit que la constante de temps de
    ! d\'ecroissance.

    use dimensions, only: llm, nqmx
    use dimphy, only: klon
    use jumble, only: assert

    REAL, intent(in):: tr_seri(:, :, :) ! (klon, llm, nqmx - 2)
    REAL, intent(in):: pdtphys
    REAL, intent(in):: tautr(:) ! (nqmx - 2)
    real radiornpb(klon, llm, 2)

    ! Local:
    INTEGER it

    !-----------------------------------------------

    call assert(shape(tr_seri) == [klon, llm, nqmx - 2], "radiornpb tr_seri")
    call assert(size(tautr) == nqmx - 2, "radiornpb tautr")

    DO it = 1, 2
       IF (tautr(it) > 0.) THEN
          radiornpb(:, :, it) = - tr_seri(:, :, it) * pdtphys / tautr(it)
       ELSE
          radiornpb(:, :, it) = 0.
       END IF
    END DO

    ! Cas particulier radon 1 => plomb 2
    radiornpb(:, :, 2) = radiornpb(:, :, 2) - radiornpb(:, :, 1)

  END function radiornpb

end module radiornpb_m
