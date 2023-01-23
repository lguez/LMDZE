module flumass_m

  IMPLICIT NONE

contains

  SUBROUTINE flumass(massebx, masseby, vcont, ucont, pbaru, pbarv)

    ! From LMDZ4/libf/dyn3d/flumass.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Auteurs : P. Le Van, F. Hourdin
    ! Objet: calcul du flux de masse aux niveaux s 

    USE comgeom, ONLY: aire, aireu
    USE dimensions, ONLY: iim, llm
    USE paramet_m, ONLY: iip1, iip2, ip1jm, ip1jmi1, ip1jmp1

    REAL, intent(in):: massebx(ip1jmp1, llm), masseby(ip1jm, llm) ! in kg
    real, intent(in):: vcont(ip1jm, llm), ucont(ip1jmp1, llm) ! in s-1

    ! Flux de masse :
    real, intent(out):: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm) ! in kg s-1

    ! Local:
    REAL apbarun(iip1), apbarus(iip1)
    REAL sairen, saireun, saires, saireus, ctn, cts, ctn0, cts0
    INTEGER l, ij, i

    !----------------------------------------------------------------

    DO l = 1, llm
       DO ij = iip2, ip1jm
          pbaru(ij, l) = massebx(ij, l) * ucont(ij, l)
       end DO

       DO ij = 1, ip1jm
          pbarv(ij, l) = masseby(ij, l) * vcont(ij, l)
       end DO
    end DO

    ! Calcul de la composante du flux de masse en x aux pôles 
    ! par la résolution d'un système de deux équations

    ! la première équation décrivant le calcul de la divergence en un
    ! point i du pôle, ce calcul étant itéré de i = 1 à i = im.
    ! c'est-à-dire,

    ! ((0.5 * pbaru(i) - 0.5 * pbaru(i - 1) - pbarv(i)) / aire(i) =
    ! - somme de (pbarv(n)) / aire pôle

    ! l'autre équation spécifiant que la moyenne du flux de masse au
    ! pôle est nulle c'est-à-dire somme de pbaru(n) * aire locale(n) =
    ! 0.

    ! on en revient ainsi à déterminer la constante additive commune
    ! aux pbaru qui représentait pbaru(0, j, l) dans l'équation du
    ! calcul de la divergence au point i=1.

    ! i variant de 1 à im
    ! n variant de 1 à im

    sairen = SUM(aire(:iim))
    saireun= SUM(aireu(:iim))
    saires = SUM(aire(ip1jm + 1: ip1jm + iim))
    saireus= SUM(aireu(ip1jm + 1: ip1jm + iim))

    DO l = 1, llm
       ctn = SUM(pbarv(:iim, l))/ sairen
       cts = SUM(pbarv(ip1jmi1 + 1: ip1jmi1 + iim, l)) / saires

       pbaru(1, l)= pbarv(1, l) - ctn * aire(1)
       pbaru(ip1jm+1, l)= - pbarv(ip1jmi1+1, l) + cts * aire(ip1jm+1)

       DO i = 2, iim
          pbaru(i, l) = pbaru(i - 1, l) + &
               pbarv(i, l) - ctn * aire(i)

          pbaru(i+ ip1jm, l) = pbaru(i+ ip1jm-1, l) - &
               pbarv(i+ ip1jmi1, l) + cts * aire(i+ ip1jm)
       end DO
       DO i = 1, iim
          apbarun(i) = aireu(i) * pbaru(i, l)
          apbarus(i) = aireu(i +ip1jm) * pbaru(i +ip1jm, l)
       end DO
       ctn0 = - SUM(apbarun(:iim)) / saireun
       cts0 = - SUM(apbarus(:iim)) / saireus
       DO i = 1, iim
          pbaru(i, l) = 2. * (pbaru(i, l) + ctn0)
          pbaru(i+ ip1jm, l) = 2. * (pbaru(i +ip1jm, l) + cts0)
       end DO

       pbaru(iip1, l) = pbaru(1, l)
       pbaru(ip1jmp1, l) = pbaru(ip1jm +1, l)
    end DO

  END SUBROUTINE flumass

end module flumass_m
