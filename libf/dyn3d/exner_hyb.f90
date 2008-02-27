module exner_hyb_m

  IMPLICIT NONE

contains

  SUBROUTINE exner_hyb(ps, p, pks, pk, pkf)

    ! From dyn3d/exner_hyb.F, v 1.1.1.1 2004/05/19 12:53:07

    ! Auteurs :  P. Le Van, F. Hourdin.

    ! Calcule la fonction d'Exner :
    ! pk = Cp * p ** kappa
    ! aux milieux des couches.
    ! Pk(l) est calculé aux milieux des couches "l", entre les pressions
    ! "p(l)" et "p(l+1)", définies aux interfaces des "llm" couches.

    ! Au sommet de l'atmosphère :
    ! p(llm+1) = 0.
    ! et ps et pks sont la pression et la fonction d'Exner au sol.

    ! À partir des relations :

    !   -------- z
    !(1) p*dz(pk) = kappa *pk*dz(p)

    !(2) pk(l) = alpha(l)+ beta(l)*pk(l-1)

    ! (voir note de F. Hourdin), on determine successivement, du haut
    ! vers le bas des couches, les coefficients :
    ! alpha(llm), beta(llm)..., alpha(l), beta(l)..., alpha(2), beta(2)
    ! puis "pk(ij, 1)".
    ! Ensuite, on calcule, du bas vers le haut des couches, "pk(ij, l)"
    ! donné par la relation (2), pour l = 2 à l = llm.

    use dimens_m, only: iim, jjm, llm
    use comconst, only: kappa, cpp
    use comvert, only: preff
    use comgeom, only: aire_2d, apoln, apols

    REAL, intent(in):: ps((iim + 1) * (jjm + 1))
    REAL, intent(in):: p((iim + 1) * (jjm + 1), llm + 1)

    real, intent(out):: pks((iim + 1) * (jjm + 1))
    real, intent(out):: pk((iim + 1) * (jjm + 1), llm)
    real, intent(out), optional:: pkf((iim + 1) * (jjm + 1), llm)

    ! Variables locales

    real alpha((iim + 1) * (jjm + 1), llm), beta((iim + 1) * (jjm + 1), llm)
    INTEGER l, ij
    REAL unpl2k, dellta

    REAL ppn(iim), pps(iim)
    REAL xpn, xps
    REAL SSUM

    !-------------------------------------

    pks(:) = cpp * (ps(:) / preff)**kappa
    ppn(:) = aire_2d(:iim, 1) * pks(:iim)
    pps(:) = aire_2d(:iim, jjm + 1) &
         * pks(1 + (iim + 1) * jjm: iim + (iim + 1) * jjm)
    xpn = SSUM(iim, ppn, 1) /apoln
    xps = SSUM(iim, pps, 1) /apols
    pks(:iim + 1) = xpn
    pks(1+(iim + 1) * jjm:) = xps

    unpl2k = 1. + 2 * kappa

    ! Calcul des coeff. alpha et beta  pour la couche l = llm
    alpha(:, llm) = 0.
    beta (:, llm) = 1./ unpl2k

    ! Calcul des coeff. alpha et beta  pour l = llm-1  à l = 2
    DO l = llm -1 , 2 , -1
       DO ij = 1, (iim + 1) * (jjm + 1)
          dellta = p(ij, l)* unpl2k + p(ij, l+1)* ( beta(ij, l+1)-unpl2k )
          alpha(ij, l)  = - p(ij, l+1) / dellta * alpha(ij, l+1)
          beta (ij, l)  =   p(ij, l  ) / dellta   
       ENDDO
    ENDDO

    ! Calcul de pk pour la couche 1, près du sol :
    pk(:, 1) = (p(:, 1) * pks(:) - 0.5 * alpha(:, 2) * p(:, 2))  &
         / (p(:, 1) * (1. + kappa) + 0.5 * (beta(:, 2) - unpl2k) * p(:, 2))

    ! Calcul de pk(ij, l) , pour l = 2 à l = llm
    DO l = 2, llm
       DO   ij   = 1, (iim + 1) * (jjm + 1)
          pk(ij, l) = alpha(ij, l) + beta(ij, l) * pk(ij, l-1)
       ENDDO
    ENDDO

    if (present(pkf)) then
       pkf(:, :) = pk(:, :)
       CALL filtreg(pkf, jjm + 1, llm, 2, 1, .TRUE., 1)
    end if

  END SUBROUTINE exner_hyb

end module exner_hyb_m
