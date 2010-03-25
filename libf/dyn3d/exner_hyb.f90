module exner_hyb_m

  IMPLICIT NONE

contains

  SUBROUTINE exner_hyb(ps, p, pks, pk, pkf)

    ! From dyn3d/exner_hyb.F, v 1.1.1.1 2004/05/19 12:53:07

    ! Authors : P. Le Van, F. Hourdin

    ! Calcule la fonction d'Exner :
    ! pk = Cp * p ** kappa
    ! aux milieux des "llm" couches.
    ! "Pk(l)" est calcul� au milieu de la couche "l", entre les pressions
    ! "p(l)" et "p(l+1)", d�finies aux interfaces des couches.

    ! Au sommet de l'atmosph�re :
    ! p(llm+1) = 0.
    ! "ps" et "pks" sont la pression et la fonction d'Exner au sol.

    ! � partir des relations :

    !   -------- z
    !(1) p*dz(pk) = kappa * pk * dz(p)

    !(2) pk(l) = alpha(l)+ beta(l) * pk(l-1)

    ! (voir note de F. Hourdin), on d�termine successivement, du haut
    ! vers le bas des couches, les coefficients :
    ! alpha(llm), beta(llm)..., alpha(l), beta(l)..., alpha(2), beta(2)
    ! puis "pk(:, 1)".
    ! Ensuite, on calcule, du bas vers le haut des couches, "pk(:, l)"
    ! donn� par la relation (2), pour l = 2 � l = llm.

    use dimens_m, only: iim, jjm, llm
    use comconst, only: kappa, cpp
    use comvert, only: preff
    use comgeom, only: aire_2d, apoln, apols
    use filtreg_m, only: filtreg

    REAL, intent(in):: ps((iim + 1) * (jjm + 1))
    REAL, intent(in):: p((iim + 1) * (jjm + 1), llm + 1)

    real, intent(out):: pks((iim + 1) * (jjm + 1))
    real, intent(out):: pk((iim + 1) * (jjm + 1), llm)
    real, intent(out), optional:: pkf((iim + 1) * (jjm + 1), llm)

    ! Variables locales

    real alpha((iim + 1) * (jjm + 1), llm), beta((iim + 1) * (jjm + 1), llm)
    INTEGER l
    REAL unpl2k, dellta((iim + 1) * (jjm + 1))

    REAL ppn(iim), pps(iim)

    !-------------------------------------

    pks = cpp * (ps / preff)**kappa
    ppn = aire_2d(:iim, 1) * pks(:iim)
    pps = aire_2d(:iim, jjm + 1) &
         * pks(1 + (iim + 1) * jjm: iim + (iim + 1) * jjm)
    pks(:iim + 1) = SUM(ppn) /apoln
    pks(1+(iim + 1) * jjm:) = SUM(pps) /apols

    unpl2k = 1. + 2 * kappa

    ! Calcul des coefficients alpha et beta pour la couche l = llm :
    alpha(:, llm) = 0.
    beta(:, llm) = 1./ unpl2k

    ! Calcul des coefficients alpha et beta pour l = llm-1 � l = 2 :
    DO l = llm - 1, 2, -1
       dellta = p(:, l) * unpl2k + p(:, l+1) * (beta(:, l+1) - unpl2k)
       alpha(:, l) = - p(:, l+1) / dellta * alpha(:, l+1)
       beta(:, l) = p(:, l) / dellta   
    ENDDO

    ! Calcul de pk pour la couche 1, pr�s du sol :
    pk(:, 1) = (p(:, 1) * pks - 0.5 * alpha(:, 2) * p(:, 2))  &
         / (p(:, 1) * (1. + kappa) + 0.5 * (beta(:, 2) - unpl2k) * p(:, 2))

    ! Calcul de pk(:, l) pour l = 2 � l = llm :
    DO l = 2, llm
       pk(:, l) = alpha(:, l) + beta(:, l) * pk(:, l-1)
    ENDDO

    if (present(pkf)) then
       pkf = pk
       CALL filtreg(pkf, jjm + 1, llm, 2, 1, .TRUE., 1)
    end if

  END SUBROUTINE exner_hyb

end module exner_hyb_m
