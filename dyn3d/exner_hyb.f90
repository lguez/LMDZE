module exner_hyb_m

  IMPLICIT NONE

contains

  SUBROUTINE exner_hyb(ps, p, pks, pk)

    ! From dyn3d/exner_hyb.F, version 1.1.1.1, 2004/05/19 12:53:07
    ! Authors: P. Le Van, F. Hourdin

    ! Calcule la fonction d'Exner :
    ! pk = Cp * p ** kappa
    ! aux milieux des "llm" couches.
    ! "Pk(l)" est calcul\'e au milieu de la couche "l", entre les pressions
    ! "p(l)" et "p(l+1)", d\'efinies aux interfaces des couches.

    ! Au sommet de l'atmosph\`ere :
    ! p(llm+1) = 0
    ! "ps" et "pks" sont la pression et la fonction d'Exner au sol.

    ! \`A partir des relations :
    !(1) \overline{p * \delta_z pk}^z = kappa * pk * \delta_z p
    !(2) pk(l) = beta(l) * pk(l-1)
    ! (cf. documentation), on d\'etermine successivement, du haut vers
    ! le bas des couches, les coefficients : beta(llm), ..., beta(l), ...,
    ! beta(2) puis "pk(:, :, 1)". Ensuite, on calcule, du bas vers le
    ! haut des couches, "pk(:, :, l)" donn\'e par la relation (2), pour
    ! l = 2 \`a l = llm.

    use dimensions, only: llm
    use disvert_m, only: preff
    use suphec_m, only: rcpd, rkappa

    REAL, intent(in):: ps(:, :) ! (longitude, latitude)
    REAL, intent(in):: p(:, :, :) ! (longitude, latitude, llm + 1)

    real, intent(out):: pks(:, :) ! (longitude, latitude)
    real, intent(out):: pk(:, :, :) ! (longitude, latitude, llm)

    ! Variables locales :
    real beta(size(ps, 1), size(ps, 2), 2:llm)
    INTEGER l
    REAL unpl2k

    !-------------------------------------

    pks = rcpd * (ps / preff)**rkappa
    unpl2k = 1. + 2 * rkappa

    beta(:, :, llm) = 1. / unpl2k
    DO l = llm - 1, 2, -1
       beta(:, :, l) = p(:, :, l) &
            / (p(:, :, l) * unpl2k + p(:, :, l+1) * (beta(:, :, l+1) - unpl2k))
    ENDDO

    pk(:, :, 1) = ps * pks &
         / (ps * (1. + rkappa) + 0.5 * (beta(:, :, 2) - unpl2k) * p(:, :, 2))
    DO l = 2, llm
       pk(:, :, l) = beta(:, :, l) * pk(:, :, l - 1)
    ENDDO

  END SUBROUTINE exner_hyb

end module exner_hyb_m
