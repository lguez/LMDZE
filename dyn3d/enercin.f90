module enercin_m

  IMPLICIT NONE

contains

  SUBROUTINE enercin(vcov, ucov, vcont, ucont, ecin)

    ! From LMDZ4/libf/dyn3d/enercin.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Auteur: P. Le Van

    ! Objet: calcul de l'\'energie cin\'etique aux niveaux s

    USE comgeom, only: apoln, apols, alpha2p3, alpha1p4, alpha1p2, alpha3p4, &
         aire
    USE dimensions, only: iim, llm
    USE paramet_m, only: ip1jm, ip1jmp1, iip1, iip2, ip1jmi1

    REAL, INTENT(IN):: vcov(ip1jm, llm), ucov(ip1jmp1, llm)
    REAL, INTENT(IN):: vcont(ip1jm, llm), ucont(ip1jmp1, llm)
    real, INTENT(out):: ecin(ip1jmp1, llm)

    ! Local:
    REAL ecinni(iip1), ecinsi(iip1)
    REAL ecinpn, ecinps
    INTEGER l, ij, i
    REAL ssum

    !------------------------------------------------------------------------

    ! L'energie cinetique au point scalaire P(i, j), autre que les
    ! p\^oles, est :

    ! Ecin = 0.5 * U(i - 1, j)**2 * (alpha3 + alpha4) + 0.5 * U(i,
    ! j)**2 * (alpha1 + alpha2) + 0.5 * V(i, j - 1)**2 * (alpha1 +
    ! alpha4) + 0.5 * V(i, j)**2 * (alpha2 + alpha3)

    DO l = 1, llm
       DO ij = iip2, ip1jm - 1
          ecin(ij + 1, l) = 0.5 * (ucov(ij, l) * ucont(ij, l) &
               * alpha3p4(ij + 1) + ucov(ij + 1, l) * ucont(ij + 1, l) &
               * alpha1p2(ij + 1) + vcov(ij - iim, l) * vcont(ij - iim, l) &
               * alpha1p4 (ij + 1) + vcov(ij + 1, l) * vcont(ij + 1, l) &
               * alpha2p3(ij + 1))
       END DO

       ! correction pour ecin(1, j, l)
       ! ecin(1, j, l)= ecin(iip1, j, l)
       DO ij = iip2, ip1jm, iip1
          ecin(ij, l) = ecin(ij + iim, l)
       END DO

       ! calcul aux poles

       DO i = 1, iim
          ecinni(i) = vcov(i, l) * vcont(i, l) * aire(i)
          ecinsi(i) = vcov(i + ip1jmi1, l) * vcont(i + ip1jmi1, l) &
               * aire(i + ip1jm)
       END DO

       ecinpn = 0.5 * ssum(iim, ecinni, 1)/apoln
       ecinps = 0.5 * ssum(iim, ecinsi, 1)/apols

       DO ij = 1, iip1
          ecin(ij, l) = ecinpn
          ecin(ij + ip1jm, l) = ecinps
       END DO
    END DO

  END SUBROUTINE enercin

end module enercin_m
