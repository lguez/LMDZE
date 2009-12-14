module caldyn0_m

  IMPLICIT NONE

contains

  SUBROUTINE caldyn0(ucov, vcov, teta, ps, masse, pk, phis, phi, w, &
       pbaru, pbarv)

    ! From dyn3d/caldyn0.F, v 1.1.1.1 2004/05/19 12:53:07
    ! Auteur :  P. Le Van
    ! Objet: calcul des tendances dynamiques
    ! Modif 04/93 F.Forget

    USE dimens_m, ONLY : llm
    USE paramet_m, ONLY : iip1, ip1jm, ip1jmp1, jjp1, llmp1
    USE comvert, ONLY : ap, bp
    USE comgeom, ONLY : airesurg
    USE pression_m, ONLY : pression

    !   Arguments:
    REAL, INTENT (IN) :: vcov(ip1jm, llm), ucov(ip1jmp1, llm)
    REAL :: teta(ip1jmp1, llm)
    REAL, INTENT (IN) :: ps(ip1jmp1)
    REAL, INTENT (IN) :: phis(ip1jmp1)
    REAL, INTENT (IN) :: pk(iip1, jjp1, llm)
    REAL :: vcont(ip1jm, llm), ucont(ip1jmp1, llm)
    REAL :: phi(ip1jmp1, llm), masse(ip1jmp1, llm)
    REAL :: pbaru(ip1jmp1, llm), pbarv(ip1jm, llm)

    !   Local:

    REAL :: p(ip1jmp1, llmp1)
    REAL :: massebx(ip1jmp1, llm), masseby(ip1jm, llm), psexbarxy(ip1jm)
    REAL :: vorpot(ip1jm, llm)
    REAL :: w(ip1jmp1, llm), ecin(ip1jmp1, llm), convm(ip1jmp1, llm)
    REAL :: bern(ip1jmp1, llm)
    REAL :: massebxy(ip1jm, llm), dp(ip1jmp1)

    INTEGER :: ij

    !-----------------------------------------------------------------------

    PRINT *, 'Call sequence information: caldyn0'

    !   Calcul des tendances dynamiques:

    CALL covcont(llm, ucov, vcov, ucont, vcont)
    CALL pression(ip1jmp1, ap, bp, ps, p)
    CALL psextbar(ps, psexbarxy)
    CALL massdair(p, masse)
    CALL massbar(masse, massebx, masseby)
    CALL massbarxy(masse, massebxy)
    CALL flumass(massebx, masseby, vcont, ucont, pbaru, pbarv)
    CALL convmas(pbaru, pbarv, convm)

    DO ij = 1, ip1jmp1
       dp(ij) = convm(ij, 1)/airesurg(ij)
    END DO

    CALL vitvert(convm, w)
    CALL tourpot(vcov, ucov, massebxy, vorpot)
    CALL enercin(vcov, ucov, vcont, ucont, ecin)
    CALL bernoui(ip1jmp1, llm, phi, ecin, bern)
    CALL sortvarc0(ucov, teta, ps, masse, pk, phis, vorpot, phi, bern, dp)

  END SUBROUTINE caldyn0

end module caldyn0_m
