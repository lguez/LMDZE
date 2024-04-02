module lwbv_m

  IMPLICIT NONE

contains

  SUBROUTINE lwbv(klim, pdt0, pemis, ppmb, ptl, ptave, pabcu, pfluc, pbint, &
       pbsui, pcts, pcntrb)

    USE dimensions, only: llm
    USE dimphy, only: klon
    use lwv_m, only: lwv
    USE raddimlw, only: nua, ninter, ntra

    ! PURPOSE. TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
    ! VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY SAVING

    ! METHOD.
    ! 1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
    ! GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
    ! 2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
    ! TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
    ! BOUNDARIES.
    ! 3. COMPUTES THE CLEAR-SKY COOLING RATES.

    ! REFERENCE. SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    ! AUTHOR. JEAN-JACQUES MORCRETTE *ECMWF*

    ! MODIFICATIONS.
    ! ORIGINAL : 89-07-14
    ! MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
    ! MEMORY)

    ! ARGUMENTS:
    INTEGER klim
    DOUBLE PRECISION pdt0(klon)
    DOUBLE PRECISION pemis(klon)
    DOUBLE PRECISION ppmb(klon, llm+1)
    DOUBLE PRECISION ptl(klon, llm+1)
    DOUBLE PRECISION ptave(klon, llm)
    DOUBLE PRECISION pfluc(klon, 2, llm+1)
    DOUBLE PRECISION pabcu(klon, nua, 3*llm+1)
    DOUBLE PRECISION pbint(klon, llm+1)
    DOUBLE PRECISION pbsui(klon)
    DOUBLE PRECISION pcts(klon, llm)
    DOUBLE PRECISION pcntrb(klon, llm+1, llm+1)

    ! LOCAL VARIABLES:
    DOUBLE PRECISION zb(klon, ninter, llm+1)
    DOUBLE PRECISION zbsur(klon, ninter)
    DOUBLE PRECISION zbtop(klon, ninter)
    DOUBLE PRECISION zdbsl(klon, ninter, llm*2)
    DOUBLE PRECISION zga(klon, 8, 2, llm)
    DOUBLE PRECISION zgb(klon, 8, 2, llm)
    DOUBLE PRECISION zgasur(klon, 8, 2)
    DOUBLE PRECISION zgbsur(klon, 8, 2)
    DOUBLE PRECISION zgatop(klon, 8, 2)
    DOUBLE PRECISION zgbtop(klon, 8, 2)
    INTEGER nuaer, ntraer

    !------------------------------------------------------------------

    ! COMPUTES PLANCK FUNCTIONS:
    CALL lwb(pdt0, ptave, ptl, zb, pbint, pbsui, zbsur, zbtop, zdbsl, zga, &
         zgb, zgasur, zgbsur, zgatop, zgbtop)

    ! PERFORMS THE VERTICAL INTEGRATION:
    nuaer = nua
    ntraer = ntra
    CALL lwv(nuaer, ntraer, klim, pabcu, zb, pbint, pbsui, zbsur, zbtop, &
         zdbsl, pemis, ppmb, zga, zgb, zgasur, zgbsur, zgatop, zgbtop, pcntrb, &
         pcts, pfluc)

  END SUBROUTINE lwbv

end module lwbv_m
