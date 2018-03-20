module lwbv_m

  IMPLICIT NONE

contains

  SUBROUTINE lwbv(klim, pdt0, pemis, ppmb, ptl, ptave, pabcu, pfluc, &
       pbint, pbsui, pcts, pcntrb)

    USE dimensions
    USE dimphy
    use lwv_m, only: lwv
    USE suphec_m
    USE raddim
    USE raddimlw

    ! PURPOSE.
    ! --------
    ! TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
    ! VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
    ! SAVING

    ! METHOD.
    ! -------

    ! 1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
    ! GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
    ! 2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
    ! TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
    ! BOUNDARIES.
    ! 3. COMPUTES THE CLEAR-SKY COOLING RATES.

    ! REFERENCE.
    ! ----------

    ! SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    ! ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    ! AUTHOR.
    ! -------
    ! JEAN-JACQUES MORCRETTE  *ECMWF*

    ! MODIFICATIONS.
    ! --------------
    ! ORIGINAL : 89-07-14
    ! MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
    ! MEMORY)
    ! -----------------------------------------------------------------------
    ! * ARGUMENTS:
    INTEGER klim

    DOUBLE PRECISION pdt0(kdlon)
    DOUBLE PRECISION pemis(kdlon)
    DOUBLE PRECISION ppmb(kdlon, kflev+1)
    DOUBLE PRECISION ptl(kdlon, kflev+1)
    DOUBLE PRECISION ptave(kdlon, kflev)

    DOUBLE PRECISION pfluc(kdlon, 2, kflev+1)

    DOUBLE PRECISION pabcu(kdlon, nua, 3*kflev+1)
    DOUBLE PRECISION pbint(kdlon, kflev+1)
    DOUBLE PRECISION pbsui(kdlon)
    DOUBLE PRECISION pcts(kdlon, kflev)
    DOUBLE PRECISION pcntrb(kdlon, kflev+1, kflev+1)

    ! -------------------------------------------------------------------------

    ! * LOCAL VARIABLES:
    DOUBLE PRECISION zb(kdlon, ninter, kflev+1)
    DOUBLE PRECISION zbsur(kdlon, ninter)
    DOUBLE PRECISION zbtop(kdlon, ninter)
    DOUBLE PRECISION zdbsl(kdlon, ninter, kflev*2)
    DOUBLE PRECISION zga(kdlon, 8, 2, kflev)
    DOUBLE PRECISION zgb(kdlon, 8, 2, kflev)
    DOUBLE PRECISION zgasur(kdlon, 8, 2)
    DOUBLE PRECISION zgbsur(kdlon, 8, 2)
    DOUBLE PRECISION zgatop(kdlon, 8, 2)
    DOUBLE PRECISION zgbtop(kdlon, 8, 2)

    INTEGER nuaer, ntraer
    ! ------------------------------------------------------------------
    ! * COMPUTES PLANCK FUNCTIONS:
    CALL lwb(pdt0, ptave, ptl, zb, pbint, pbsui, zbsur, zbtop, zdbsl, zga, &
         zgb, zgasur, zgbsur, zgatop, zgbtop)
    ! ------------------------------------------------------------------
    ! * PERFORMS THE VERTICAL INTEGRATION:
    nuaer = nua
    ntraer = ntra
    CALL lwv(nuaer, ntraer, klim, pabcu, zb, pbint, pbsui, zbsur, zbtop, &
         zdbsl, pemis, ppmb, zga, zgb, zgasur, zgbsur, zgatop, zgbtop, pcntrb, &
         pcts, pfluc)

  END SUBROUTINE lwbv

end module lwbv_m
