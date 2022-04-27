module lwvn_m

  IMPLICIT NONE

contains

  SUBROUTINE lwvn(kuaer, pabcu, pdbsl, pga, pgb, padjd, padju, pcntrb, pdbdt)
    USE dimensions
    USE dimphy
    USE raddimlw
    ! -----------------------------------------------------------------------
    ! PURPOSE.
    ! --------
    ! CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
    ! TO GIVE LONGWAVE FLUXES OR RADIANCES

    ! METHOD.
    ! -------

    ! 1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
    ! CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE

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
    ! -----------------------------------------------------------------------

    ! * ARGUMENTS:

    INTEGER kuaer

    DOUBLE PRECISION pabcu(klon, nua, 3*llm+1) ! ABSORBER AMOUNTS
    DOUBLE PRECISION pdbsl(klon, ninter, llm*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
    DOUBLE PRECISION pga(klon, 8, 2, llm) ! PADE APPROXIMANTS
    DOUBLE PRECISION pgb(klon, 8, 2, llm) ! PADE APPROXIMANTS

    DOUBLE PRECISION padjd(klon, llm+1) ! CONTRIBUTION OF ADJACENT LAYERS
    DOUBLE PRECISION padju(klon, llm+1) ! CONTRIBUTION OF ADJACENT LAYERS
    DOUBLE PRECISION pcntrb(klon, llm+1, llm+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
    DOUBLE PRECISION pdbdt(klon, ninter, llm) !  LAYER PLANCK FUNCTION GRADIENT

    ! * LOCAL ARRAYS:

    DOUBLE PRECISION zglayd(klon)
    DOUBLE PRECISION zglayu(klon)
    DOUBLE PRECISION ztt(klon, ntra)
    DOUBLE PRECISION zuu(klon, nua)

    INTEGER jk, jl, ja, im12, ind, inu, ixu, jg
    INTEGER ixd, ibs, idd, imu, jk1, jk2, jnu
    DOUBLE PRECISION zwtr

    ! * Data Block:

    DOUBLE PRECISION wg1(2)
    SAVE wg1
    DATA (wg1(jk), jk=1, 2)/1d0, 1d0/
    ! -----------------------------------------------------------------------

    ! *         1.    INITIALIZATION
    ! --------------


    ! *         1.1     INITIALIZE LAYER CONTRIBUTIONS
    ! ------------------------------


    DO jk = 1, llm + 1
       DO jl = 1, klon
          padjd(jl, jk) = 0.
          padju(jl, jk) = 0.
       END DO
    END DO

    ! *         1.2     INITIALIZE TRANSMISSION FUNCTIONS
    ! ---------------------------------


    DO ja = 1, ntra
       DO jl = 1, klon
          ztt(jl, ja) = 1.0
       END DO
    END DO

    DO ja = 1, nua
       DO jl = 1, klon
          zuu(jl, ja) = 0.
       END DO
    END DO

    ! ------------------------------------------------------------------

    ! *         2.      VERTICAL INTEGRATION
    ! --------------------



    ! *         2.1     CONTRIBUTION FROM ADJACENT LAYERS
    ! ---------------------------------


    DO jk = 1, llm

       ! *         2.1.1   DOWNWARD LAYERS
       ! ---------------


       im12 = 2*(jk-1)
       ind = (jk-1)*ng1p1 + 1
       ixd = ind
       inu = jk*ng1p1 + 1
       ixu = ind

       DO jl = 1, klon
          zglayd(jl) = 0.
          zglayu(jl) = 0.
       END DO

       DO jg = 1, ng1
          ibs = im12 + jg
          idd = ixd + jg
          DO ja = 1, kuaer
             DO jl = 1, klon
                zuu(jl, ja) = pabcu(jl, ja, ind) - pabcu(jl, ja, idd)
             END DO
          END DO


          CALL lwtt(pga(1,1,1,jk), pgb(1,1,1,jk), zuu, ztt)

          DO jl = 1, klon
             zwtr = pdbsl(jl, 1, ibs)*ztt(jl, 1)*ztt(jl, 10) + &
                  pdbsl(jl, 2, ibs)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
                  pdbsl(jl, 3, ibs)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
                  pdbsl(jl, 4, ibs)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
                  pdbsl(jl, 5, ibs)*ztt(jl, 3)*ztt(jl, 14) + &
                  pdbsl(jl, 6, ibs)*ztt(jl, 6)*ztt(jl, 15)
             zglayd(jl) = zglayd(jl) + zwtr*wg1(jg)
          END DO

          ! *         2.1.2   DOWNWARD LAYERS
          ! ---------------


          imu = ixu + jg
          DO ja = 1, kuaer
             DO jl = 1, klon
                zuu(jl, ja) = pabcu(jl, ja, imu) - pabcu(jl, ja, inu)
             END DO
          END DO


          CALL lwtt(pga(1,1,1,jk), pgb(1,1,1,jk), zuu, ztt)

          DO jl = 1, klon
             zwtr = pdbsl(jl, 1, ibs)*ztt(jl, 1)*ztt(jl, 10) + &
                  pdbsl(jl, 2, ibs)*ztt(jl, 2)*ztt(jl, 7)*ztt(jl, 11) + &
                  pdbsl(jl, 3, ibs)*ztt(jl, 4)*ztt(jl, 8)*ztt(jl, 12) + &
                  pdbsl(jl, 4, ibs)*ztt(jl, 5)*ztt(jl, 9)*ztt(jl, 13) + &
                  pdbsl(jl, 5, ibs)*ztt(jl, 3)*ztt(jl, 14) + &
                  pdbsl(jl, 6, ibs)*ztt(jl, 6)*ztt(jl, 15)
             zglayu(jl) = zglayu(jl) + zwtr*wg1(jg)
          END DO

       END DO

       DO jl = 1, klon
          padjd(jl, jk) = zglayd(jl)
          pcntrb(jl, jk, jk+1) = zglayd(jl)
          padju(jl, jk+1) = zglayu(jl)
          pcntrb(jl, jk+1, jk) = zglayu(jl)
          pcntrb(jl, jk, jk) = 0.0
       END DO
    END DO

    DO jk = 1, llm
       jk2 = 2*jk
       jk1 = jk2 - 1
       DO jnu = 1, ninter
          DO jl = 1, klon
             pdbdt(jl, jnu, jk) = pdbsl(jl, jnu, jk1) + pdbsl(jl, jnu, jk2)
          END DO
       END DO
    END DO

  END SUBROUTINE lwvn

end module lwvn_m
