module sw1s_m

  IMPLICIT NONE

contains

  SUBROUTINE sw1s(knu, palbd, palbp, pcg, pcld, pclear, pdsig, pomega, poz, &
       prmu, psec, ptau, pud, pfd, pfu)
    
    use dimensions, only: llm
    use dimphy, only: klon
    use swclr_m, only: swclr
    use swr_m, only: swr

    ! PURPOSE.
    ! THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
    ! SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    ! METHOD.
    ! 1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
    ! CONTINUUM SCATTERING
    ! 2. MULTIPLY BY OZONE TRANSMISSION FUNCTION

    ! REFERENCE.
    ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    ! AUTHOR.
    ! JEAN-JACQUES MORCRETTE  *ECMWF*

    ! MODIFICATIONS.
    ! ORIGINAL : 89-07-14
    ! 94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO

    ! * ARGUMENTS:

    INTEGER knu
    DOUBLE PRECISION, intent(in):: palbd(klon, 2)
    DOUBLE PRECISION, intent(in):: palbp(klon, 2)
    DOUBLE PRECISION, intent(in):: pcg(klon, 2, llm)
    DOUBLE PRECISION pcld(klon, llm)
    DOUBLE PRECISION pclear(klon)
    DOUBLE PRECISION pdsig(klon, llm)
    DOUBLE PRECISION, intent(in):: pomega(klon, 2, llm)
    DOUBLE PRECISION poz(klon, llm)
    DOUBLE PRECISION prmu(klon)
    DOUBLE PRECISION psec(klon)
    DOUBLE PRECISION, intent(in):: ptau(klon, 2, llm)
    DOUBLE PRECISION pud(klon, 5, llm+1)

    DOUBLE PRECISION pfd(klon, llm+1)
    DOUBLE PRECISION pfu(klon, llm+1)

    ! LOCAL VARIABLES:

    INTEGER iind(4)

    DOUBLE PRECISION zcgaz(klon, llm)
    DOUBLE PRECISION zdiff(klon)
    DOUBLE PRECISION zdirf(klon)
    DOUBLE PRECISION zpizaz(klon, llm)
    DOUBLE PRECISION zrayl(klon)
    DOUBLE PRECISION zray1(klon, llm+1)
    DOUBLE PRECISION zray2(klon, llm+1)
    DOUBLE PRECISION zrefz(klon, 2, llm+1)
    DOUBLE PRECISION zrj(klon, 6, llm+1)
    DOUBLE PRECISION zrj0(klon, 6, llm+1)
    DOUBLE PRECISION zrk(klon, 6, llm+1)
    DOUBLE PRECISION zrk0(klon, 6, llm+1)
    DOUBLE PRECISION zrmue(klon, llm+1)
    DOUBLE PRECISION zrmu0(klon, llm+1)
    DOUBLE PRECISION zr(klon, 4)
    DOUBLE PRECISION ztauaz(klon, llm)
    DOUBLE PRECISION ztra1(klon, llm+1)
    DOUBLE PRECISION ztra2(klon, llm+1)
    DOUBLE PRECISION zw(klon, 4)

    INTEGER jl, jk, k, jaj, ikm1, ikl

    ! Prescribed Data:

    DOUBLE PRECISION rsun(2)
    SAVE rsun
    DOUBLE PRECISION rray(2, 6)
    SAVE rray
    DATA rsun(1)/0.441676d0/
    DATA rsun(2)/0.558324d0/
    DATA (rray(1,k), k=1, 6)/.428937d-01, .890743d+00, -.288555d+01, &
         .522744d+01, -.469173d+01, .161645d+01/
    DATA (rray(2,k), k=1, 6)/.697200d-02, .173297d-01, -.850903d-01, &
         .248261d+00, -.302031d+00, .129662d+00/
    ! ------------------------------------------------------------------

    ! *         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
    ! ----------------------- ------------------



    ! *         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
    ! -----------------------------------------


    DO jl = 1, klon
       zrayl(jl) = rray(knu, 1) + prmu(jl)*(rray(knu,2)+prmu(jl)*(rray(knu, &
            3)+prmu(jl)*(rray(knu,4)+prmu(jl)*(rray(knu,5)+prmu(jl)*rray(knu,6)))))
    END DO


    ! ------------------------------------------------------------------

    ! *         2.    CONTINUUM SCATTERING CALCULATIONS
    ! ---------------------------------


    ! *         2.1   CLEAR-SKY FRACTION OF THE COLUMN
    ! --------------------------------


    CALL swclr(knu, palbp, pdsig, zrayl, psec, zpizaz, zray1, zray2, zrefz, &
         zrj0, zrk0, zrmu0, ztauaz, ztra1, ztra2)


    ! *         2.2   CLOUDY FRACTION OF THE COLUMN
    ! -----------------------------

    zcgaz = 0d0
    CALL swr(knu, palbd, pcg, pcld, pomega, psec, ptau, zcgaz, &
         zpizaz, zray1, zray2, zrefz, zrj, zrk, zrmue, ztauaz, ztra1, ztra2)


    ! ------------------------------------------------------------------

    ! *         3.    OZONE ABSORPTION
    ! ----------------


    iind(1) = 1
    iind(2) = 3
    iind(3) = 1
    iind(4) = 3


    ! *         3.1   DOWNWARD FLUXES
    ! ---------------


    jaj = 2

    DO jl = 1, klon
       zw(jl, 1) = 0.
       zw(jl, 2) = 0.
       zw(jl, 3) = 0.
       zw(jl, 4) = 0.
       pfd(jl, llm+1) = ((1.-pclear(jl))*zrj(jl,jaj,llm+1)+pclear(jl)*zrj0( &
            jl,jaj,llm+1))*rsun(knu)
    END DO
    DO jk = 1, llm
       ikl = llm + 1 - jk
       DO jl = 1, klon
          zw(jl, 1) = zw(jl, 1) + pud(jl, 1, ikl)/zrmue(jl, ikl)
          zw(jl, 2) = zw(jl, 2) + poz(jl, ikl)/zrmue(jl, ikl)
          zw(jl, 3) = zw(jl, 3) + pud(jl, 1, ikl)/zrmu0(jl, ikl)
          zw(jl, 4) = zw(jl, 4) + poz(jl, ikl)/zrmu0(jl, ikl)
       END DO

       CALL swtt1(knu, 4, iind, zw, zr)

       DO jl = 1, klon
          zdiff(jl) = zr(jl, 1)*zr(jl, 2)*zrj(jl, jaj, ikl)
          zdirf(jl) = zr(jl, 3)*zr(jl, 4)*zrj0(jl, jaj, ikl)
          pfd(jl, ikl) = ((1.-pclear(jl))*zdiff(jl)+pclear(jl)*zdirf(jl))* &
               rsun(knu)
       END DO
    END DO


    ! *         3.2   UPWARD FLUXES
    ! -------------


    DO jl = 1, klon
       pfu(jl, 1) = ((1.-pclear(jl))*zdiff(jl)*palbd(jl,knu)+pclear(jl)*zdirf(jl &
            )*palbp(jl,knu))*rsun(knu)
    END DO

    DO jk = 2, llm + 1
       ikm1 = jk - 1
       DO jl = 1, klon
          zw(jl, 1) = zw(jl, 1) + pud(jl, 1, ikm1)*1.66
          zw(jl, 2) = zw(jl, 2) + poz(jl, ikm1)*1.66
          zw(jl, 3) = zw(jl, 3) + pud(jl, 1, ikm1)*1.66
          zw(jl, 4) = zw(jl, 4) + poz(jl, ikm1)*1.66
       END DO

       CALL swtt1(knu, 4, iind, zw, zr)

       DO jl = 1, klon
          zdiff(jl) = zr(jl, 1)*zr(jl, 2)*zrk(jl, jaj, jk)
          zdirf(jl) = zr(jl, 3)*zr(jl, 4)*zrk0(jl, jaj, jk)
          pfu(jl, jk) = ((1.-pclear(jl))*zdiff(jl)+pclear(jl)*zdirf(jl))* &
               rsun(knu)
       END DO
    END DO

  END SUBROUTINE sw1s

end module sw1s_m
