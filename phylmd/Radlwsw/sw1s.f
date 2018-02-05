module sw1s_m

  IMPLICIT NONE

contains

  SUBROUTINE sw1s(knu, palbd, palbp, pcg, pcld, pclear, pdsig, pomega, poz, &
       prmu, psec, ptau, pud, pfd, pfu)
    
    USE dimens_m
    USE dimphy
    USE raddim
    use swclr_m, only: swclr
    use swr_m, only: swr

    ! ------------------------------------------------------------------
    ! PURPOSE.
    ! --------

    ! THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
    ! SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    ! METHOD.
    ! -------

    ! 1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
    ! CONTINUUM SCATTERING
    ! 2. MULTIPLY BY OZONE TRANSMISSION FUNCTION

    ! REFERENCE.
    ! ----------

    ! SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    ! DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    ! AUTHOR.
    ! -------
    ! JEAN-JACQUES MORCRETTE  *ECMWF*

    ! MODIFICATIONS.
    ! --------------
    ! ORIGINAL : 89-07-14
    ! 94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
    ! ------------------------------------------------------------------

    ! * ARGUMENTS:

    INTEGER knu
    DOUBLE PRECISION palbd(kdlon, 2)
    DOUBLE PRECISION palbp(kdlon, 2)
    DOUBLE PRECISION pcg(kdlon, 2, kflev)
    DOUBLE PRECISION pcld(kdlon, kflev)
    DOUBLE PRECISION pclear(kdlon)
    DOUBLE PRECISION pdsig(kdlon, kflev)
    DOUBLE PRECISION pomega(kdlon, 2, kflev)
    DOUBLE PRECISION poz(kdlon, kflev)
    DOUBLE PRECISION prmu(kdlon)
    DOUBLE PRECISION psec(kdlon)
    DOUBLE PRECISION ptau(kdlon, 2, kflev)
    DOUBLE PRECISION pud(kdlon, 5, kflev+1)

    DOUBLE PRECISION pfd(kdlon, kflev+1)
    DOUBLE PRECISION pfu(kdlon, kflev+1)

    ! * LOCAL VARIABLES:

    INTEGER iind(4)

    DOUBLE PRECISION zcgaz(kdlon, kflev)
    DOUBLE PRECISION zdiff(kdlon)
    DOUBLE PRECISION zdirf(kdlon)
    DOUBLE PRECISION zpizaz(kdlon, kflev)
    DOUBLE PRECISION zrayl(kdlon)
    DOUBLE PRECISION zray1(kdlon, kflev+1)
    DOUBLE PRECISION zray2(kdlon, kflev+1)
    DOUBLE PRECISION zrefz(kdlon, 2, kflev+1)
    DOUBLE PRECISION zrj(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrj0(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrk(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrk0(kdlon, 6, kflev+1)
    DOUBLE PRECISION zrmue(kdlon, kflev+1)
    DOUBLE PRECISION zrmu0(kdlon, kflev+1)
    DOUBLE PRECISION zr(kdlon, 4)
    DOUBLE PRECISION ztauaz(kdlon, kflev)
    DOUBLE PRECISION ztra1(kdlon, kflev+1)
    DOUBLE PRECISION ztra2(kdlon, kflev+1)
    DOUBLE PRECISION zw(kdlon, 4)

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


    DO jl = 1, kdlon
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

    DO jl = 1, kdlon
       zw(jl, 1) = 0.
       zw(jl, 2) = 0.
       zw(jl, 3) = 0.
       zw(jl, 4) = 0.
       pfd(jl, kflev+1) = ((1.-pclear(jl))*zrj(jl,jaj,kflev+1)+pclear(jl)*zrj0( &
            jl,jaj,kflev+1))*rsun(knu)
    END DO
    DO jk = 1, kflev
       ikl = kflev + 1 - jk
       DO jl = 1, kdlon
          zw(jl, 1) = zw(jl, 1) + pud(jl, 1, ikl)/zrmue(jl, ikl)
          zw(jl, 2) = zw(jl, 2) + poz(jl, ikl)/zrmue(jl, ikl)
          zw(jl, 3) = zw(jl, 3) + pud(jl, 1, ikl)/zrmu0(jl, ikl)
          zw(jl, 4) = zw(jl, 4) + poz(jl, ikl)/zrmu0(jl, ikl)
       END DO

       CALL swtt1(knu, 4, iind, zw, zr)

       DO jl = 1, kdlon
          zdiff(jl) = zr(jl, 1)*zr(jl, 2)*zrj(jl, jaj, ikl)
          zdirf(jl) = zr(jl, 3)*zr(jl, 4)*zrj0(jl, jaj, ikl)
          pfd(jl, ikl) = ((1.-pclear(jl))*zdiff(jl)+pclear(jl)*zdirf(jl))* &
               rsun(knu)
       END DO
    END DO


    ! *         3.2   UPWARD FLUXES
    ! -------------


    DO jl = 1, kdlon
       pfu(jl, 1) = ((1.-pclear(jl))*zdiff(jl)*palbd(jl,knu)+pclear(jl)*zdirf(jl &
            )*palbp(jl,knu))*rsun(knu)
    END DO

    DO jk = 2, kflev + 1
       ikm1 = jk - 1
       DO jl = 1, kdlon
          zw(jl, 1) = zw(jl, 1) + pud(jl, 1, ikm1)*1.66
          zw(jl, 2) = zw(jl, 2) + poz(jl, ikm1)*1.66
          zw(jl, 3) = zw(jl, 3) + pud(jl, 1, ikm1)*1.66
          zw(jl, 4) = zw(jl, 4) + poz(jl, ikm1)*1.66
       END DO

       CALL swtt1(knu, 4, iind, zw, zr)

       DO jl = 1, kdlon
          zdiff(jl) = zr(jl, 1)*zr(jl, 2)*zrk(jl, jaj, jk)
          zdirf(jl) = zr(jl, 3)*zr(jl, 4)*zrk0(jl, jaj, jk)
          pfu(jl, jk) = ((1.-pclear(jl))*zdiff(jl)+pclear(jl)*zdirf(jl))* &
               rsun(knu)
       END DO
    END DO

  END SUBROUTINE sw1s

end module sw1s_m
