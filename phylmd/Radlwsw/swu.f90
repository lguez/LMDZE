module swu_m

  IMPLICIT NONE

contains

  SUBROUTINE swu(psct, pcldsw, ppmb, ppsol, prmu0, fract, ptave, pwv, paki, &
       pcld, pclear, pdsig, pfact, prmu, psec, pud)

    USE clesphys, only: rco2
    use dimphy, only: klon
    USE suphec_m, only: rg
    use dimensions, only: llm
    USE radepsi, only: zepscq, zepsec
    USE radopt, only: novlp

    ! ARGUMENTS:

    DOUBLE PRECISION, intent(in):: psct
    DOUBLE PRECISION, intent(in):: pcldsw(klon, llm)
    DOUBLE PRECISION, intent(in):: ppmb(klon, llm + 1)
    DOUBLE PRECISION, intent(in):: ppsol(klon)
    DOUBLE PRECISION, intent(in):: prmu0(klon)
    DOUBLE PRECISION, intent(in):: fract(klon)
    DOUBLE PRECISION, intent(in):: ptave(klon, llm)
    DOUBLE PRECISION, intent(in):: pwv(klon, llm)

    DOUBLE PRECISION paki(klon, 2)
    DOUBLE PRECISION pcld(klon, llm)
    DOUBLE PRECISION pclear(klon)
    DOUBLE PRECISION pdsig(klon, llm)
    DOUBLE PRECISION pfact(klon)
    DOUBLE PRECISION prmu(klon)
    DOUBLE PRECISION psec(klon)
    DOUBLE PRECISION pud(klon, 5, llm + 1)

    ! Local:

    INTEGER iind(2)
    DOUBLE PRECISION zc1j(klon, llm + 1)
    DOUBLE PRECISION zclear(klon)
    DOUBLE PRECISION zcloud(klon)
    DOUBLE PRECISION zn175(klon)
    DOUBLE PRECISION zn190(klon)
    DOUBLE PRECISION zo175(klon)
    DOUBLE PRECISION zo190(klon)
    DOUBLE PRECISION zsign(klon)
    DOUBLE PRECISION zr(klon, 2)
    DOUBLE PRECISION zsigo(klon)
    DOUBLE PRECISION zud(klon, 2)
    DOUBLE PRECISION zrth, zrtu, zwh2o, zdsco2, zdsh2o, zfppw
    INTEGER jl, jk, jkp1, jkl, ja

    ! Prescribed Data:

    DOUBLE PRECISION zpdh2o, zpdumg
    SAVE zpdh2o, zpdumg
    DOUBLE PRECISION zprh2o, zprumg
    SAVE zprh2o, zprumg
    DOUBLE PRECISION rtdh2o, rtdumg
    SAVE rtdh2o, rtdumg
    DOUBLE PRECISION rth2o, rtumg
    SAVE rth2o, rtumg
    DATA zpdh2o, zpdumg /0.8d0, 0.75d0/
    DATA zprh2o, zprumg /30000.d0, 30000.d0/
    DATA rtdh2o, rtdumg /0.40d0, 0.375d0/
    DATA rth2o, rtumg /240.d0, 240.d0/

    !------------------------------------------------------------------

    ! 1. COMPUTES AMOUNTS OF ABSORBERS

    iind(1) = 1
    iind(2) = 2

    ! 1.1 INITIALIZES QUANTITIES

    DO jl = 1, klon
       pud(jl, 1, llm + 1) = 0.
       pud(jl, 2, llm + 1) = 0.
       pud(jl, 3, llm + 1) = 0.
       pud(jl, 4, llm + 1) = 0.
       pud(jl, 5, llm + 1) = 0.
       pfact(jl) = prmu0(jl) * fract(jl) * psct
       prmu(jl) = sqrt(1224. * prmu0(jl) * prmu0(jl) + 1.) / 35.
       psec(jl) = 1. / prmu(jl)
       zc1j(jl, llm + 1) = 0.
    END DO

    ! 1.3 AMOUNTS OF ABSORBERS

    DO jl = 1, klon
       zud(jl, 1) = 0.
       zud(jl, 2) = 0.
       zo175(jl) = ppsol(jl)**(zpdumg + 1.)
       zo190(jl) = ppsol(jl)**(zpdh2o + 1.)
       zsigo(jl) = ppsol(jl)
       zclear(jl) = 1.
       zcloud(jl) = 0.
    END DO

    DO jk = 1, llm
       jkp1 = jk + 1
       jkl = llm + 1 - jk
       DO jl = 1, klon
          zrth = (rth2o / ptave(jl, jk))**rtdh2o
          zrtu = (rtumg / ptave(jl, jk))**rtdumg
          zwh2o = max(pwv(jl, jk), zepscq)
          zsign(jl) = 100. * ppmb(jl, jkp1)
          pdsig(jl, jk) = (zsigo(jl) - zsign(jl)) / ppsol(jl)
          zn175(jl) = zsign(jl)**(zpdumg + 1.)
          zn190(jl) = zsign(jl)**(zpdh2o + 1.)
          zdsco2 = zo175(jl) - zn175(jl)
          zdsh2o = zo190(jl) - zn190(jl)
          pud(jl, 1, jk) = 1. / (10. * rg * (zpdh2o + 1.)) / zprh2o**zpdh2o &
               * zdsh2o * zwh2o * zrth
          pud(jl, 2, jk) = 1. / (10. * rg * (zpdumg + 1.)) / zprumg**zpdumg &
               * zdsco2 * rco2 * zrtu
          zfppw = 1.6078 * zwh2o / (1. + 0.608 * zwh2o)
          pud(jl, 4, jk) = pud(jl, 1, jk) * zfppw
          pud(jl, 5, jk) = pud(jl, 1, jk) * (1. - zfppw)
          zud(jl, 1) = zud(jl, 1) + pud(jl, 1, jk)
          zud(jl, 2) = zud(jl, 2) + pud(jl, 2, jk)
          zsigo(jl) = zsign(jl)
          zo175(jl) = zn175(jl)
          zo190(jl) = zn190(jl)

          IF (novlp==1) THEN
             zclear(jl) = zclear(jl) &
                  * (1. - max(pcldsw(jl, jkl), zcloud(jl))) &
                  / (1. - min(zcloud(jl), 1. - zepsec))
             zc1j(jl, jkl) = 1.0 - zclear(jl)
             zcloud(jl) = pcldsw(jl, jkl)
          ELSE IF (novlp==2) THEN
             zcloud(jl) = max(pcldsw(jl, jkl), zcloud(jl))
             zc1j(jl, jkl) = zcloud(jl)
          ELSE IF (novlp==3) THEN
             zclear(jl) = zclear(jl) * (1. - pcldsw(jl, jkl))
             zcloud(jl) = 1.0 - zclear(jl)
             zc1j(jl, jkl) = zcloud(jl)
          END IF
       END DO
    END DO
    DO jl = 1, klon
       pclear(jl) = 1. - zc1j(jl, 1)
    END DO
    DO jk = 1, llm
       DO jl = 1, klon
          IF (pclear(jl)<1.) THEN
             pcld(jl, jk) = pcldsw(jl, jk) / (1. - pclear(jl))
          ELSE
             pcld(jl, jk) = 0.
          END IF
       END DO
    END DO

    ! 1.4 COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS

    DO ja = 1, 2
       DO jl = 1, klon
          zud(jl, ja) = zud(jl, ja) * psec(jl)
       END DO
    END DO

    CALL swtt1(2, 2, iind, zud, zr)

    DO ja = 1, 2
       DO jl = 1, klon
          paki(jl, ja) = - log(zr(jl, ja)) / zud(jl, ja)
       END DO
    END DO

  END SUBROUTINE swu

end module swu_m
