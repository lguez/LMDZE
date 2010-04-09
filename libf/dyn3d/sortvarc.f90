module sortvarc_m

  IMPLICIT NONE

contains

  SUBROUTINE sortvarc(itau, ucov, teta, ps, masse, pk, phis, vorpot, phi, &
       bern, dp, time_0)

    ! From dyn3d/sortvarc.F, version 1.1.1.1 2004/05/19 12:53:07
    ! Author: P. Le Van
    ! Objet : sortie des variables de contrôle

    USE conf_gcm_m, ONLY: day_step
    USE dimens_m, ONLY : iim, jjm, llm
    USE paramet_m, ONLY : iip1, iip2, ijp1llm, ip1jm, ip1jmp1, jjp1
    USE comconst, ONLY : daysec, dtvr, g, omeg, rad
    USE comgeom, ONLY : aire, cu, rlatu
    USE dynetat0_m, ONLY : day_ini
    USE ener, ONLY : ang, ang0, etot, etot0, ptot, ptot0, rmsdpdt, rmsv, &
         stot, stot0, ztot, ztot0
    use filtreg_m, only: filtreg

    ! Arguments:
    INTEGER, INTENT (IN) :: itau
    REAL :: ucov(ip1jmp1, llm), teta(ip1jmp1, llm), masse(ip1jmp1, llm)
    REAL :: ps(ip1jmp1), phis(ip1jmp1)
    REAL :: vorpot(ip1jm, llm)
    REAL :: phi(ip1jmp1, llm), bern(ip1jmp1, llm)
    REAL :: dp(ip1jmp1)
    REAL, INTENT (IN):: time_0
    REAL, INTENT (IN):: pk(ip1jmp1, llm)

    ! Local:
    REAL :: vor(ip1jm), bernf(ip1jmp1, llm), ztotl(llm)
    REAL :: etotl(llm), stotl(llm), rmsvl(llm), angl(llm), ge(ip1jmp1)
    REAL :: cosphi(ip1jm), omegcosp(ip1jm)
    REAL :: dtvrs1j, rjour, heure, radsg, radomeg
    real massebxy(ip1jm, llm)
    INTEGER :: l, ij, imjmp1
    REAL :: ssum
    real time

    !-----------------------------------------------------------------------

    time = real(itau) / day_step + time_0
    dtvrs1j = dtvr/daysec
    rjour = real(int(itau*dtvrs1j))
    heure = (itau*dtvrs1j-rjour)*24.
    imjmp1 = iim*jjp1
    IF (abs(heure-24.)<=0.0001) heure = 0.

    CALL massbarxy(masse, massebxy)

    ! Calcul  de  rmsdpdt
    ge(:) = dp(:)*dp(:)
    rmsdpdt = ssum(ip1jmp1, ge, 1) - ssum(jjp1, ge, iip1)
    rmsdpdt = daysec*1.E-2*sqrt(rmsdpdt/imjmp1)
    CALL scopy(ijp1llm, bern, 1, bernf, 1)
    CALL filtreg(bernf, jjp1, llm, -2, 2, .TRUE., 1)

    ! Calcul du moment  angulaire
    radsg = rad/g
    radomeg = rad*omeg
    DO ij = iip2, ip1jm
       cosphi(ij) = cos(rlatu((ij-1)/iip1+1))
       omegcosp(ij) = radomeg*cosphi(ij)
    END DO

    ! Calcul  de l'energie, de l'enstrophie, de l'entropie et de rmsv
    DO l = 1, llm
       DO ij = 1, ip1jm
          vor(ij) = vorpot(ij, l)*vorpot(ij, l)*massebxy(ij, l)
       END DO
       ztotl(l) = (ssum(ip1jm, vor, 1)-ssum(jjm, vor, iip1))

       DO ij = 1, ip1jmp1
          ge(ij) = masse(ij, l) * (phis(ij) + teta(ij, l) * pk(ij, l) &
               + bernf(ij, l)-phi(ij, l))
       END DO
       etotl(l) = ssum(ip1jmp1, ge, 1) - ssum(jjp1, ge, iip1)

       DO ij = 1, ip1jmp1
          ge(ij) = masse(ij, l)*teta(ij, l)
       END DO
       stotl(l) = ssum(ip1jmp1, ge, 1) - ssum(jjp1, ge, iip1)

       DO ij = 1, ip1jmp1
          ge(ij) = masse(ij, l)*amax1(bernf(ij, l)-phi(ij, l), 0.)
       END DO
       rmsvl(l) = 2.*(ssum(ip1jmp1, ge, 1)-ssum(jjp1, ge, iip1))

       DO ij = iip2, ip1jm
          ge(ij) = (ucov(ij, l)/cu(ij)+omegcosp(ij))*masse(ij, l)*cosphi(ij)
       END DO
       angl(l) = radsg * (ssum(ip1jm-iip1, ge(iip2), 1) &
            - ssum(jjm-1, ge(iip2), iip1))
    END DO

    DO ij = 1, ip1jmp1
       ge(ij) = ps(ij)*aire(ij)
    END DO
    ptot = ssum(ip1jmp1, ge, 1) - ssum(jjp1, ge, iip1)
    etot = ssum(llm, etotl, 1)
    ztot = ssum(llm, ztotl, 1)
    stot = ssum(llm, stotl, 1)
    rmsv = ssum(llm, rmsvl, 1)
    ang = ssum(llm, angl, 1)

    IF (ptot0 == 0.) THEN
       PRINT *, 'WARNING!!! On recalcule les valeurs initiales de :'
       PRINT *, 'ptot, rmsdpdt, etot, ztot, stot, rmsv, ang'
       PRINT *, ptot, rmsdpdt, etot, ztot, stot, rmsv, ang
       etot0 = etot
       ptot0 = ptot
       ztot0 = ztot
       stot0 = stot
       ang0 = ang
    END IF

    etot = etot/etot0
    rmsv = sqrt(rmsv/ptot)
    ptot = ptot/ptot0
    ztot = ztot/ztot0
    stot = stot/stot0
    ang = ang/ang0

    PRINT 3500, itau, int(day_ini + time), heure, time
    PRINT 4000, ptot, rmsdpdt, etot, ztot, stot, rmsv, ang

3500 FORMAT (4X, 'pas', I7, 5X, 'jour', i5, 1X, 'heure', F5.1, 4X, 'date', &
          F10.5)
4000 FORMAT (10X, 'masse', 4X, 'rmsdpdt', 7X, 'energie', 2X, 'enstrophie', &
          2X, 'entropie', 3X, 'rmsv', 4X, 'mt.ang', /, 'GLOB  ', F10.6, &
          E13.6, 5F10.3/)

  END SUBROUTINE sortvarc

end module sortvarc_m
