SUBROUTINE sortvarc0(ucov, teta, ps, masse, pk, phis, vorpot, phi, bern, dp)

  ! From dyn3d/sortvarc0.F, v 1.1.1.1 2004/05/19 12:53:07
  ! Auteur : P. Le Van
  ! Objet : sortie des variables de contrôle

  USE dimens_m, ONLY : iim, jjm, llm
  USE paramet_m, ONLY : iip1, iip2, ijp1llm, ip1jm, ip1jmp1, jjp1
  USE comconst, ONLY : daysec, g, omeg, rad
  USE comgeom, ONLY : aire, cu, rlatu
  USE ener, ONLY : ang0, etot0, ptot0, rmsdpdt, rmsv, stot0, ztot0
  use filtreg_m, only: filtreg

  IMPLICIT NONE

  !   Arguments:

  REAL, INTENT (IN) :: ucov(ip1jmp1, llm)
  REAL, INTENT(IN):: teta(ip1jmp1, llm)
  real masse(ip1jmp1, llm)
  REAL, INTENT (IN) :: ps(ip1jmp1)
  REAL, INTENT (IN) :: phis(ip1jmp1)
  REAL :: vorpot(ip1jm, llm)
  REAL :: phi(ip1jmp1, llm), bern(ip1jmp1, llm)
  REAL :: dp(ip1jmp1)
  REAL, INTENT (IN) :: pk(ip1jmp1, llm)

  !   Local:

  REAL :: vor(ip1jm), bernf(ip1jmp1, llm), ztotl(llm)
  REAL :: etotl(llm), stotl(llm), rmsvl(llm), angl(llm), ge(ip1jmp1)
  REAL :: cosphi(ip1jm), omegcosp(ip1jm)
  REAL radomeg
  REAL massebxy(ip1jm, llm)
  INTEGER l, ij

  REAL :: ssum

  !-----------------------------------------------------------------------

  PRINT *, 'Call sequence information: sortvarc0'

  CALL massbarxy(masse, massebxy)

  ! Calcul  de  rmsdpdt
  ge = dp*dp
  rmsdpdt = ssum(ip1jmp1, ge, 1) - ssum(jjp1, ge, iip1)
  rmsdpdt = daysec*1.E-2*sqrt(rmsdpdt / (iim * jjp1))
  CALL scopy(ijp1llm, bern, 1, bernf, 1)
  CALL filtreg(bernf, jjp1, llm, -2, 2, .TRUE.)

  !  Calcul du moment  angulaire

  radomeg = rad * omeg

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
        ge(ij) = masse(ij, l) &
             *(phis(ij)+teta(ij, l)*pk(ij, l)+bernf(ij, l)-phi(ij, l))
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
     angl(l) = rad / g &
          * (ssum(ip1jm-iip1, ge(iip2), 1)-ssum(jjm-1, ge(iip2), iip1))
  END DO

  DO ij = 1, ip1jmp1
     ge(ij) = ps(ij)*aire(ij)
  END DO
  ptot0 = ssum(ip1jmp1, ge, 1) - ssum(jjp1, ge, iip1)
  etot0 = ssum(llm, etotl, 1)
  ztot0 = ssum(llm, ztotl, 1)
  stot0 = ssum(llm, stotl, 1)
  rmsv = ssum(llm, rmsvl, 1)
  ang0 = ssum(llm, angl, 1)

  PRINT *, 'ptot0 = ', ptot0
  PRINT *, 'etot0 = ', etot0
  PRINT *, 'ztot0 = ', ztot0
  PRINT *, 'stot0 = ', stot0
  PRINT *, 'ang0 = ', ang0

END SUBROUTINE sortvarc0
