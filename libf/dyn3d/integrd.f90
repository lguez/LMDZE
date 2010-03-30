SUBROUTINE integrd(nq,vcovm1,ucovm1,tetam1,psm1,massem1,dv,du,dteta,dq,dp, &
     vcov,ucov,teta,q,ps,masse,phis,finvmaold,leapf, dt)

  ! From dyn3d/integrd.F,v 1.1.1.1 2004/05/19 12:53:05
  !   Auteur:  P. Le Van                                                  
  !   objet:                                                              
  !   Incrementation des tendances dynamiques                             

  USE dimens_m, ONLY : iim, llm
  USE paramet_m, ONLY : iip1, iip2, ijp1llm, ip1jm, ip1jmp1, jjp1, llmp1
  USE comvert, ONLY : ap, bp
  USE comgeom, ONLY : aire, apoln, apols
  USE pression_m, ONLY : pression
  USE filtreg_m, ONLY : filtreg

  IMPLICIT NONE

  !   Arguments:                                                          

  INTEGER nq

  REAL vcov(ip1jm,llm), ucov(ip1jmp1,llm), teta(ip1jmp1,llm)
  REAL q(ip1jmp1,llm,nq)
  REAL ps(ip1jmp1), masse(ip1jmp1,llm), phis(ip1jmp1)

  REAL vcovm1(ip1jm,llm), ucovm1(ip1jmp1,llm)
  REAL tetam1(ip1jmp1,llm), psm1(ip1jmp1), massem1(ip1jmp1,llm)

  REAL dv(ip1jm,llm), du(ip1jmp1,llm)
  REAL dteta(ip1jmp1,llm), dp(ip1jmp1)
  REAL dq(ip1jmp1,llm,nq), finvmaold(ip1jmp1,llm)
  LOGICAL, INTENT (IN) :: leapf
  real, intent(in):: dt

  !   Local:                                                              

  REAL vscr(ip1jm), uscr(ip1jmp1), hscr(ip1jmp1), pscr(ip1jmp1)
  REAL massescr(ip1jmp1,llm), finvmasse(ip1jmp1,llm)
  REAL p(ip1jmp1,llmp1)
  REAL tpn, tps, tppn(iim), tpps(iim)
  REAL qpn, qps, qppn(iim), qpps(iim)
  REAL deltap(ip1jmp1,llm)

  INTEGER l, ij, iq

  REAL ssum

  !-----------------------------------------------------------------------

  DO l = 1, llm
     DO ij = 1, iip1
        ucov(ij,l) = 0.
        ucov(ij+ip1jm,l) = 0.
        uscr(ij) = 0.
        uscr(ij+ip1jm) = 0.
     END DO
  END DO


  !    ............    integration  de       ps         ..............    

  CALL scopy(ip1jmp1*llm,masse,1,massescr,1)

  DO ij = 1, ip1jmp1
     pscr(ij) = ps(ij)
     ps(ij) = psm1(ij) + dt*dp(ij)
  END DO

  DO ij = 1, ip1jmp1
     IF (ps(ij)<0.) THEN
        PRINT *, ' Au point ij = ', ij, ' , pression sol neg. ', ps(ij)
        STOP 'integrd'
     END IF
  END DO

  DO ij = 1, iim
     tppn(ij) = aire(ij)*ps(ij)
     tpps(ij) = aire(ij+ip1jm)*ps(ij+ip1jm)
  END DO
  tpn = ssum(iim,tppn,1)/apoln
  tps = ssum(iim,tpps,1)/apols
  DO ij = 1, iip1
     ps(ij) = tpn
     ps(ij+ip1jm) = tps
  END DO

  !  ... Calcul  de la nouvelle masse d'air au dernier temps integre t+1 .

  CALL pression(ip1jmp1,ap,bp,ps,p)
  CALL massdair(p,masse)

  CALL scopy(ijp1llm,masse,1,finvmasse,1)
  CALL filtreg(finvmasse,jjp1,llm,-2,2,.TRUE.,1)


  !    ............   integration  de  ucov, vcov,  h     ..............  

  DO  l = 1, llm

     DO ij = iip2, ip1jm
        uscr(ij) = ucov(ij,l)
        ucov(ij,l) = ucovm1(ij,l) + dt*du(ij,l)
     END DO

     DO ij = 1, ip1jm
        vscr(ij) = vcov(ij,l)
        vcov(ij,l) = vcovm1(ij,l) + dt*dv(ij,l)
     END DO

     DO ij = 1, ip1jmp1
        hscr(ij) = teta(ij,l)
        teta(ij,l) = tetam1(ij,l)*massem1(ij,l)/masse(ij,l) + &
             dt*dteta(ij,l)/masse(ij,l)
     END DO

     !   ....  Calcul de la valeur moyenne, unique  aux poles pour  teta    .


     DO ij = 1, iim
        tppn(ij) = aire(ij)*teta(ij,l)
        tpps(ij) = aire(ij+ip1jm)*teta(ij+ip1jm,l)
     END DO
     tpn = ssum(iim,tppn,1)/apoln
     tps = ssum(iim,tpps,1)/apols

     DO ij = 1, iip1
        teta(ij,l) = tpn
        teta(ij+ip1jm,l) = tps
     END DO


     IF (leapf) THEN
        CALL scopy(ip1jmp1,uscr(1),1,ucovm1(1,l),1)
        CALL scopy(ip1jm,vscr(1),1,vcovm1(1,l),1)
        CALL scopy(ip1jmp1,hscr(1),1,tetam1(1,l),1)
     END IF

  END DO

  DO l = 1, llm
     DO ij = 1, ip1jmp1
        deltap(ij,l) = p(ij,l) - p(ij,l+1)
     END DO
  END DO

  CALL qminimum(q,nq,deltap)

  !    .....  Calcul de la valeur moyenne, unique  aux poles pour  q .....


  DO iq = 1, nq
     DO l = 1, llm

        DO ij = 1, iim
           qppn(ij) = aire(ij)*q(ij,l,iq)
           qpps(ij) = aire(ij+ip1jm)*q(ij+ip1jm,l,iq)
        END DO
        qpn = ssum(iim,qppn,1)/apoln
        qps = ssum(iim,qpps,1)/apols

        DO ij = 1, iip1
           q(ij,l,iq) = qpn
           q(ij+ip1jm,l,iq) = qps
        END DO

     END DO
  END DO


  CALL scopy(ijp1llm,finvmasse,1,finvmaold,1)


  !     .....   FIN  de l'integration  de   q    .......                  

  IF (leapf) THEN
     CALL scopy(ip1jmp1,pscr,1,psm1,1)
     CALL scopy(ip1jmp1*llm,massescr,1,massem1,1)
  END IF

END SUBROUTINE integrd