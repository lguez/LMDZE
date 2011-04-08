SUBROUTINE advect(ucov, vcov, teta, w, massebx, masseby, du, dv, dteta, conser)

  ! From dyn3d/advect.F, v 1.1.1.1 2004/05/19 12:53:06
  !   Auteurs:  P. Le Van , Fr. Hourdin  .                                
  !   Objet:                                                              
  !   .... calcul des termes d'advection vertic.pour u, v, teta, q ...       
  !        ces termes sont ajoutes a du, dv, dteta et dq .                  
  !  Modif F.Forget 03/94 : on retire q de advect                         

  USE dimens_m
  USE paramet_m
  USE comconst
  USE comvert
  USE comgeom
  USE ener

  IMPLICIT NONE

  !   Arguments:                                                          

  REAL vcov(ip1jm, llm), ucov(ip1jmp1, llm)
  real, intent(in):: teta(ip1jmp1, llm)
  REAL massebx(ip1jmp1, llm), masseby(ip1jm, llm)
  real, INTENT (IN):: w(ip1jmp1, llm)
  REAL dv(ip1jm, llm), du(ip1jmp1, llm), dteta(ip1jmp1, llm)
  LOGICAL, INTENT (IN):: conser

  !   Local:                                                              

  REAL uav(ip1jmp1, llm), vav(ip1jm, llm), wsur2(ip1jmp1)
  REAL unsaire2(ip1jmp1), ge(ip1jmp1)
  REAL deuxjour, ww, gt, uu, vv

  INTEGER ij, l

  REAL ssum

  !-----------------------------------------------------------------------
  !   2. Calculs preliminaires:                                           
  !   -------------------------                                           

  IF (conser) THEN
     deuxjour = 2.*daysec

     DO ij = 1, ip1jmp1
        unsaire2(ij) = unsaire(ij)*unsaire(ij)
     end DO
  END IF


  !------------------  -yy ---------------------------------------------- 
  !   .  Calcul de     u                                                  

  DO l = 1, llm
     DO ij = iip2, ip1jmp1
        uav(ij, l) = 0.25*(ucov(ij, l)+ucov(ij-iip1, l))
     END DO
     DO ij = iip2, ip1jm
        uav(ij, l) = uav(ij, l) + uav(ij+iip1, l)
     END DO
     DO ij = 1, iip1
        uav(ij, l) = 0.
        uav(ip1jm+ij, l) = 0.
     END DO
  END DO

  !------------------  -xx ---------------------------------------------- 
  !   .  Calcul de     v                                                  

  DO l = 1, llm
     DO ij = 2, ip1jm
        vav(ij, l) = 0.25*(vcov(ij, l)+vcov(ij-1, l))
     END DO
     DO ij = 1, ip1jm, iip1
        vav(ij, l) = vav(ij+iim, l)
     END DO
     DO ij = 1, ip1jm - 1
        vav(ij, l) = vav(ij, l) + vav(ij+1, l)
     END DO
     DO ij = 1, ip1jm, iip1
        vav(ij+iim, l) = vav(ij, l)
     END DO
  END DO

  !-----------------------------------------------------------------------


  DO l = 1, llmm1


     !       ......   calcul de  - w/2.    au niveau  l+1   .......          

     DO ij = 1, ip1jmp1
        wsur2(ij) = -0.5*w(ij, l+1)
     END DO


     !     .....................     calcul pour  du     ..................  

     DO ij = iip2, ip1jm - 1
        ww = wsur2(ij) + wsur2(ij+1)
        uu = 0.5*(ucov(ij, l)+ucov(ij, l+1))
        du(ij, l) = du(ij, l) - ww*(uu-uav(ij, l))/massebx(ij, l)
        du(ij, l+1) = du(ij, l+1) + ww*(uu-uav(ij, l+1))/massebx(ij, l+1)
     END DO

     !     .....  correction pour  du(iip1, j, l)  ........                    
     !     .....     du(iip1, j, l)= du(1, j, l)   .....                         

     !DIR$ IVDEP                                                             
     DO ij = iip1 + iip1, ip1jm, iip1
        du(ij, l) = du(ij-iim, l)
        du(ij, l+1) = du(ij-iim, l+1)
     END DO

     !     .................    calcul pour   dv      .....................  

     DO ij = 1, ip1jm
        ww = wsur2(ij+iip1) + wsur2(ij)
        vv = 0.5*(vcov(ij, l)+vcov(ij, l+1))
        dv(ij, l) = dv(ij, l) - ww*(vv-vav(ij, l))/masseby(ij, l)
        dv(ij, l+1) = dv(ij, l+1) + ww*(vv-vav(ij, l+1))/masseby(ij, l+1)
     END DO



     !     ............................................................      
     !     ...............    calcul pour   dh      ...................      
     !     ............................................................      

     !                       ---z                                            
     !       calcul de  - d( teta  * w )      qu'on ajoute a   dh            
     !                   ...............                                     

     DO ij = 1, ip1jmp1
        ww = wsur2(ij)*(teta(ij, l)+teta(ij, l+1))
        dteta(ij, l) = dteta(ij, l) - ww
        dteta(ij, l+1) = dteta(ij, l+1) + ww
     end DO

     IF (conser) THEN
        DO ij = 1, ip1jmp1
           ge(ij) = wsur2(ij)*wsur2(ij)*unsaire2(ij)
        end DO
        gt = ssum(ip1jmp1, ge, 1)
        gtot(l) = deuxjour*sqrt(gt/ip1jmp1)
     END IF

  END DO

END SUBROUTINE advect
