SUBROUTINE dissip(vcov,ucov,teta,p,dv,du,dh)

  ! From dyn3d/dissip.F,v 1.1.1.1 2004/05/19 12:53:05
  ! Avec nouveaux operateurs star :  gradiv2 , divgrad2, nxgraro2  ...
  ! Auteur:  P. Le Van                                                  
  ! Objet: dissipation horizontale                                             

  USE dimens_m, ONLY : llm
  USE paramet_m, ONLY : iip1, iip2, ip1jm, ip1jmp1, llmp1
  USE comdissnew, ONLY : lstardis, nitergdiv, nitergrot, niterh
  USE inidissip_m, ONLY : dtdiss, tetah, tetaudiv, tetaurot

  IMPLICIT NONE

  !   Arguments:                                                          
  REAL :: vcov(ip1jm,llm), ucov(ip1jmp1,llm), teta(ip1jmp1,llm)
  REAL, INTENT (IN) :: p(ip1jmp1,llmp1)
  REAL :: dv(ip1jm,llm), du(ip1jmp1,llm), dh(ip1jmp1,llm)

  !   Local:                                                              
  REAL :: gdx(ip1jmp1,llm), gdy(ip1jm,llm)
  REAL :: grx(ip1jmp1,llm), gry(ip1jm,llm)
  REAL :: te1dt(llm), te2dt(llm), te3dt(llm)
  REAL :: deltapres(ip1jmp1,llm)

  INTEGER :: l, ij

  !-----------------------------------------------------------------------

  !   initialisations:                                                    

  DO l = 1, llm
     te1dt(l) = tetaudiv(l)*dtdiss
     te2dt(l) = tetaurot(l)*dtdiss
     te3dt(l) = tetah(l)*dtdiss
  END DO
  du = 0.
  dv = 0.
  dh = 0.

  !   Calcul de la dissipation:                                           

  !   Calcul de la partie   grad  ( div ) :                               

  IF (lstardis) THEN
     CALL gradiv2(llm,ucov,vcov,nitergdiv,gdx,gdy)
  ELSE
     CALL gradiv(llm,ucov,vcov,nitergdiv,gdx,gdy)
  END IF

  DO l = 1, llm

     DO ij = 1, iip1
        gdx(ij,l) = 0.
        gdx(ij+ip1jm,l) = 0.
     END DO

     DO ij = iip2, ip1jm
        du(ij,l) = du(ij,l) - te1dt(l)*gdx(ij,l)
     END DO
     DO ij = 1, ip1jm
        dv(ij,l) = dv(ij,l) - te1dt(l)*gdy(ij,l)
     END DO
  END DO

  !   calcul de la partie   n X grad ( rot ):                             

  IF (lstardis) THEN
     CALL nxgraro2(llm,ucov,vcov,nitergrot,grx,gry)
  ELSE
     CALL nxgrarot(llm,ucov,vcov,nitergrot,grx,gry)
  END IF


  DO l = 1, llm
     DO ij = 1, iip1
        grx(ij,l) = 0.
     END DO

     DO ij = iip2, ip1jm
        du(ij,l) = du(ij,l) - te2dt(l)*grx(ij,l)
     END DO
     DO ij = 1, ip1jm
        dv(ij,l) = dv(ij,l) - te2dt(l)*gry(ij,l)
     END DO
  END DO

  !   calcul de la partie   div ( grad ):                                 

  IF (lstardis) THEN

     DO l = 1, llm
        DO ij = 1, ip1jmp1
           deltapres(ij,l) = amax1(0.,p(ij,l)-p(ij,l+1))
        END DO
     END DO

     CALL divgrad2(llm,teta,deltapres,niterh,gdx)
  ELSE
     CALL divgrad(llm,teta,niterh,gdx)
  END IF

  DO l = 1, llm
     DO ij = 1, ip1jmp1
        dh(ij,l) = dh(ij,l) - te3dt(l)*gdx(ij,l)
     END DO
  END DO

END SUBROUTINE dissip
