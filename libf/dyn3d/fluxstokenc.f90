SUBROUTINE fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,time_step,itau)

  !     Auteur :  F. Hourdin

  USE ioipsl
  USE dimens_m
  USE paramet_m
  USE comconst
  USE comvert
  USE comgeom
  USE temps
  USE tracstoke

  IMPLICIT NONE

  REAL, intent(in):: time_step
  real t_wrt, t_ops
  REAL pbaru(ip1jmp1,llm), pbarv(ip1jm,llm)
  REAL masse(ip1jmp1,llm), teta(ip1jmp1,llm), phi(ip1jmp1,llm)
  REAL phis(ip1jmp1)

  REAL pbaruc(ip1jmp1,llm), pbarvc(ip1jm,llm)
  REAL massem(ip1jmp1,llm), tetac(ip1jmp1,llm), phic(ip1jmp1,llm)

  REAL pbarug(ip1jmp1,llm), pbarvg(iip1,jjm,llm), wg(ip1jmp1,llm)

  REAL pbarvst(iip1,jjp1,llm), zistdyn
  REAL dtcum

  INTEGER iadvtr, ndex(1)
  INTEGER nscal
  REAL tst(1), ist(1), istp(1)
  INTEGER ij, l, irec, i, j
  INTEGER, INTENT (IN) :: itau
  INTEGER fluxid, fluxvid, fluxdid

  SAVE iadvtr, massem, pbaruc, pbarvc, irec
  SAVE phic, tetac
  LOGICAL first
  SAVE first
  DATA first/ .TRUE./
  DATA iadvtr/0/

  !-------------------------------------------------------------

  IF (first) THEN
     CALL initfluxsto(time_step,istdyn*time_step,istdyn*time_step,nqmx, &
          fluxid,fluxvid,fluxdid)

     ndex(1) = 0
     CALL histwrite(fluxid,'phis',1,phis)
     CALL histwrite(fluxid,'aire',1,aire)

     ndex(1) = 0
     nscal = 1
     tst(1) = time_step
     CALL histwrite(fluxdid,'dtvr',1,tst)
     ist(1) = istdyn
     CALL histwrite(fluxdid,'istdyn',1,ist)
     istp(1) = istphy
     CALL histwrite(fluxdid,'istphy',1,istp)

     first = .FALSE.
  END IF


  IF (iadvtr==0) THEN
     CALL initial0(ijp1llm,phic)
     CALL initial0(ijp1llm,tetac)
     CALL initial0(ijp1llm,pbaruc)
     CALL initial0(ijmllm,pbarvc)
  END IF

  !   accumulation des flux de masse horizontaux
  DO l = 1, llm
     DO ij = 1, ip1jmp1
        pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
        tetac(ij,l) = tetac(ij,l) + teta(ij,l)
        phic(ij,l) = phic(ij,l) + phi(ij,l)
     END DO
     DO ij = 1, ip1jm
        pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
     END DO
  END DO

  !   selection de la masse instantannee des mailles avant le transport.
  IF (iadvtr==0) THEN
     CALL scopy(ip1jmp1*llm,masse,1,massem,1)
  END IF

  iadvtr = iadvtr + 1


  !   Test pour savoir si on advecte a ce pas de temps
  IF (iadvtr==istdyn) THEN
     !    normalisation
     DO l = 1, llm
        DO ij = 1, ip1jmp1
           pbaruc(ij,l) = pbaruc(ij,l)/float(istdyn)
           tetac(ij,l) = tetac(ij,l)/float(istdyn)
           phic(ij,l) = phic(ij,l)/float(istdyn)
        END DO
        DO ij = 1, ip1jm
           pbarvc(ij,l) = pbarvc(ij,l)/float(istdyn)
        END DO
     END DO

     !   traitement des flux de masse avant advection.
     !     1. calcul de w
     !     2. groupement des mailles pres du pole.

     CALL groupe(massem,pbaruc,pbarvc,pbarug,pbarvg,wg)

     DO l = 1, llm
        DO j = 1, jjm
           DO i = 1, iip1
              pbarvst(i,j,l) = pbarvg(i,j,l)
           END DO
        END DO
        DO i = 1, iip1
           pbarvst(i,jjp1,l) = 0.
        END DO
     END DO

     iadvtr = 0
     PRINT *, 'ITAU auqel on stoke les fluxmasses', itau

     CALL histwrite(fluxid,'masse',itau,massem)
     CALL histwrite(fluxid,'pbaru',itau,pbarug)
     CALL histwrite(fluxvid,'pbarv',itau,pbarvg)
     CALL histwrite(fluxid,'w',itau,wg)
     CALL histwrite(fluxid,'teta',itau,tetac)
     CALL histwrite(fluxid,'phi',itau,phic)
  END IF ! if iadvtr.EQ.istdyn                                       

END SUBROUTINE fluxstokenc
