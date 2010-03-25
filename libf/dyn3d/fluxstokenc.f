      SUBROUTINE fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,
     . time_step,itau )

       USE IOIPSL
c
c     Auteur :  F. Hourdin
c
c
ccc   ..   Modif. P. Le Van  ( 20/12/97 )  ...
c
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      use temps
      use tracstoke
      IMPLICIT NONE
c

      REAL time_step,t_wrt, t_ops
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL masse(ip1jmp1,llm),teta(ip1jmp1,llm),phi(ip1jmp1,llm)
      REAL phis(ip1jmp1)

      REAL pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
      REAL massem(ip1jmp1,llm),tetac(ip1jmp1,llm),phic(ip1jmp1,llm)

      REAL pbarug(ip1jmp1,llm),pbarvg(iip1,jjm,llm),wg(ip1jmp1,llm)

      REAL pbarvst(iip1,jjp1,llm),zistdyn
	real dtcum

      INTEGER iadvtr,ndex(1) 
      integer nscal
      real tst(1),ist(1),istp(1)
      INTEGER ij,l,irec,i,j
      integer, intent(in):: itau
      INTEGER fluxid, fluxvid,fluxdid
 
      SAVE iadvtr, massem,pbaruc,pbarvc,irec
      SAVE phic,tetac
      logical first
      save first
      data first/.true./
      DATA iadvtr/0/

      if(first) then

	CALL initfluxsto(
     .  time_step,istdyn* time_step,istdyn* time_step,
     . nqmx, fluxid,fluxvid,fluxdid) 
	
	ndex(1) = 0
        call histwrite(fluxid, 'phis', 1, phis)
        call histwrite(fluxid, 'aire', 1, aire)
	
	ndex(1) = 0
        nscal = 1
        tst(1) = time_step
        call histwrite(fluxdid, 'dtvr', 1, tst)
        ist(1)=istdyn
        call histwrite(fluxdid, 'istdyn', 1, ist)
        istp(1)= istphy
        call histwrite(fluxdid, 'istphy', 1, istp)
	
	first = .false.

      endif


      IF(iadvtr.EQ.0) THEN
         CALL initial0(ijp1llm,phic)
         CALL initial0(ijp1llm,tetac)
         CALL initial0(ijp1llm,pbaruc)
         CALL initial0(ijmllm,pbarvc)
      ENDIF

c   accumulation des flux de masse horizontaux
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
            tetac(ij,l) = tetac(ij,l) + teta(ij,l)
            phic(ij,l) = phic(ij,l) + phi(ij,l)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
         ENDDO
      ENDDO

c   selection de la masse instantannee des mailles avant le transport.
      IF(iadvtr.EQ.0) THEN
         CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
      ENDIF

      iadvtr   = iadvtr+1


c   Test pour savoir si on advecte a ce pas de temps
      IF ( iadvtr.EQ.istdyn ) THEN
c    normalisation
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l)/float(istdyn)
            tetac(ij,l) = tetac(ij,l)/float(istdyn)
            phic(ij,l) = phic(ij,l)/float(istdyn)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l)/float(istdyn)
         ENDDO
      ENDDO

c   traitement des flux de masse avant advection.
c     1. calcul de w
c     2. groupement des mailles pres du pole.

        CALL groupe( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )

        do l=1,llm
           do j=1,jjm
              do i=1,iip1
                 pbarvst(i,j,l)=pbarvg(i,j,l)
              enddo
           enddo
           do i=1,iip1
              pbarvst(i,jjp1,l)=0.
           enddo
        enddo

         iadvtr=0
	Print*,'ITAU auqel on stoke les fluxmasses',itau
	
	call histwrite(fluxid, 'masse', itau, massem)
	
	call histwrite(fluxid, 'pbaru', itau, pbarug)
	
	call histwrite(fluxvid, 'pbarv', itau, pbarvg)
	
        call histwrite(fluxid, 'w' ,itau, wg) 
	
	call histwrite(fluxid, 'teta' ,itau, tetac) 
	
	call histwrite(fluxid, 'phi' ,itau, phic) 
	
C

      ENDIF ! if iadvtr.EQ.istdyn

      RETURN
      END
