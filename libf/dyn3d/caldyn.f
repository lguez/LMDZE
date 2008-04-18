!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/caldyn.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
c
c
      SUBROUTINE caldyn
     $ (itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
     $  phi,conser,du,dv,dteta,dp,w,pbaru,pbarv,time )

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      use pression_m, only: pression

      IMPLICIT NONE

c=======================================================================
c
c  Auteur :  P. Le Van
c
c   Objet:
c   ------
c
c   Calcul des tendances dynamiques.
c
c Modif 04/93 F.Forget
c=======================================================================

c-----------------------------------------------------------------------
c   0. Declarations:
c   ----------------


c   Arguments:
c   ----------

      LOGICAL conser

      INTEGER, intent(in):: itau
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
      REAL ps(ip1jmp1),phis(ip1jmp1)
      REAL, intent(in):: pk(iip1,jjp1,llm)
      real pkf(ip1jmp1,llm)
      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
      REAL phi(ip1jmp1,llm),masse(ip1jmp1,llm)
      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL dteta(ip1jmp1,llm),dp(ip1jmp1)
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL time

c   Local:
c   ------

      REAL ang(ip1jmp1,llm),p(ip1jmp1,llmp1)
      REAL massebx(ip1jmp1,llm),masseby(ip1jm,llm),psexbarxy(ip1jm)
      REAL vorpot(ip1jm,llm)
      REAL w(ip1jmp1,llm),ecin(ip1jmp1,llm),convm(ip1jmp1,llm)
      REAL bern(ip1jmp1,llm)
      REAL massebxy(ip1jm,llm)
    
      INTEGER   ij,l

c-----------------------------------------------------------------------
c   Calcul des tendances dynamiques:
c   --------------------------------

      CALL covcont  ( llm    , ucov    , vcov , ucont, vcont        )
      CALL pression ( ip1jmp1, ap      , bp   ,  ps  , p            )
      CALL psextbar (   ps   , psexbarxy                            )
      CALL massdair (    p   , masse                                )
      CALL massbar  (   masse, massebx , masseby                    )
      call massbarxy(   masse, massebxy                             )
      CALL flumass  ( massebx, masseby , vcont, ucont ,pbaru, pbarv )
      CALL dteta1   (   teta , pbaru   , pbarv, dteta               )
      CALL convmas  (   pbaru, pbarv   , convm                      )

      DO ij =1, ip1jmp1
         dp( ij ) = convm( ij,1 ) / airesurg( ij )
      ENDDO

      CALL vitvert ( convm  , w                                  )
      CALL tourpot ( vcov   , ucov  , massebxy  , vorpot         )
      CALL dudv1   ( vorpot , pbaru , pbarv     , du     , dv    )
      CALL enercin ( vcov   , ucov  , vcont     , ucont  , ecin  )
      CALL bernoui ( ip1jmp1, llm   , phi       , ecin   , bern  )
      CALL dudv2   ( teta   , pkf   , bern      , du     , dv    )


      DO l=1,llm
         DO ij=1,ip1jmp1
            ang(ij,l) = ucov(ij,l) + constang(ij)
      ENDDO
      ENDDO


      CALL advect( ang, vcov, teta, w, massebx, masseby, du, dv,dteta ) 

C  WARNING probleme de peridocite de dv sur les PC/linux. Pb d'arrondi 
C          probablement. Observe sur le code compile avec pgf90 3.0-1 

      DO l = 1, llm
         DO ij = 1, ip1jm, iip1
           IF( dv(ij,l).NE.dv(ij+iim,l) )  THEN
c         PRINT *,'!!!ATTENTION!!! probleme de periodicite sur vcov',  
c    ,   ' dans caldyn'
c         PRINT *,' l,  ij = ', l, ij, ij+iim,dv(ij+iim,l),dv(ij,l)
          dv(ij+iim,l) = dv(ij,l)
          endif
         enddo
      enddo
c-----------------------------------------------------------------------
c   Sorties eventuelles des variables de controle:
c   ----------------------------------------------

      IF( conser )  THEN
        CALL sortvarc
     $ ( itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time,vcov )

      ENDIF

      RETURN
      END
