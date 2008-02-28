      SUBROUTINE caldyn0
     $ (itau,ucov,vcov,teta,ps,masse,pk,phis ,
     $  phi,w,pbaru,pbarv,time )

!     From dyn3d/caldyn0.F,v 1.1.1.1 2004/05/19 12:53:07

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

      INTEGER, intent(in):: itau
      REAL, intent(in):: vcov(ip1jm,llm),ucov(ip1jmp1,llm)
      real teta(ip1jmp1,llm)
      REAL, intent(in):: ps(ip1jmp1)
      real, intent(in):: phis(ip1jmp1)
      REAL, intent(in):: pk(iip1,jjp1,llm)
      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
      REAL phi(ip1jmp1,llm),masse(ip1jmp1,llm)
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL time

c   Local:
c   ------

      REAL ang(ip1jmp1,llm),p(ip1jmp1,llmp1)
      REAL massebx(ip1jmp1,llm),masseby(ip1jm,llm),psexbarxy(ip1jm)
      REAL vorpot(ip1jm,llm)
      REAL w(ip1jmp1,llm),ecin(ip1jmp1,llm),convm(ip1jmp1,llm)
      REAL bern(ip1jmp1,llm)
      REAL massebxy(ip1jm,llm), dp(ip1jmp1)
    

      INTEGER   ij,l

c-----------------------------------------------------------------------
      print *, "Call sequence information: caldyn0"
c   Calcul des tendances dynamiques:
c   --------------------------------

      CALL covcont  ( llm    , ucov    , vcov , ucont, vcont        )
      CALL pression ( ip1jmp1, ap      , bp   ,  ps  , p            )
      CALL psextbar (   ps   , psexbarxy                            )
      CALL massdair (    p   , masse                                )
      CALL massbar  (   masse, massebx , masseby                    )
      CALL massbarxy(   masse, massebxy                             )
      CALL flumass  ( massebx, masseby , vcont, ucont ,pbaru, pbarv )
      CALL convmas  (   pbaru, pbarv   , convm                      )

      DO ij =1, ip1jmp1
         dp( ij ) = convm( ij,1 ) / airesurg( ij )
      ENDDO

      CALL vitvert ( convm  , w                                  )
      CALL tourpot ( vcov   , ucov  , massebxy  , vorpot         )
      CALL enercin ( vcov   , ucov  , vcont     , ucont  , ecin  )
      CALL bernoui ( ip1jmp1, llm   , phi       , ecin   , bern  )

      DO l=1,llm
         DO ij=1,ip1jmp1
            ang(ij,l) = ucov(ij,l) + constang(ij)
         ENDDO
      ENDDO

        CALL sortvarc0
     $ ( itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time,vcov )

      RETURN
      END
