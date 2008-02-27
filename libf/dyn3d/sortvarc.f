!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/sortvarc.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE sortvarc
     $(itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time ,
     $ vcov )

c=======================================================================
c
c   Auteur:    P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   sortie des variables de controle
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      use temps
      use ener
      IMPLICIT NONE

c   Arguments:
c   ----------

      INTEGER, intent(in):: itau
      REAL ucov(ip1jmp1,llm),teta(ip1jmp1,llm),masse(ip1jmp1,llm)
      REAL vcov(ip1jm,llm)
      REAL ps(ip1jmp1),phis(ip1jmp1)
      REAL vorpot(ip1jm,llm)
      REAL phi(ip1jmp1,llm),bern(ip1jmp1,llm)
      REAL dp(ip1jmp1)
      REAL time
      REAL, intent(in):: pk(ip1jmp1,llm)

c   Local:
c   ------

      REAL vor(ip1jm),bernf(ip1jmp1,llm),ztotl(llm)
      REAL etotl(llm),stotl(llm),rmsvl(llm),angl(llm),ge(ip1jmp1)
      REAL cosphi(ip1jm),omegcosp(ip1jm)
      REAL dtvrs1j,rjour,heure,radsg,radomeg
      REAL rday, massebxy(ip1jm,llm)
      INTEGER  l, ij, imjmp1

      REAL       SSUM

c-----------------------------------------------------------------------

       dtvrs1j   = dtvr/daysec
       rjour     = FLOAT( INT( itau * dtvrs1j ))
       heure     = ( itau*dtvrs1j-rjour ) * 24.
       imjmp1    = iim * jjp1
       IF(ABS(heure - 24.).LE.0.0001 ) heure = 0.
c
       CALL massbarxy ( masse, massebxy )

c   .....  Calcul  de  rmsdpdt  .....

       ge(:)=dp(:)*dp(:)

       rmsdpdt = SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)
c
       rmsdpdt = daysec* 1.e-2 * SQRT(rmsdpdt/imjmp1) 

       CALL SCOPY( ijp1llm,bern,1,bernf,1 )
       CALL filtreg(bernf,jjp1,llm,-2,2,.TRUE.,1)

c   .....  Calcul du moment  angulaire   .....

       radsg    = rad /g
       radomeg  = rad * omeg
c
       DO ij=iip2,ip1jm
          cosphi( ij ) = COS(rlatu((ij-1)/iip1+1))
          omegcosp(ij) = radomeg   * cosphi(ij)
       ENDDO

c  ...  Calcul  de l'energie,de l'enstrophie,de l'entropie et de rmsv  .

       DO l=1,llm
          DO ij = 1,ip1jm
             vor(ij)=vorpot(ij,l)*vorpot(ij,l)*massebxy(ij,l)
          ENDDO
          ztotl(l)=(SSUM(ip1jm,vor,1)-SSUM(jjm,vor,iip1))

          DO ij = 1,ip1jmp1
             ge(ij)= masse(ij,l)*(phis(ij)+teta(ij,l)*pk(ij,l)  +
     s        bernf(ij,l)-phi(ij,l))
          ENDDO
          etotl(l) = SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)

          DO   ij   = 1, ip1jmp1
             ge(ij) = masse(ij,l)*teta(ij,l)
          ENDDO
          stotl(l)= SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)

          DO ij=1,ip1jmp1
             ge(ij)=masse(ij,l)*AMAX1(bernf(ij,l)-phi(ij,l),0.)
          ENDDO
          rmsvl(l)=2.*(SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1))

          DO ij =iip2,ip1jm
             ge(ij)=(ucov(ij,l)/cu(ij)+omegcosp(ij))*masse(ij,l) *
     *               cosphi(ij)
          ENDDO
          angl(l) = radsg *
     s    (SSUM(ip1jm-iip1,ge(iip2),1)-SSUM(jjm-1,ge(iip2),iip1))
      ENDDO

          DO ij=1,ip1jmp1
            ge(ij)= ps(ij)*aire(ij)
          ENDDO
      ptot  = SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1)
      etot  = SSUM(     llm, etotl, 1 )
      ztot  = SSUM(     llm, ztotl, 1 )
      stot  = SSUM(     llm, stotl, 1 )
      rmsv  = SSUM(     llm, rmsvl, 1 )
      ang   = SSUM(     llm,  angl, 1 )

      rday = FLOAT(INT ( day_ini + time ))
c
      IF(ptot0.eq.0.)  THEN
         PRINT 3500, itau, rday, heure,time
         PRINT*,'WARNING!!! On recalcule les valeurs initiales de :'
         PRINT*,'ptot,rmsdpdt,etot,ztot,stot,rmsv,ang'
         PRINT *, ptot,rmsdpdt,etot,ztot,stot,rmsv,ang
         etot0 = etot
         ptot0 = ptot
         ztot0 = ztot
         stot0 = stot
         ang0  = ang
      END IF

      etot= etot/etot0
      rmsv= SQRT(rmsv/ptot)
      ptot= ptot/ptot0
      ztot= ztot/ztot0
      stot= stot/stot0
      ang = ang /ang0


      PRINT 3500, itau, rday, heure, time
      PRINT 4000, ptot,rmsdpdt,etot,ztot,stot,rmsv,ang

      RETURN

3500   FORMAT(4x,'pas',i7,5x,'jour',f5.0,'heure',f5.1,4x 
     *   ,'date',f10.5)
4000   FORMAT(10x,'masse',4x,'rmsdpdt',7x,'energie',2x,'enstrophie'
     * ,2x,'entropie',3x,'rmsv',4x,'mt.ang',/,'GLOB  '
     .  ,f10.6,e13.6,5f10.3/
     * )
      END

