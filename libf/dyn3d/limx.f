!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/limx.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE limx(s0,sx,sm,pente_max)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real pente_max
      REAL s0(ip1jmp1,llm),sm(ip1jmp1,llm)
      real sx(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      integer n0,iadvplus(ip1jmp1,llm),nl(llm)
c
      REAL q(ip1jmp1,llm)
      real dxq(ip1jmp1,llm)


      REAL new_m,zm
      real dxqu(ip1jmp1)
      real adxqu(ip1jmp1),dxqmax(ip1jmp1)

      Logical extremum,first
      save first

      REAL      SSUM,CVMGP,CVMGT
      integer ismax,ismin
      EXTERNAL  SSUM, convflu,ismin,ismax
      EXTERNAL filtreg

      data first/.true./


       DO  l = 1,llm
         DO  ij=1,ip1jmp1
               q(ij,l) = s0(ij,l) / sm ( ij,l )
               dxq(ij,l) = sx(ij,l) /sm(ij,l)
         ENDDO
       ENDDO

c   calcul de la pente a droite et a gauche de la maille

      do l = 1, llm
         do ij=iip2,ip1jm-1
            dxqu(ij)=q(ij+1,l)-q(ij,l)
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            dxqu(ij)=dxqu(ij-iim)
         enddo

         do ij=iip2,ip1jm
            adxqu(ij)=abs(dxqu(ij))
         enddo

c   calcul de la pente maximum dans la maille en valeur absolue

         do ij=iip2+1,ip1jm
            dxqmax(ij)=pente_max*min(adxqu(ij-1),adxqu(ij))
         enddo

         do ij=iip1+iip1,ip1jm,iip1
            dxqmax(ij-iim)=dxqmax(ij)
         enddo

c   calcul de la pente avec limitation

         do ij=iip2+1,ip1jm
            if(     dxqu(ij-1)*dxqu(ij).gt.0.
     &         .and. dxq(ij,l)*dxqu(ij).gt.0.) then
              dxq(ij,l)=
     &         sign(min(abs(dxq(ij,l)),dxqmax(ij)),dxq(ij,l))
            else
c   extremum local
               dxq(ij,l)=0.
            endif
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            dxq(ij-iim,l)=dxq(ij,l)
         enddo

         DO  ij=1,ip1jmp1
               sx(ij,l) = dxq(ij,l)*sm(ij,l)
         ENDDO

       ENDDO

      RETURN
      END
