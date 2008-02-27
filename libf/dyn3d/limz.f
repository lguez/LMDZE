!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/limz.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE limz(s0,sz,sm,pente_max)
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
      real sz(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      integer n0,iadvplus(ip1jmp1,llm),nl(llm)
c
      REAL q(ip1jmp1,llm)
      real dzq(ip1jmp1,llm)


      REAL new_m,zm
      real dzqw(ip1jmp1)
      real adzqw(ip1jmp1),dzqmax(ip1jmp1)

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
               dzq(ij,l) = sz(ij,l) /sm(ij,l)
         ENDDO
       ENDDO

c   calcul de la pente en haut et en bas de la maille
       do ij=1,ip1jmp1
       do l = 1, llm-1
            dzqw(l)=q(ij,l+1)-q(ij,l)
         enddo
            dzqw(llm)=0.

         do  l=1,llm
            adzqw(l)=abs(dzqw(l))
         enddo

c   calcul de la pente maximum dans la maille en valeur absolue

         do l=2,llm-1
            dzqmax(l)=pente_max*min(adzqw(l-1),adzqw(l))
         enddo

c   calcul de la pente avec limitation

         do l=2,llm-1
            if(     dzqw(l-1)*dzqw(l).gt.0.
     &         .and. dzq(ij,l)*dzqw(l).gt.0.) then
              dzq(ij,l)=
     &         sign(min(abs(dzq(ij,l)),dzqmax(l)),dzq(ij,l))
            else
c   extremum local
               dzq(ij,l)=0.
            endif
         enddo

         DO  l=1,llm
               sz(ij,l) = dzq(ij,l)*sm(ij,l)
         ENDDO

       ENDDO

      RETURN
      END
