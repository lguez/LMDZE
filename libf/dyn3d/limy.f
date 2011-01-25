!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/limy.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE limy(s0,sy,sm,pente_max)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      USE nr_util, ONLY : pi
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real pente_max
      real s0(ip1jmp1,llm),sy(ip1jmp1,llm),sm(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL q(ip1jmp1,llm)
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      real sigv,dyq(ip1jmp1),dyqv(ip1jm)
      real adyqv(ip1jm),dyqmax(ip1jmp1)
      REAL qbyv(ip1jm,llm)

      REAL qpns,qpsn,apn,aps,dyn1,dys1,dyn2,dys2
      Logical extremum,first
      save first

      real convpn,convps,convmpn,convmps
      real sinlon(iip1),sinlondlon(iip1)
      real coslon(iip1),coslondlon(iip1)
      save sinlon,coslon,sinlondlon,coslondlon
c
c
      REAL      SSUM
      integer ismax,ismin
      EXTERNAL  SSUM, convflu,ismin,ismax

      data first/.true./

      if(first) then
         print*,'SCHEMA AMONT NOUVEAU'
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
         enddo
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
      endif

c

      do l = 1, llm
c
         DO ij=1,ip1jmp1
               q(ij,l) = s0(ij,l) / sm ( ij,l )
               dyq(ij) = sy(ij,l) / sm ( ij,l )
         ENDDO
c
c   --------------------------------
c      CALCUL EN LATITUDE
c   --------------------------------

c   On commence par calculer la valeur du traceur moyenne sur le premier cercle
c   de latitude autour du pole (qpns pour le pole nord et qpsn pour
c    le pole nord) qui sera utilisee pour evaluer les pentes au pole.

      airej2 = SSUM( iim, aire(iip2), 1 )
      airejjm= SSUM( iim, aire(ip1jm -iim), 1 ) 
      DO i = 1, iim
      airescb(i) = aire(i+ iip1) * q(i+ iip1,l)
      airesch(i) = aire(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
      ENDDO
      qpns   = SSUM( iim,  airescb ,1 ) / airej2
      qpsn   = SSUM( iim,  airesch ,1 ) / airejjm

c   calcul des pentes aux points v

      do ij=1,ip1jm
         dyqv(ij)=q(ij,l)-q(ij+iip1,l)
         adyqv(ij)=abs(dyqv(ij))
      ENDDO

c   calcul des pentes aux points scalaires

      do ij=iip2,ip1jm
         dyqmax(ij)=min(adyqv(ij-iip1),adyqv(ij))
         dyqmax(ij)=pente_max*dyqmax(ij)
      enddo

c   calcul des pentes aux poles

c   calcul des pentes limites aux poles

c   cas ou on a un extremum au pole

c     if(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
c    &   apn=0.
c     if(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
c    &   dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
c    &   aps=0.

c   limitation des pentes aux poles
c     do ij=1,iip1
c        dyq(ij)=apn*dyq(ij)
c        dyq(ip1jm+ij)=aps*dyq(ip1jm+ij)
c     enddo

c   test
c      do ij=1,iip1
c         dyq(iip1+ij)=0.
c         dyq(ip1jm+ij-iip1)=0.
c      enddo
c      do ij=1,ip1jmp1
c         dyq(ij)=dyq(ij)*cos(rlatu((ij-1)/iip1+1))
c      enddo

      if(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
     &   then
         do ij=1,iip1
            dyqmax(ij)=0.
         enddo
      else
         do ij=1,iip1
            dyqmax(ij)=pente_max*abs(dyqv(ij))
         enddo
      endif

      if(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
     & dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
     &then
         do ij=ip1jm+1,ip1jmp1
            dyqmax(ij)=0.
         enddo
      else
         do ij=ip1jm+1,ip1jmp1
            dyqmax(ij)=pente_max*abs(dyqv(ij-iip1))
         enddo
      endif

c   calcul des pentes limitees

      do ij=1,ip1jmp1
         if(dyqv(ij)*dyqv(ij-iip1).gt.0.) then
            dyq(ij)=sign(min(abs(dyq(ij)),dyqmax(ij)),dyq(ij))
         else
            dyq(ij)=0.
         endif
      enddo

         DO ij=1,ip1jmp1
               sy(ij,l) = dyq(ij) * sm ( ij,l )
        ENDDO

      enddo ! fin de la boucle sur les couches verticales

      RETURN
      END
