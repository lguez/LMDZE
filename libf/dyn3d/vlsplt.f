!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/vlsplt.F,v 1.2 2005/02/24 12:16:57 fairhead Exp $
!
c
c

      SUBROUTINE vlsplt(q,pente_max,masse,w,pbaru,pbarv,pdt)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c   pente_max facteur de limitation des pentes: 2 en general
c                                               0 pour un schema amont
c   pbaru,pbarv,w flux de masse en u ,v ,w
c   pdt pas de temps
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      IMPLICIT NONE
c

c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
c      REAL masse(iip1,jjp1,llm),pente_max
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm)
c      REAL q(iip1,jjp1,llm)
      REAL w(ip1jmp1,llm),pdt
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
      INTEGER ijlqmin,iqmin,jqmin,lqmin
c
      REAL zm(ip1jmp1,llm),newmasse
      REAL mu(ip1jmp1,llm)
      REAL mv(ip1jm,llm)
      REAL mw(ip1jmp1,llm+1)
      REAL zq(ip1jmp1,llm),zz
      REAL dqx(ip1jmp1,llm),dqy(ip1jmp1,llm),dqz(ip1jmp1,llm)
      REAL second,temps0,temps1,temps2,temps3
      REAL ztemps1,ztemps2,ztemps3
      REAL zzpbar, zzw
      LOGICAL testcpu
      SAVE testcpu
      SAVE temps1,temps2,temps3
      INTEGER iminn,imaxx

      REAL qmin,qmax
      DATA qmin,qmax/0.,1.e33/
      DATA testcpu/.false./
      DATA temps1,temps2,temps3/0.,0.,0./


        zzpbar = 0.5 * pdt
        zzw    = pdt
      DO l=1,llm
        DO ij = iip2,ip1jm
            mu(ij,l)=pbaru(ij,l) * zzpbar
         ENDDO
         DO ij=1,ip1jm
            mv(ij,l)=pbarv(ij,l) * zzpbar
         ENDDO
         DO ij=1,ip1jmp1
            mw(ij,l)=w(ij,l) * zzw
         ENDDO
      ENDDO

      DO ij=1,ip1jmp1
         mw(ij,llm+1)=0.
      ENDDO
      
      CALL SCOPY(ijp1llm,q,1,zq,1)
      CALL SCOPY(ijp1llm,masse,1,zm,1)

cprint*,'Entree vlx1'
c	call minmaxq(zq,qmin,qmax,'avant vlx     ')
      call vlx(zq,pente_max,zm,mu)
cprint*,'Sortie vlx1'
c	call minmaxq(zq,qmin,qmax,'apres vlx1    ')

c print*,'Entree vly1'
      call vly(zq,pente_max,zm,mv)
c	call minmaxq(zq,qmin,qmax,'apres vly1     ')
cprint*,'Sortie vly1'
      call vlz(zq,pente_max,zm,mw)
c	call minmaxq(zq,qmin,qmax,'apres vlz     ')


      call vly(zq,pente_max,zm,mv)
c	call minmaxq(zq,qmin,qmax,'apres vly     ')


      call vlx(zq,pente_max,zm,mu)
c	call minmaxq(zq,qmin,qmax,'apres vlx2    ')
	

      DO l=1,llm
         DO ij=1,ip1jmp1
           q(ij,l)=zq(ij,l)
         ENDDO
         DO ij=1,ip1jm+1,iip1
            q(ij+iim,l)=q(ij,l)
         ENDDO
      ENDDO

      RETURN
      END
      SUBROUTINE vlx(q,pente_max,masse,u_m)

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
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL u_m( ip1jmp1,llm ),pbarv( iip1,jjm,llm)
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      INTEGER n0,iadvplus(ip1jmp1,llm),nl(llm)
c
      REAL new_m,zu_m,zdum(ip1jmp1,llm)
      REAL sigu(ip1jmp1),dxq(ip1jmp1,llm),dxqu(ip1jmp1)
      REAL zz(ip1jmp1)
      REAL adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
      REAL u_mq(ip1jmp1,llm)

      Logical extremum,first,testcpu
      SAVE first,testcpu

      REAL      SSUM
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5

      REAL z1,z2,z3

      DATA first,testcpu/.true.,.false./

      IF(first) THEN
         temps1=0.
         temps2=0.
         temps3=0.
         temps4=0.
         temps5=0.
         first=.false.
      ENDIF

c   calcul de la pente a droite et a gauche de la maille


      IF (pente_max.gt.-1.e-5) THEN
c       IF (pente_max.gt.10) THEN

c   calcul des pentes avec limitation, Van Leer scheme I:
c   -----------------------------------------------------

c   calcul de la pente aux points u
         DO l = 1, llm
            DO ij=iip2,ip1jm-1
               dxqu(ij)=q(ij+1,l)-q(ij,l)
c              IF(u_m(ij,l).lt.0.) stop'limx n admet pas les U<0'
c              sigu(ij)=u_m(ij,l)/masse(ij,l)
            ENDDO
            DO ij=iip1+iip1,ip1jm,iip1
               dxqu(ij)=dxqu(ij-iim)
c              sigu(ij)=sigu(ij-iim)
            ENDDO

            DO ij=iip2,ip1jm
               adxqu(ij)=abs(dxqu(ij))
            ENDDO

c   calcul de la pente maximum dans la maille en valeur absolue

            DO ij=iip2+1,ip1jm
               dxqmax(ij,l)=pente_max*
     ,      min(adxqu(ij-1),adxqu(ij))
c limitation subtile
c    ,      min(adxqu(ij-1)/sigu(ij-1),adxqu(ij)/(1.-sigu(ij)))
          

            ENDDO

            DO ij=iip1+iip1,ip1jm,iip1
               dxqmax(ij-iim,l)=dxqmax(ij,l)
            ENDDO

            DO ij=iip2+1,ip1jm
               IF(dxqu(ij-1)*dxqu(ij).gt.0) THEN
                  dxq(ij,l)=dxqu(ij-1)+dxqu(ij)
               ELSE
c   extremum local
                  dxq(ij,l)=0.
               ENDIF
               dxq(ij,l)=0.5*dxq(ij,l)
               dxq(ij,l)=
     ,         sign(min(abs(dxq(ij,l)),dxqmax(ij,l)),dxq(ij,l))
            ENDDO

         ENDDO ! l=1,llm
cprint*,'Ok calcul des pentes'

      ELSE ! (pente_max.lt.-1.e-5)

c   Pentes produits:
c   ----------------

         DO l = 1, llm
            DO ij=iip2,ip1jm-1
               dxqu(ij)=q(ij+1,l)-q(ij,l)
            ENDDO
            DO ij=iip1+iip1,ip1jm,iip1
               dxqu(ij)=dxqu(ij-iim)
            ENDDO

            DO ij=iip2+1,ip1jm
               zz(ij)=dxqu(ij-1)*dxqu(ij)
               zz(ij)=zz(ij)+zz(ij)
               IF(zz(ij).gt.0) THEN
                  dxq(ij,l)=zz(ij)/(dxqu(ij-1)+dxqu(ij))
               ELSE
c   extremum local
                  dxq(ij,l)=0.
               ENDIF
            ENDDO

         ENDDO

      ENDIF ! (pente_max.lt.-1.e-5)

c   bouclage de la pente en iip1:
c   -----------------------------

      DO l=1,llm
         DO ij=iip1+iip1,ip1jm,iip1
            dxq(ij-iim,l)=dxq(ij,l)
         ENDDO
         DO ij=1,ip1jmp1
            iadvplus(ij,l)=0
         ENDDO

      ENDDO

c print*,'Bouclage en iip1'

c   calcul des flux a gauche et a droite

c   on cumule le flux correspondant a toutes les mailles dont la masse
c   au travers de la paroi pENDant le pas de temps.
cprint*,'Cumule ....'

      DO l=1,llm
       DO ij=iip2,ip1jm-1
c	print*,'masse(',ij,')=',masse(ij,l)
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
             u_mq(ij,l)=u_m(ij,l)*(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
             u_mq(ij,l)=u_m(ij,l)*(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l))
          ENDIF
       ENDDO
      ENDDO
c	stop

c	go to 9999
c   detection des points ou on advecte plus que la masse de la
c   maille
      DO l=1,llm
         DO ij=iip2,ip1jm-1
            IF(zdum(ij,l).lt.0) THEN
               iadvplus(ij,l)=1
               u_mq(ij,l)=0.
            ENDIF
         ENDDO
      ENDDO
cprint*,'Ok test 1'
      DO l=1,llm
       DO ij=iip1+iip1,ip1jm,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO
c print*,'Ok test 2'


c   traitement special pour le cas ou on advecte en longitude plus que le
c   contenu de la maille.
c   cette partie est mal vectorisee.

c  calcul du nombre de maille sur lequel on advecte plus que la maille.

      n0=0
      DO l=1,llm
         nl(l)=0
         DO ij=iip2,ip1jm
            nl(l)=nl(l)+iadvplus(ij,l)
         ENDDO
         n0=n0+nl(l)
      ENDDO

      IF(n0.gt.0) THEN
CC      PRINT*,'Nombre de points pour lesquels on advect plus que le'
CC     &       ,'contenu de la maille : ',n0

         DO l=1,llm
            IF(nl(l).gt.0) THEN
               iju=0
c   indicage des mailles concernees par le traitement special
               DO ij=iip2,ip1jm
                  IF(iadvplus(ij,l).eq.1.and.mod(ij,iip1).ne.0) THEN
                     iju=iju+1
                     indu(iju)=ij
                  ENDIF
               ENDDO
               niju=iju
c              PRINT*,'niju,nl',niju,nl(l)

c  traitement des mailles
               DO iju=1,niju
                  ij=indu(iju)
                  j=(ij-1)/iip1+1
                  zu_m=u_m(ij,l)
                  u_mq(ij,l)=0.
                  IF(zu_m.gt.0.) THEN
                     ijq=ij
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)+q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m-masse(ijq,l)
                        i=mod(i-2+iim,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*
     &               (q(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq(ijq,l))
                  ELSE
                     ijq=ij+1
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(-zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)-q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m+masse(ijq,l)
                        i=mod(i,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l)-
     &               0.5*(1.+zu_m/masse(ijq,l))*dxq(ijq,l))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF  ! n0.gt.0 
9999    continue


c   bouclage en latitude
cprint*,'cvant bouclage en latitude'
      DO l=1,llm
        DO ij=iip1+iip1,ip1jm,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO


c   calcul des tENDances

      DO l=1,llm
         DO ij=iip2+1,ip1jm
            new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+
     &      u_mq(ij-1,l)-u_mq(ij,l))
     &      /new_m
            masse(ij,l)=new_m
         ENDDO
c   ModIF Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         DO ij=iip1+iip1,ip1jm,iip1
            q(ij-iim,l)=q(ij,l)
            masse(ij-iim,l)=masse(ij,l)
         ENDDO
      ENDDO
c     CALL SCOPY((jjm-1)*llm,q(iip1+iip1,1),iip1,q(iip2,1),iip1)
c     CALL SCOPY((jjm-1)*llm,masse(iip1+iip1,1),iip1,masse(iip2,1),iip1)


      RETURN
      END
      SUBROUTINE vly(q,pente_max,masse,masse_adv_v)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
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
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL masse_adv_v( ip1jm,llm)
      REAL q(ip1jmp1,llm), dq( ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      REAL dyq(ip1jmp1,llm),dyqv(ip1jm),zdvm(ip1jmp1,llm)
      REAL adyqv(ip1jm),dyqmax(ip1jmp1)
      REAL qbyv(ip1jm,llm)

      REAL qpns,qpsn,apn,aps,dyn1,dys1,dyn2,dys2,newmasse,fn,fs
c     REAL newq,oldmasse
      Logical extremum,first,testcpu
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
      SAVE first,testcpu

      REAL convpn,convps,convmpn,convmps
      real massepn,masseps,qpn,qps
      REAL sinlon(iip1),sinlondlon(iip1)
      REAL coslon(iip1),coslondlon(iip1)
      SAVE sinlon,coslon,sinlondlon,coslondlon
      SAVE airej2,airejjm
c
c
      REAL      SSUM

      DATA first,testcpu/.true.,.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./

      IF(first) THEN
         PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
         ENDDO
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
         airej2 = SSUM( iim, aire(iip2), 1 )
         airejjm= SSUM( iim, aire(ip1jm -iim), 1 ) 
      ENDIF

c
cPRINT*,'CALCUL EN LATITUDE'

      DO l = 1, llm
c
c   --------------------------------
c      CALCUL EN LATITUDE
c   --------------------------------

c   On commence par calculer la valeur du traceur moyenne sur le premier cercle
c   de latitude autour du pole (qpns pour le pole nord et qpsn pour
c    le pole nord) qui sera utilisee pour evaluer les pentes au pole.

      DO i = 1, iim
      airescb(i) = aire(i+ iip1) * q(i+ iip1,l)
      airesch(i) = aire(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
      ENDDO
      qpns   = SSUM( iim,  airescb ,1 ) / airej2
      qpsn   = SSUM( iim,  airesch ,1 ) / airejjm

c   calcul des pentes aux points v

      DO ij=1,ip1jm
         dyqv(ij)=q(ij,l)-q(ij+iip1,l)
         adyqv(ij)=abs(dyqv(ij))
      ENDDO

c   calcul des pentes aux points scalaires

      DO ij=iip2,ip1jm
         dyq(ij,l)=.5*(dyqv(ij-iip1)+dyqv(ij))
         dyqmax(ij)=min(adyqv(ij-iip1),adyqv(ij))
         dyqmax(ij)=pente_max*dyqmax(ij)
      ENDDO

c   calcul des pentes aux poles

      DO ij=1,iip1
         dyq(ij,l)=qpns-q(ij+iip1,l)
         dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l)-qpsn
      ENDDO

c   filtrage de la derivee
      dyn1=0.
      dys1=0.
      dyn2=0.
      dys2=0.
      DO ij=1,iim
         dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
         dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
         dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
         dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
      ENDDO
      DO ij=1,iip1
         dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
         dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
      ENDDO

c   calcul des pentes limites aux poles

      goto 8888
      fn=1.
      fs=1.
      DO ij=1,iim
         IF(pente_max*adyqv(ij).lt.abs(dyq(ij,l))) THEN
            fn=min(pente_max*adyqv(ij)/abs(dyq(ij,l)),fn)
         ENDIF
      IF(pente_max*adyqv(ij+ip1jm-iip1).lt.abs(dyq(ij+ip1jm,l))) THEN
         fs=min(pente_max*adyqv(ij+ip1jm-iip1)/abs(dyq(ij+ip1jm,l)),fs)
         ENDIF
      ENDDO
      DO ij=1,iip1
         dyq(ij,l)=fn*dyq(ij,l)
         dyq(ip1jm+ij,l)=fs*dyq(ip1jm+ij,l)
      ENDDO
8888    continue
      DO ij=1,iip1
         dyq(ij,l)=0.
         dyq(ip1jm+ij,l)=0.
      ENDDO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  En memoire de dIFferents tests sur la 
C  limitation des pentes aux poles.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT*,dyq(1)
C     PRINT*,dyqv(iip1+1)
C     apn=abs(dyq(1)/dyqv(iip1+1))
C     PRINT*,dyq(ip1jm+1)
C     PRINT*,dyqv(ip1jm-iip1+1)
C     aps=abs(dyq(ip1jm+1)/dyqv(ip1jm-iip1+1))
C     DO ij=2,iim
C        apn=amax1(abs(dyq(ij)/dyqv(ij)),apn)
C        aps=amax1(abs(dyq(ip1jm+ij)/dyqv(ip1jm-iip1+ij)),aps)
C     ENDDO
C     apn=min(pente_max/apn,1.)
C     aps=min(pente_max/aps,1.)
C
C
C   cas ou on a un extremum au pole
C
C     IF(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
C    &   apn=0.
C     IF(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
C    &   dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
C    &   aps=0.
C
C   limitation des pentes aux poles
C     DO ij=1,iip1
C        dyq(ij)=apn*dyq(ij)
C        dyq(ip1jm+ij)=aps*dyq(ip1jm+ij)
C     ENDDO
C
C   test
C      DO ij=1,iip1
C         dyq(iip1+ij)=0.
C         dyq(ip1jm+ij-iip1)=0.
C      ENDDO
C      DO ij=1,ip1jmp1
C         dyq(ij)=dyq(ij)*cos(rlatu((ij-1)/iip1+1))
C      ENDDO
C
C changement 10 07 96
C     IF(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
C    &   THEN
C        DO ij=1,iip1
C           dyqmax(ij)=0.
C        ENDDO
C     ELSE
C        DO ij=1,iip1
C           dyqmax(ij)=pente_max*abs(dyqv(ij))
C        ENDDO
C     ENDIF
C
C     IF(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
C    & dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
C    &THEN
C        DO ij=ip1jm+1,ip1jmp1
C           dyqmax(ij)=0.
C        ENDDO
C     ELSE
C        DO ij=ip1jm+1,ip1jmp1
C           dyqmax(ij)=pente_max*abs(dyqv(ij-iip1))
C        ENDDO
C     ENDIF
C   fin changement 10 07 96
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c   calcul des pentes limitees

      DO ij=iip2,ip1jm
         IF(dyqv(ij)*dyqv(ij-iip1).gt.0.) THEN
            dyq(ij,l)=sign(min(abs(dyq(ij,l)),dyqmax(ij)),dyq(ij,l))
         ELSE
            dyq(ij,l)=0.
         ENDIF
      ENDDO

      ENDDO

      DO l=1,llm
       DO ij=1,ip1jm
          IF(masse_adv_v(ij,l).gt.0) THEN
              qbyv(ij,l)=q(ij+iip1,l)+dyq(ij+iip1,l)*
     ,                   0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
          ELSE
              qbyv(ij,l)=q(ij,l)-dyq(ij,l)*
     ,                   0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
          ENDIF
          qbyv(ij,l)=masse_adv_v(ij,l)*qbyv(ij,l)
       ENDDO
      ENDDO


      DO l=1,llm
         DO ij=iip2,ip1jm
            newmasse=masse(ij,l)
     &      +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l))
     &         /newmasse
            masse(ij,l)=newmasse
         ENDDO
c.-. ancienne version
c        convpn=SSUM(iim,qbyv(1,l),1)/apoln
c        convmpn=ssum(iim,masse_adv_v(1,l),1)/apoln

         convpn=SSUM(iim,qbyv(1,l),1)
         convmpn=ssum(iim,masse_adv_v(1,l),1)
         massepn=ssum(iim,masse(1,l),1)
         qpn=0.
         do ij=1,iim
            qpn=qpn+masse(ij,l)*q(ij,l)
         enddo
         qpn=(qpn+convpn)/(massepn+convmpn)
         do ij=1,iip1
            q(ij,l)=qpn
         enddo

c        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)/apols
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)/apols

         convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
         convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
         masseps=ssum(iim, masse(ip1jm+1,l),1)
         qps=0.
         do ij = ip1jm+1,ip1jmp1-1
            qps=qps+masse(ij,l)*q(ij,l)
         enddo
         qps=(qps+convps)/(masseps+convmps)
         do ij=ip1jm+1,ip1jmp1
            q(ij,l)=qps
         enddo

c.-. fin ancienne version

c._. nouvelle version
c        convpn=SSUM(iim,qbyv(1,l),1)
c        convmpn=ssum(iim,masse_adv_v(1,l),1)
c        oldmasse=ssum(iim,masse(1,l),1)
c        newmasse=oldmasse+convmpn
c        newq=(q(1,l)*oldmasse+convpn)/newmasse
c        newmasse=newmasse/apoln
c        DO ij = 1,iip1
c           q(ij,l)=newq
c           masse(ij,l)=newmasse*aire(ij)
c        ENDDO
c        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
c        oldmasse=ssum(iim,masse(ip1jm-iim,l),1)
c        newmasse=oldmasse+convmps
c        newq=(q(ip1jmp1,l)*oldmasse+convps)/newmasse
c        newmasse=newmasse/apols
c        DO ij = ip1jm+1,ip1jmp1
c           q(ij,l)=newq
c           masse(ij,l)=newmasse*aire(ij)
c        ENDDO
c._. fin nouvelle version
      ENDDO

      RETURN
      END
      SUBROUTINE vlz(q,pente_max,masse,w)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c    q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm+1)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
c
      REAL wq(ip1jmp1,llm+1),newmasse

      REAL dzq(ip1jmp1,llm),dzqw(ip1jmp1,llm),adzqw(ip1jmp1,llm),dzqmax
      REAL sigw

      LOGICAL testcpu
      SAVE testcpu

      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
      REAL      SSUM

      DATA testcpu/.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./

c    On oriente tout dans le sens de la pression c'est a dire dans le
c    sens de W

      DO l=2,llm
         DO ij=1,ip1jmp1
            dzqw(ij,l)=q(ij,l-1)-q(ij,l)
            adzqw(ij,l)=abs(dzqw(ij,l))
         ENDDO
      ENDDO

      DO l=2,llm-1
         DO ij=1,ip1jmp1
            IF(dzqw(ij,l)*dzqw(ij,l+1).gt.0.) THEN
                dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
            ELSE
                dzq(ij,l)=0.
            ENDIF
            dzqmax=pente_max*min(adzqw(ij,l),adzqw(ij,l+1))
            dzq(ij,l)=sign(min(abs(dzq(ij,l)),dzqmax),dzq(ij,l))
         ENDDO
      ENDDO

      DO ij=1,ip1jmp1
         dzq(ij,1)=0.
         dzq(ij,llm)=0.
      ENDDO

c ---------------------------------------------------------------
c   .... calcul des termes d'advection verticale  .......
c ---------------------------------------------------------------

c calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour calculer dq

       DO l = 1,llm-1
         do  ij = 1,ip1jmp1
          IF(w(ij,l+1).gt.0.) THEN
             sigw=w(ij,l+1)/masse(ij,l+1)
             wq(ij,l+1)=w(ij,l+1)*(q(ij,l+1)+0.5*(1.-sigw)*dzq(ij,l+1))
          ELSE
             sigw=w(ij,l+1)/masse(ij,l)
             wq(ij,l+1)=w(ij,l+1)*(q(ij,l)-0.5*(1.+sigw)*dzq(ij,l))
          ENDIF
         ENDDO
       ENDDO

       DO ij=1,ip1jmp1
          wq(ij,llm+1)=0.
          wq(ij,1)=0.
       ENDDO

      DO l=1,llm
         DO ij=1,ip1jmp1
            newmasse=masse(ij,l)+w(ij,l+1)-w(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l))
     &         /newmasse
            masse(ij,l)=newmasse
         ENDDO
      ENDDO

      END
