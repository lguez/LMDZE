!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/vlspltqs.F,v 1.2 2005/02/24 12:16:57 fairhead Exp $
!
       SUBROUTINE vlspltqs ( q,pente_max,masse,w,pbaru,pbarv,pdt,
     ,                                  p,pk,teta                 )
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget, F.Codron 
c
c    ********************************************************************
c          Shema  d'advection " pseudo amont " .
c      + test sur humidite specifique: Q advecte< Qsat aval
c                   (F. Codron, 10/99)
c    ********************************************************************
c     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c     pente_max facteur de limitation des pentes: 2 en general
c                                                0 pour un schema amont
c     pbaru,pbarv,w flux de masse en u ,v ,w
c     pdt pas de temps
c
c     teta temperature potentielle, p pression aux interfaces,
c     pk exner au milieu des couches necessaire pour calculer Qsat
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
      REAL, intent(in):: pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm)
      real, intent(in):: pdt
      REAL, intent(in):: p(ip1jmp1,llmp1)
      real teta(ip1jmp1,llm),pk(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
c
      REAL qsat(ip1jmp1,llm)
      REAL zm(ip1jmp1,llm)
      REAL mu(ip1jmp1,llm)
      REAL mv(ip1jm,llm)
      REAL mw(ip1jmp1,llm+1)
      REAL zq(ip1jmp1,llm)
      REAL temps1,temps2,temps3
      REAL zzpbar, zzw
      LOGICAL testcpu
      SAVE testcpu
      SAVE temps1,temps2,temps3

      REAL qmin,qmax
      DATA qmin,qmax/0.,1.e33/
      DATA testcpu/.false./
      DATA temps1,temps2,temps3/0.,0.,0./

c--pour rapport de melange saturant--

      REAL rtt,retv,r2es,r3les,r3ies,r4les,r4ies,play
      REAL ptarg,pdelarg,foeew,zdelta
      REAL tempe(ip1jmp1)

c    fonction psat(T)

       FOEEW ( PTARG,PDELARG ) = EXP (
     *          (R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT)
     * / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )

        r2es  = 380.11733 
        r3les = 17.269
        r3ies = 21.875
        r4les = 35.86
        r4ies = 7.66
        retv = 0.6077667
        rtt  = 273.16

c-- Calcul de Qsat en chaque point
c-- approximation: au milieu des couches play(l)=(p(l)+p(l+1))/2
c   pour eviter une exponentielle.
        DO l = 1, llm
         DO ij = 1, ip1jmp1
          tempe(ij) = teta(ij,l) * pk(ij,l) /cpp
         ENDDO
         DO ij = 1, ip1jmp1
          zdelta = MAX( 0., SIGN(1., rtt - tempe(ij)) )
          play   = 0.5*(p(ij,l)+p(ij,l+1))
          qsat(ij,l) = MIN(0.5, r2es* FOEEW(tempe(ij),zdelta) / play )
          qsat(ij,l) = qsat(ij,l) / ( 1. - retv * qsat(ij,l) )
         ENDDO
        ENDDO

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

c      call minmaxq(zq,qmin,qmax,'avant vlxqs     ')
      call vlxqs(zq,pente_max,zm,mu,qsat)


c     call minmaxq(zq,qmin,qmax,'avant vlyqs     ')

      call vlyqs(zq,pente_max,zm,mv,qsat)


c      call minmaxq(zq,qmin,qmax,'avant vlz     ')

      call vlz(zq,pente_max,zm,mw)


c     call minmaxq(zq,qmin,qmax,'avant vlyqs     ')
c     call minmaxq(zm,qmin,qmax,'M avant vlyqs     ')

      call vlyqs(zq,pente_max,zm,mv,qsat)


c     call minmaxq(zq,qmin,qmax,'avant vlxqs     ')
c     call minmaxq(zm,qmin,qmax,'M avant vlxqs     ')

      call vlxqs(zq,pente_max,zm,mu,qsat)

c     call minmaxq(zq,qmin,qmax,'apres vlxqs     ')
c     call minmaxq(zm,qmin,qmax,'M apres vlxqs     ')


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
      SUBROUTINE vlxqs(q,pente_max,masse,u_m,qsat)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
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
      REAL u_m( ip1jmp1,llm )
      REAL q(ip1jmp1,llm)
      REAL qsat(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      INTEGER n0,iadvplus(ip1jmp1,llm),nl(llm)
c
      REAL new_m,zu_m,zdum(ip1jmp1,llm)
      REAL dxq(ip1jmp1,llm),dxqu(ip1jmp1)
      REAL zz(ip1jmp1)
      REAL adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
      REAL u_mq(ip1jmp1,llm)

      Logical first,testcpu
      SAVE first,testcpu

      REAL      SSUM
      REAL temps0,temps1,temps2,temps3,temps4,temps5
      SAVE temps0,temps1,temps2,temps3,temps4,temps5


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
c     IF (pente_max.gt.10) THEN

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


c   calcul des flux a gauche et a droite

c   on cumule le flux correspondant a toutes les mailles dont la masse
c   au travers de la paroi pENDant le pas de temps.
c   le rapport de melange de l'air advecte est min(q_vanleer, Qsat_downwind)
      DO l=1,llm
       DO ij=iip2,ip1jm-1
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
             u_mq(ij,l)=u_m(ij,l)*
     $         min(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l),qsat(ij+1,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
             u_mq(ij,l)=u_m(ij,l)*
     $         min(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l),qsat(ij,l))
          ENDIF
       ENDDO
      ENDDO


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
      DO l=1,llm
       DO ij=iip1+iip1,ip1jm,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO



c   traitement special pour le cas ou on advecte en longitude plus que le
c   contenu de la maille.
c   cette partie est mal vectorisee.

c   pas d'influence de la pression saturante (pour l'instant)

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



c   bouclage en latitude

      DO l=1,llm
        DO ij=iip1+iip1,ip1jm,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO


c   calcul des tendances

      DO l=1,llm
         DO ij=iip2+1,ip1jm
            new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+
     &      u_mq(ij-1,l)-u_mq(ij,l))
     &      /new_m
            masse(ij,l)=new_m
         ENDDO
c   Modif Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         DO ij=iip1+iip1,ip1jm,iip1
            q(ij-iim,l)=q(ij,l)
            masse(ij-iim,l)=masse(ij,l)
         ENDDO
      ENDDO

c     CALL SCOPY((jjm-1)*llm,q(iip1+iip1,1),iip1,q(iip2,1),iip1)
c     CALL SCOPY((jjm-1)*llm,masse(iip1+iip1,1),iip1,masse(iip2,1),iip1)


      RETURN
      END
      SUBROUTINE vlyqs(q,pente_max,masse,masse_adv_v,qsat)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
c     qsat 	       est   un argument de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
c
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      IMPLICIT NONE
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL masse_adv_v( ip1jm,llm)
      REAL q(ip1jmp1,llm)
      REAL qsat(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      REAL dyq(ip1jmp1,llm),dyqv(ip1jm)
      REAL adyqv(ip1jm),dyqmax(ip1jmp1)
      REAL qbyv(ip1jm,llm)

      REAL qpns,qpsn,dyn1,dys1,dyn2,dys2,newmasse,fn,fs
c     REAL newq,oldmasse
      Logical first,testcpu
      REAL temps0,temps1,temps2,temps3,temps4,temps5
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
      SAVE first,testcpu

      REAL convpn,convps,convmpn,convmps
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
         IF( masse_adv_v(ij,l).GT.0. ) THEN
           qbyv(ij,l)= MIN( qsat(ij+iip1,l), q(ij+iip1,l )  +
     ,      dyq(ij+iip1,l)*0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l)))
         ELSE
              qbyv(ij,l)= MIN( qsat(ij,l), q(ij,l) - dyq(ij,l) *
     ,                   0.5*(1.+masse_adv_v(ij,l)/masse(ij,l)) )
         ENDIF
          qbyv(ij,l) = masse_adv_v(ij,l)*qbyv(ij,l)
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
         convpn=SSUM(iim,qbyv(1,l),1)/apoln
         convmpn=ssum(iim,masse_adv_v(1,l),1)/apoln
         DO ij = 1,iip1
            newmasse=masse(ij,l)+convmpn*aire(ij)
            q(ij,l)=(q(ij,l)*masse(ij,l)+convpn*aire(ij))/
     &               newmasse
            masse(ij,l)=newmasse
         ENDDO
         convps  = -SSUM(iim,qbyv(ip1jm-iim,l),1)/apols
         convmps = -SSUM(iim,masse_adv_v(ip1jm-iim,l),1)/apols
         DO ij = ip1jm+1,ip1jmp1
            newmasse=masse(ij,l)+convmps*aire(ij)
            q(ij,l)=(q(ij,l)*masse(ij,l)+convps*aire(ij))/
     &               newmasse
            masse(ij,l)=newmasse
         ENDDO
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
