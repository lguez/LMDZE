!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/pentes_ini.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE pentes_ini (q,w,masse,pbaru,pbarv,mode)
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      IMPLICIT NONE

c=======================================================================
c   Adaptation LMDZ:  A.Armengaud (LGGE)
c   ----------------
c
c   ********************************************************************
c   Transport des traceurs par la methode des pentes
c   ********************************************************************
c   Reference possible : Russel. G.L., Lerner J.A.:
c         A new Finite-Differencing Scheme for Traceur Transport 
c         Equation , Journal of Applied Meteorology, pp 1483-1498,dec. 81 
c   ********************************************************************
c   q,w,masse,pbaru et pbarv 
c                      sont des arguments d'entree  pour le s-pg ....
c
c=======================================================================



c   Arguments:
c   ----------
      integer mode
      REAL, intent(in):: pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm )
      REAL q( iip1,jjp1,llm,0:3)
      REAL w( ip1jmp1,llm )
      REAL masse( iip1,jjp1,llm)
c   Local:
c   ------
      LOGICAL limit
      REAL sm ( iip1,jjp1, llm )
      REAL s0( iip1,jjp1,llm ),  sx( iip1,jjp1,llm )
      REAL sy( iip1,jjp1,llm ),  sz( iip1,jjp1,llm )
      real masn,mass,zz
      INTEGER i,j,l,iq

c  modif Fred 24 03 96

      real sinlon(iip1),sinlondlon(iip1)
      real coslon(iip1),coslondlon(iip1)
      save sinlon,coslon,sinlondlon,coslondlon
      real dyn1,dyn2,dys1,dys2
      real qpn,qps,dqzpn,dqzps
      real smn,sms,s0n,s0s,sxn(iip1),sxs(iip1)
      real qmin,zq,pente_max
c
      REAL      SSUM
      integer ismax,ismin,lati,latf
      EXTERNAL  SSUM, convflu,ismin,ismax
      logical first
      save first
c   fin modif

c      EXTERNAL masskg
      EXTERNAL advx
      EXTERNAL advy
      EXTERNAL advz

c  modif Fred 24 03 96
      data first/.true./

      limit = .TRUE.
      pente_max=2
c     if (mode.eq.1.or.mode.eq.3) then
c     if (mode.eq.1) then
      if (mode.ge.1) then
        lati=2
        latf=jjm
      else
        lati=1
        latf=jjp1
      endif

      qmin=0.4995
      qmin=0.
      if(first) then
         print*,'SCHEMA AMONT NOUVEAU'
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
            print*,coslondlon(i),sinlondlon(i)
         enddo
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
         print*,'sum sinlondlon ',ssum(iim,sinlondlon,1)/sinlondlon(1)
         print*,'sum coslondlon ',ssum(iim,coslondlon,1)/coslondlon(1)
        DO l = 1,llm
        DO j = 1,jjp1
         DO i = 1,iip1
         q ( i,j,l,1 )=0.
         q ( i,j,l,2 )=0.
         q ( i,j,l,3 )=0.  
         ENDDO
         ENDDO
        ENDDO
        
      endif
c   Fin modif Fred

c *** q contient les qqtes de traceur avant l'advection 

c *** Affectation des tableaux S a partir de Q
c *** Rem : utilisation de SCOPY ulterieurement
 
       DO l = 1,llm
        DO j = 1,jjp1
         DO i = 1,iip1
             s0( i,j,llm+1-l ) = q ( i,j,l,0 )
             sx( i,j,llm+1-l ) = q ( i,j,l,1 )
             sy( i,j,llm+1-l ) = q ( i,j,l,2 )
             sz( i,j,llm+1-l ) = q ( i,j,l,3 )
         ENDDO
        ENDDO
       ENDDO

c      PRINT*,'----- S0 just before conversion -------'
c      PRINT*,'S0(16,12,1)=',s0(16,12,1) 
c      PRINT*,'Q(16,12,1,4)=',q(16,12,1,4)

c *** On calcule la masse d'air en kg

       DO  l = 1,llm
         DO  j = 1,jjp1
           DO  i = 1,iip1
            sm ( i,j,llm+1-l)=masse( i,j,l )
          ENDDO
         ENDDO
       ENDDO

c *** On converti les champs S en atome (resp. kg) 
c *** Les routines d'advection traitent les champs
c *** a advecter si ces derniers sont en atome (resp. kg)
c *** A optimiser !!!

       DO  l = 1,llm
         DO  j = 1,jjp1
           DO  i = 1,iip1
               s0(i,j,l) = s0(i,j,l) * sm ( i,j,l )
               sx(i,j,l) = sx(i,j,l) * sm ( i,j,l )
               sy(i,j,l) = sy(i,j,l) * sm ( i,j,l )
               sz(i,j,l) = sz(i,j,l) * sm ( i,j,l )
           ENDDO
         ENDDO
       ENDDO

c       ss0 = 0.
c       DO l = 1,llm
c        DO j = 1,jjp1
c         DO i = 1,iim
c            ss0 = ss0 + s0 ( i,j,l )
c         ENDDO
c        ENDDO
c       ENDDO
c       PRINT*, 'valeur tot s0 avant advection=',ss0

c *** Appel des subroutines d'advection en X, en Y et en Z
c *** Advection avec "time-splitting"
      
c-----------------------------------------------------------
c      PRINT*,'----- S0 just before ADVX -------'
c      PRINT*,'S0(16,12,1)=',s0(16,12,1)

c-----------------------------------------------------------
c      do l=1,llm
c         do j=1,jjp1
c          do i=1,iip1
c             zq=s0(i,j,l)/sm(i,j,l)
c            if(zq.lt.qmin)
c    ,       print*,'avant advx1, s0(',i,',',j,',',l,')=',zq
c          enddo
c         enddo
c      enddo
CCC
       if(mode.eq.2) then
          do l=1,llm
            s0s=0.
            s0n=0.
            dyn1=0.
            dys1=0.
            dyn2=0.
            dys2=0.
            smn=0.
            sms=0.
            do i=1,iim
               smn=smn+sm(i,1,l)
               sms=sms+sm(i,jjp1,l)
               s0n=s0n+s0(i,1,l)
               s0s=s0s+s0(i,jjp1,l)
               zz=sy(i,1,l)/sm(i,1,l)
               dyn1=dyn1+sinlondlon(i)*zz
               dyn2=dyn2+coslondlon(i)*zz
               zz=sy(i,jjp1,l)/sm(i,jjp1,l)
               dys1=dys1+sinlondlon(i)*zz
               dys2=dys2+coslondlon(i)*zz
            enddo
            do i=1,iim
               sy(i,1,l)=dyn1*sinlon(i)+dyn2*coslon(i)
               sy(i,jjp1,l)=dys1*sinlon(i)+dys2*coslon(i)
            enddo
            do i=1,iim
               s0(i,1,l)=s0n/smn+sy(i,1,l)
               s0(i,jjp1,l)=s0s/sms-sy(i,jjp1,l)
            enddo

            s0(iip1,1,l)=s0(1,1,l)
            s0(iip1,jjp1,l)=s0(1,jjp1,l)

            do i=1,iim
               sxn(i)=s0(i+1,1,l)-s0(i,1,l)
               sxs(i)=s0(i+1,jjp1,l)-s0(i,jjp1,l)
c   on rerentre les masses
            enddo
            do i=1,iim
               sy(i,1,l)=sy(i,1,l)*sm(i,1,l)
               sy(i,jjp1,l)=sy(i,jjp1,l)*sm(i,jjp1,l)
               s0(i,1,l)=s0(i,1,l)*sm(i,1,l)
               s0(i,jjp1,l)=s0(i,jjp1,l)*sm(i,jjp1,l)
            enddo
            sxn(iip1)=sxn(1)
            sxs(iip1)=sxs(1)
            do i=1,iim
               sx(i+1,1,l)=0.25*(sxn(i)+sxn(i+1))*sm(i+1,1,l)
               sx(i+1,jjp1,l)=0.25*(sxs(i)+sxs(i+1))*sm(i+1,jjp1,l)
            enddo
            s0(iip1,1,l)=s0(1,1,l)
            s0(iip1,jjp1,l)=s0(1,jjp1,l)
            sy(iip1,1,l)=sy(1,1,l)
            sy(iip1,jjp1,l)=sy(1,jjp1,l)
            sx(1,1,l)=sx(iip1,1,l)
            sx(1,jjp1,l)=sx(iip1,jjp1,l)
          enddo
      endif

      if (mode.eq.4) then
         do l=1,llm
            do i=1,iip1
               sx(i,1,l)=0.
               sx(i,jjp1,l)=0.
               sy(i,1,l)=0.
               sy(i,jjp1,l)=0.
            enddo
         enddo
      endif
      call limx(s0,sx,sm,pente_max)
c     call minmaxq(zq,1.e33,-1.e33,'avant advx     ')
       call advx( limit,.5*dtvr,pbaru,sm,s0,sx,sy,sz,lati,latf)
c     call minmaxq(zq,1.e33,-1.e33,'avant advy     ')
      if (mode.eq.4) then
         do l=1,llm
            do i=1,iip1
               sx(i,1,l)=0.
               sx(i,jjp1,l)=0.
               sy(i,1,l)=0.
               sy(i,jjp1,l)=0.
            enddo
         enddo
      endif
       call   limy(s0,sy,sm,pente_max)
       call advy( limit,.5*dtvr,pbarv,sm,s0,sx,sy,sz ) 
c     call minmaxq(zq,1.e33,-1.e33,'avant advz     ')
       do j=1,jjp1
          do i=1,iip1
             sz(i,j,1)=0.
             sz(i,j,llm)=0.
          enddo
       enddo
       call limz(s0,sz,sm,pente_max)
       call advz( limit,dtvr,w,sm,s0,sx,sy,sz )
      if (mode.eq.4) then
         do l=1,llm
            do i=1,iip1
               sx(i,1,l)=0.
               sx(i,jjp1,l)=0.
               sy(i,1,l)=0.
               sy(i,jjp1,l)=0.
            enddo
         enddo
      endif
        call limy(s0,sy,sm,pente_max)
       call advy( limit,.5*dtvr,pbarv,sm,s0,sx,sy,sz ) 
       do l=1,llm
          do j=1,jjp1
             sm(iip1,j,l)=sm(1,j,l)
             s0(iip1,j,l)=s0(1,j,l)
             sx(iip1,j,l)=sx(1,j,l)
             sy(iip1,j,l)=sy(1,j,l)
             sz(iip1,j,l)=sz(1,j,l)
          enddo
       enddo


c     call minmaxq(zq,1.e33,-1.e33,'avant advx     ')
      if (mode.eq.4) then
         do l=1,llm
            do i=1,iip1
               sx(i,1,l)=0.
               sx(i,jjp1,l)=0.
               sy(i,1,l)=0.
               sy(i,jjp1,l)=0.
            enddo
         enddo
      endif
       call limx(s0,sx,sm,pente_max)
       call advx( limit,.5*dtvr,pbaru,sm,s0,sx,sy,sz,lati,latf) 
c     call minmaxq(zq,1.e33,-1.e33,'apres advx     ')
c      do l=1,llm
c         do j=1,jjp1
c          do i=1,iip1
c             zq=s0(i,j,l)/sm(i,j,l)
c            if(zq.lt.qmin)
c    ,       print*,'apres advx2, s0(',i,',',j,',',l,')=',zq
c          enddo
c         enddo
c      enddo
c ***   On repasse les S dans la variable q directement 14/10/94
c   On revient a des rapports de melange en divisant par la masse

c En dehors des poles:

       DO  l = 1,llm
        DO  j = 1,jjp1
         DO  i = 1,iim
             q(i,j,llm+1-l,0)=s0(i,j,l)/sm(i,j,l)
             q(i,j,llm+1-l,1)=sx(i,j,l)/sm(i,j,l)
             q(i,j,llm+1-l,2)=sy(i,j,l)/sm(i,j,l)
             q(i,j,llm+1-l,3)=sz(i,j,l)/sm(i,j,l)
         ENDDO
        ENDDO
      ENDDO

c Traitements specifiques au pole

      if(mode.ge.1) then
      DO l=1,llm
c   filtrages aux poles
         masn=ssum(iim,sm(1,1,l),1)
         mass=ssum(iim,sm(1,jjp1,l),1)
         qpn=ssum(iim,s0(1,1,l),1)/masn
         qps=ssum(iim,s0(1,jjp1,l),1)/mass
         dqzpn=ssum(iim,sz(1,1,l),1)/masn
         dqzps=ssum(iim,sz(1,jjp1,l),1)/mass
         do i=1,iip1
            q( i,1,llm+1-l,3)=dqzpn
            q( i,jjp1,llm+1-l,3)=dqzps
            q( i,1,llm+1-l,0)=qpn
            q( i,jjp1,llm+1-l,0)=qps
         enddo
         if(mode.eq.3) then
            dyn1=0.
            dys1=0.
            dyn2=0.
            dys2=0.
            do i=1,iim
               dyn1=dyn1+sinlondlon(i)*sy(i,1,l)/sm(i,1,l)
               dyn2=dyn2+coslondlon(i)*sy(i,1,l)/sm(i,1,l)
               dys1=dys1+sinlondlon(i)*sy(i,jjp1,l)/sm(i,jjp1,l)
               dys2=dys2+coslondlon(i)*sy(i,jjp1,l)/sm(i,jjp1,l)
            enddo
            do i=1,iim
               q(i,1,llm+1-l,2)=
     s          (sinlon(i)*dyn1+coslon(i)*dyn2)
               q(i,1,llm+1-l,0)=q(i,1,llm+1-l,0)+q(i,1,llm+1-l,2)
               q(i,jjp1,llm+1-l,2)=
     s          (sinlon(i)*dys1+coslon(i)*dys2)
               q(i,jjp1,llm+1-l,0)=q(i,jjp1,llm+1-l,0)
     s         -q(i,jjp1,llm+1-l,2)
            enddo
         endif
         if(mode.eq.1) then
c   on filtre les valeurs au bord de la "grande maille pole"
            dyn1=0.
            dys1=0.
            dyn2=0.
            dys2=0.
            do i=1,iim
               zz=s0(i,2,l)/sm(i,2,l)-q(i,1,llm+1-l,0)
               dyn1=dyn1+sinlondlon(i)*zz
               dyn2=dyn2+coslondlon(i)*zz
               zz=q(i,jjp1,llm+1-l,0)-s0(i,jjm,l)/sm(i,jjm,l)
               dys1=dys1+sinlondlon(i)*zz
               dys2=dys2+coslondlon(i)*zz
            enddo
            do i=1,iim
               q(i,1,llm+1-l,2)=
     s          (sinlon(i)*dyn1+coslon(i)*dyn2)/2.
               q(i,1,llm+1-l,0)=q(i,1,llm+1-l,0)+q(i,1,llm+1-l,2)
               q(i,jjp1,llm+1-l,2)=
     s          (sinlon(i)*dys1+coslon(i)*dys2)/2.
               q(i,jjp1,llm+1-l,0)=q(i,jjp1,llm+1-l,0)
     s         -q(i,jjp1,llm+1-l,2)
            enddo
            q(iip1,1,llm+1-l,0)=q(1,1,llm+1-l,0)
            q(iip1,jjp1,llm+1-l,0)=q(1,jjp1,llm+1-l,0)

            do i=1,iim
               sxn(i)=q(i+1,1,llm+1-l,0)-q(i,1,llm+1-l,0)
               sxs(i)=q(i+1,jjp1,llm+1-l,0)-q(i,jjp1,llm+1-l,0)
            enddo
            sxn(iip1)=sxn(1)
            sxs(iip1)=sxs(1)
            do i=1,iim
               q(i+1,1,llm+1-l,1)=0.25*(sxn(i)+sxn(i+1))
               q(i+1,jjp1,llm+1-l,1)=0.25*(sxs(i)+sxs(i+1))
            enddo
            q(1,1,llm+1-l,1)=q(iip1,1,llm+1-l,1)
            q(1,jjp1,llm+1-l,1)=q(iip1,jjp1,llm+1-l,1)

         endif

       ENDDO
       endif

c bouclage en longitude
      do iq=0,3
         do l=1,llm
            do j=1,jjp1
               q(iip1,j,l,iq)=q(1,j,l,iq)
            enddo
         enddo
      enddo

c       PRINT*, ' SORTIE DE PENTES ---  ca peut glisser ....'

        DO l = 1,llm
    	 DO j = 1,jjp1
    	  DO i = 1,iip1
                IF (q(i,j,l,0).lt.0.)  THEN
c                    PRINT*,'------------ BIP-----------' 
c                    PRINT*,'Q0(',i,j,l,')=',q(i,j,l,0)
c                    PRINT*,'QX(',i,j,l,')=',q(i,j,l,1)
c                    PRINT*,'QY(',i,j,l,')=',q(i,j,l,2)
c                    PRINT*,'QZ(',i,j,l,')=',q(i,j,l,3)
c       		     PRINT*,' PBL EN SORTIE DE PENTES'
                     q(i,j,l,0)=0.
c                    STOP
                 ENDIF
          ENDDO
         ENDDO
        ENDDO

c       PRINT*, '-------------------------------------------'
        
       do l=1,llm
          do j=1,jjp1
           do i=1,iip1
             if(q(i,j,l,0).lt.qmin)
     ,       print*,'apres pentes, s0(',i,',',j,',',l,')=',q(i,j,l,0)
           enddo
          enddo
       enddo
      RETURN
      END












