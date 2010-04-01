!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/prather.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      SUBROUTINE prather (q,w,masse,pbaru,pbarv,nt,dt)
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
c   ************************************************
c   Transport des traceurs par la methode de prather
c   Ref : 
c
c   ************************************************
c   q,w,pext,pbaru et pbarv : arguments d'entree  pour le s-pg
c
c=======================================================================



c   Arguments:
c   ----------
      INTEGER iq,nt
      REAL, intent(in):: pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm )
      REAL masse(iip1,jjp1,llm)
      REAL q( iip1,jjp1,llm,0:9)
      REAL w( ip1jmp1,llm )
      integer ordre,ilim

c   Local:
c   ------
      LOGICAL limit
      real zq(iip1,jjp1,llm)
      REAL sm ( iip1,jjp1, llm )
      REAL s0( iip1,jjp1,llm ),  sx( iip1,jjp1,llm )
      REAL sy( iip1,jjp1,llm ),  sz( iip1,jjp1,llm )
      REAL sxx( iip1,jjp1,llm)
      REAL sxy( iip1,jjp1,llm)
      REAL sxz( iip1,jjp1,llm)
      REAL syy( iip1,jjp1,llm )
      REAL syz( iip1,jjp1,llm )
      REAL szz( iip1,jjp1,llm ),zz
      INTEGER i,j,l,indice
      real sxn(iip1),sxs(iip1)

      real sinlon(iip1),sinlondlon(iip1)
      real coslon(iip1),coslondlon(iip1)
      real qmin,qmax
      save qmin,qmax
      save sinlon,coslon,sinlondlon,coslondlon
      real dyn1,dyn2,dys1,dys2,qpn,qps,dqzpn,dqzps
      real masn,mass
c
      REAL      SSUM
      integer ismax,ismin
      EXTERNAL  SSUM, convflu,ismin,ismax
      logical first
      save first
      EXTERNAL advxp,advyp,advzp 


      data first/.true./
      data qmin,qmax/-1.e33,1.e33/


c==========================================================================
c==========================================================================
c     MODIFICATION POUR PAS DE TEMPS ADAPTATIF, dtvr remplace par dt
c==========================================================================
c==========================================================================
      REAL dt
c==========================================================================
      limit = .TRUE.
 
      if(first) then
         print*,'SCHEMA PRATHER'
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

        DO l = 1,llm
        DO j = 1,jjp1
        DO i = 1,iip1
        q( i,j,l,1 )=0.
        q( i,j,l,2)=0.
        q( i,j,l,3)=0.
        q( i,j,l,4)=0.
        q( i,j,l,5)=0.
        q( i,j,l,6)=0.
        q( i,j,l,7)=0.
        q( i,j,l,8)=0.
        q( i,j,l,9)=0.
        ENDDO
        ENDDO
        ENDDO
      endif
c   Fin modif Fred

c *** On calcule la masse d'air en kg

       DO l = 1,llm
        DO j = 1,jjp1
         DO i = 1,iip1
         sm( i,j,llm+1-l ) =masse(i,j,l)
         ENDDO
        ENDDO
       ENDDO

c *** q contient les qqtes de traceur avant l'advection 

c *** Affectation des tableaux S a partir de Q
 
       DO l = 1,llm
        DO j = 1,jjp1
         DO i = 1,iip1
       s0( i,j,l) = q ( i,j,llm+1-l,0 )*sm(i,j,l)
       sx( i,j,l) = q( i,j,llm+1-l,1 )*sm(i,j,l)
       sy( i,j,l) = q( i,j,llm+1-l,2)*sm(i,j,l)
       sz( i,j,l) = q( i,j,llm+1-l,3)*sm(i,j,l)
       sxx( i,j,l) = q( i,j,llm+1-l,4)*sm(i,j,l)
       sxy( i,j,l) = q( i,j,llm+1-l,5)*sm(i,j,l)
       sxz( i,j,l) = q( i,j,llm+1-l,6)*sm(i,j,l)
       syy( i,j,l) = q( i,j,llm+1-l,7)*sm(i,j,l)
       syz( i,j,l) = q( i,j,llm+1-l,8)*sm(i,j,l)
       szz( i,j,l) = q( i,j,llm+1-l,9)*sm(i,j,l)
         ENDDO
        ENDDO
       ENDDO
c *** Appel des subroutines d'advection en X, en Y et en Z
c *** Advection avec "time-splitting"
      
c-----------------------------------------------------------
       do indice =1,nt
       call advxp( limit,0.5*dt,pbaru,sm,s0,sx,sy,sz
     .             ,sxx,sxy,sxz,syy,syz,szz,1 )
        end do
        do l=1,llm
        do i=1,iip1
        sy(i,1,l)=0.
        sy(i,jjp1,l)=0.
        enddo
        enddo
c---------------------------------------------------------
       call advyp( limit,.5*dt*nt,pbarv,sm,s0,sx,sy,sz
     .             ,sxx,sxy,sxz,syy,syz,szz,1 )
c---------------------------------------------------------

c---------------------------------------------------------
       do j=1,jjp1
          do i=1,iip1
             sz(i,j,1)=0.
             sz(i,j,llm)=0.
             sxz(i,j,1)=0.
             sxz(i,j,llm)=0.
             syz(i,j,1)=0.
             syz(i,j,llm)=0.
             szz(i,j,1)=0.
             szz(i,j,llm)=0.
          enddo
       enddo
       call advzp( limit,dt*nt,w,sm,s0,sx,sy,sz 
     .             ,sxx,sxy,sxz,syy,syz,szz,1 )
        do l=1,llm
        do i=1,iip1
        sy(i,1,l)=0.
        sy(i,jjp1,l)=0.
        enddo
        enddo

c---------------------------------------------------------

c---------------------------------------------------------
       call advyp( limit,.5*dt*nt,pbarv,sm,s0,sx,sy,sz
     .             ,sxx,sxy,sxz,syy,syz,szz,1 )
c---------------------------------------------------------
       DO l = 1,llm
        DO j = 1,jjp1
             s0( iip1,j,l)=s0( 1,j,l )
             sx( iip1,j,l)=sx( 1,j,l )
             sy( iip1,j,l)=sy( 1,j,l )
             sz( iip1,j,l)=sz( 1,j,l )
             sxx( iip1,j,l)=sxx( 1,j,l )
             sxy( iip1,j,l)=sxy( 1,j,l) 
             sxz( iip1,j,l)=sxz( 1,j,l )
             syy( iip1,j,l)=syy( 1,j,l )
             syz( iip1,j,l)=syz( 1,j,l)
             szz( iip1,j,l)=szz( 1,j,l )
        ENDDO
       ENDDO
       do indice=1,nt
       call advxp( limit,0.5*dt,pbaru,sm,s0,sx,sy,sz
     .             ,sxx,sxy,sxz,syy,syz,szz,1 )
        end do
c---------------------------------------------------------
c---------------------------------------------------------
c ***   On repasse les S dans la variable qpr
c ***   On repasse les S dans la variable q directement 14/10/94

       DO  l = 1,llm
        DO  j = 1,jjp1
         DO  i = 1,iip1
      q( i,j,llm+1-l,0 )=s0( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,1 ) = sx( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,2 ) = sy( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,3 ) = sz( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,4 ) = sxx( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,5 ) = sxy( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,6 ) = sxz( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,7 ) = syy( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,8 ) = syz( i,j,l )/sm(i,j,l)
      q( i,j,llm+1-l,9 ) = szz( i,j,l )/sm(i,j,l)
      ENDDO
      ENDDO
      ENDDO

c---------------------------------------------------------
c      go to  777
c   filtrages aux poles

c Traitements specifiques au pole

c   filtrages aux poles
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
c       enddo
c         print*,'qpn',qpn,'qps',qps
c          print*,'dqzpn',dqzpn,'dqzps',dqzps
c       enddo
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
     $   (sinlon(i)*dyn1+coslon(i)*dyn2)/2.
         q(i,1,llm+1-l,0)=q(i,1,llm+1-l,0)
     $          +q(i,1,llm+1-l,2)
         q(i,jjp1,llm+1-l,2)=
     $   (sinlon(i)*dys1+coslon(i)*dys2)/2.
         q(i,jjp1,llm+1-l,0)=q(i,jjp1,llm+1-l,0)
     $      -q(i,jjp1,llm+1-l,2)
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
      END DO
      q(1,1,llm+1-l,1)=q(iip1,1,llm+1-l,1)
      q(1,jjp1,llm+1-l,1)=
     $   q(iip1,jjp1,llm+1-l,1)
        enddo
         do l=1,llm
           do i=1,iim
            q( i,1,llm+1-l,4)=0.
            q( i,jjp1,llm+1-l,4)=0.
            q( i,1,llm+1-l,5)=0.
            q( i,jjp1,llm+1-l,5)=0.
            q( i,1,llm+1-l,6)=0.
            q( i,jjp1,llm+1-l,6)=0.
            q( i,1,llm+1-l,7)=0.
            q( i,jjp1,llm+1-l,7)=0.
            q( i,1,llm+1-l,8)=0.
            q( i,jjp1,llm+1-l,8)=0.
            q( i,1,llm+1-l,9)=0.
            q( i,jjp1,llm+1-l,9)=0.
          enddo
         ENDDO

777      continue
c
c   bouclage en longitude
      do l=1,llm
      do j=1,jjp1
      q(iip1,j,l,0)=q(1,j,l,0)
      q(iip1,j,llm+1-l,0)=q(1,j,llm+1-l,0)
      q(iip1,j,llm+1-l,1)=q(1,j,llm+1-l,1)
      q(iip1,j,llm+1-l,2)=q(1,j,llm+1-l,2)
      q(iip1,j,llm+1-l,3)=q(1,j,llm+1-l,3)
      q(iip1,j,llm+1-l,4)=q(1,j,llm+1-l,4)
      q(iip1,j,llm+1-l,5)=q(1,j,llm+1-l,5)
      q(iip1,j,llm+1-l,6)=q(1,j,llm+1-l,6)
      q(iip1,j,llm+1-l,7)=q(1,j,llm+1-l,7)
      q(iip1,j,llm+1-l,8)=q(1,j,llm+1-l,8)
      q(iip1,j,llm+1-l,9)=q(1,j,llm+1-l,9)
      enddo
      enddo
        DO l = 1,llm
    	 DO j = 2,jjm
           DO i = 1,iip1
         IF (q(i,j,l,0).lt.0.)  THEN
         PRINT*,'------------ BIP-----------' 
         PRINT*,'S0(',i,j,l,')=',q(i,j,l,0),
     $          q(i,j-1,l,0)
         PRINT*,'SX(',i,j,l,')=',q(i,j,l,1)
         PRINT*,'SY(',i,j,l,')=',q(i,j,l,2),
     $   q(i,j-1,l,2)   
         PRINT*,'SZ(',i,j,l,')=',q(i,j,l,3)
c    		     PRINT*,' PBL EN SORTIE D'' ADVZP'
                     q(i,j,l,0)=0.
c                  STOP
               ENDIF
           ENDDO
         ENDDO
         do j=1,jjp1,jjm
         do i=1,iip1
               IF (q(i,j,l,0).lt.0.)  THEN
               PRINT*,'------------ BIP 2-----------'
         PRINT*,'S0(',i,j,l,')=',q(i,j,l,0)
         PRINT*,'SX(',i,j,l,')=',q(i,j,l,1)
         PRINT*,'SY(',i,j,l,')=',q(i,j,l,2)
         PRINT*,'SZ(',i,j,l,')=',q(i,j,l,3)

                     q(i,j,l,0)=0.
c                  STOP
               ENDIF
         enddo
         enddo
        ENDDO
      RETURN
      END
