!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/advn.F,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
      SUBROUTINE advn(q,masse,w,pbaru,pbarv,pdt,mode)
c
c     Auteur : F. Hourdin
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c   pbaru,pbarv,w flux de masse en u ,v ,w
c   pdt pas de temps
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use conf_gcm_m
      use comgeom
      IMPLICIT NONE
c

c
c   Arguments:
c   ----------
      integer mode
      real masse(ip1jmp1,llm)
      REAL, intent(in):: pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm),pdt
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
      integer ijlqmin,iqmin,jqmin,lqmin
      integer ismin
c
      real zm(ip1jmp1,llm),newmasse
      real mu(ip1jmp1,llm)
      real mv(ip1jm,llm)
      real mw(ip1jmp1,llm+1)
      real zq(ip1jmp1,llm),zz,qpn,qps
      real zqg(ip1jmp1,llm),zqd(ip1jmp1,llm)
      real zqs(ip1jmp1,llm),zqn(ip1jmp1,llm)
      real zqh(ip1jmp1,llm),zqb(ip1jmp1,llm)
      real temps0,temps1,temps2,temps3
      real ztemps1,ztemps2,ztemps3,ssum
      logical testcpu
      save testcpu
      save temps1,temps2,temps3
      real zzpbar,zzw

      real qmin,qmax
      data qmin,qmax/0.,1./
      data testcpu/.false./
      data temps1,temps2,temps3/0.,0.,0./

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

      do l=1,llm
         qpn=0.
         qps=0.
         do ij=1,iim
            qpn=qpn+q(ij,l)*masse(ij,l)
            qps=qps+q(ip1jm+ij,l)*masse(ip1jm+ij,l)
         enddo
         qpn=qpn/ssum(iim,masse(1,l),1)
         qps=qps/ssum(iim,masse(ip1jm+1,l),1)
         do ij=1,iip1
            q(ij,l)=qpn
            q(ip1jm+ij,l)=qps
         enddo
      enddo

      do ij=1,ip1jmp1
         mw(ij,llm+1)=0.
      enddo
      do l=1,llm
         do ij=1,ip1jmp1
            zq(ij,l)=q(ij,l)
            zm(ij,l)=masse(ij,l)
         enddo
      enddo

c     call minmaxq(zq,qmin,qmax,'avant vlx     ')
      call advnqx(zq,zqg,zqd)
      call advnx(zq,zqg,zqd,zm,mu,mode)
      call advnqy(zq,zqs,zqn)
      call advny(zq,zqs,zqn,zm,mv)
      call advnqz(zq,zqh,zqb)
      call advnz(zq,zqh,zqb,zm,mw)
c     call vlz(zq,0.,zm,mw)
      call advnqy(zq,zqs,zqn)
      call advny(zq,zqs,zqn,zm,mv)
      call advnqx(zq,zqg,zqd)
      call advnx(zq,zqg,zqd,zm,mu,mode)
c     call minmaxq(zq,qmin,qmax,'apres vlx     ')

      do l=1,llm
         do ij=1,ip1jmp1
           q(ij,l)=zq(ij,l)
         enddo
         do ij=1,ip1jm+1,iip1
            q(ij+iim,l)=q(ij,l)
         enddo
      enddo

      RETURN
      END

      SUBROUTINE advnqx(q,qg,qd)
c
c     Auteurs:   Calcul des valeurs de q aux point u.
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use conf_gcm_m
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real q(ip1jmp1,llm),qg(ip1jmp1,llm),qd(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l
c
      real dxqu(ip1jmp1),zqu(ip1jmp1)
      real zqmax(ip1jmp1),zqmin(ip1jmp1)
      logical extremum(ip1jmp1)

      integer mode
      save mode
      data mode/1/

c   calcul des pentes en u:
c   -----------------------
      if (mode.eq.0) then
         do l=1,llm
            do ij=1,ip1jm
               qd(ij,l)=q(ij,l)
               qg(ij,l)=q(ij,l)
            enddo
         enddo
      else
      do l = 1, llm
         do ij=iip2,ip1jm-1
            dxqu(ij)=q(ij+1,l)-q(ij,l)
            zqu(ij)=0.5*(q(ij+1,l)+q(ij,l))
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            dxqu(ij)=dxqu(ij-iim)
            zqu(ij)=zqu(ij-iim)
         enddo
         do ij=iip2,ip1jm-1
            zqu(ij)=zqu(ij)-dxqu(ij+1)/12.
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            zqu(ij)=zqu(ij-iim)
         enddo
         do ij=iip2+1,ip1jm
            zqu(ij)=zqu(ij)+dxqu(ij-1)/12.
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            zqu(ij-iim)=zqu(ij)
         enddo

c   calcul des valeurs max et min acceptees aux interfaces

         do ij=iip2,ip1jm-1
            zqmax(ij)=max(q(ij+1,l),q(ij,l))
            zqmin(ij)=min(q(ij+1,l),q(ij,l))
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            zqmax(ij)=zqmax(ij-iim)
            zqmin(ij)=zqmin(ij-iim)
         enddo
         do ij=iip2+1,ip1jm
            extremum(ij)=dxqu(ij)*dxqu(ij-1).le.0.
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            extremum(ij-iim)=extremum(ij)
         enddo
         do ij=iip2,ip1jm
            zqu(ij)=min(max(zqmin(ij),zqu(ij)),zqmax(ij))
         enddo
         do ij=iip2+1,ip1jm
            if(extremum(ij)) then
               qg(ij,l)=q(ij,l)
               qd(ij,l)=q(ij,l)
            else
               qd(ij,l)=zqu(ij)
               qg(ij,l)=zqu(ij-1)
            endif
         enddo
         do ij=iip1+iip1,ip1jm,iip1
            qd(ij-iim,l)=qd(ij,l)
            qg(ij-iim,l)=qg(ij,l)
         enddo

         goto 8888

         do ij=iip2+1,ip1jm
            if(extremum(ij).and..not.extremum(ij-1))
     s         qd(ij-1,l)=q(ij,l)
         enddo

         do ij=iip1+iip1,ip1jm,iip1
            qd(ij-iim,l)=qd(ij,l)
         enddo
         do ij=iip2,ip1jm-1
            if (extremum(ij).and..not.extremum(ij+1))
     s         qg(ij+1,l)=q(ij,l)
         enddo

         do ij=iip1+iip1,ip1jm,iip1
            qg(ij,l)=qg(ij-iim,l)
         enddo
8888     continue
      enddo
      endif
      RETURN
      END
      SUBROUTINE advnqy(q,qs,qn)
c
c     Auteurs:   Calcul des valeurs de q aux point v.
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use conf_gcm_m
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real q(ip1jmp1,llm),qs(ip1jmp1,llm),qn(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l
c
      real dyqv(ip1jm),zqv(ip1jm,llm)
      real zqmax(ip1jm),zqmin(ip1jm)
      logical extremum(ip1jmp1)

      integer mode
      save mode
      data mode/1/

      if (mode.eq.0) then
         do l=1,llm
            do ij=1,ip1jmp1
               qn(ij,l)=q(ij,l)
               qs(ij,l)=q(ij,l)
            enddo
         enddo
      else

c   calcul des pentes en u:
c   -----------------------
      do l = 1, llm
         do ij=1,ip1jm
            dyqv(ij)=q(ij,l)-q(ij+iip1,l)
         enddo

         do ij=iip2,ip1jm-iip1
            zqv(ij,l)=0.5*(q(ij+iip1,l)+q(ij,l))
            zqv(ij,l)=zqv(ij,l)+(dyqv(ij+iip1)-dyqv(ij-iip1))/12.
         enddo

         do ij=iip2,ip1jm
            extremum(ij)=dyqv(ij)*dyqv(ij-iip1).le.0.
         enddo

c Pas de pentes aux poles
         do ij=1,iip1
            zqv(ij,l)=q(ij,l)
            zqv(ip1jm-iip1+ij,l)=q(ip1jm+ij,l)
            extremum(ij)=.true.
            extremum(ip1jmp1-iip1+ij)=.true.
         enddo

c   calcul des valeurs max et min acceptees aux interfaces
         do ij=1,ip1jm
            zqmax(ij)=max(q(ij+iip1,l),q(ij,l))
            zqmin(ij)=min(q(ij+iip1,l),q(ij,l))
         enddo

         do ij=1,ip1jm
            zqv(ij,l)=min(max(zqmin(ij),zqv(ij,l)),zqmax(ij))
         enddo

         do ij=iip2,ip1jm
            if(extremum(ij)) then
               qs(ij,l)=q(ij,l)
               qn(ij,l)=q(ij,l)
c              if (.not.extremum(ij-iip1)) qs(ij-iip1,l)=q(ij,l)
c              if (.not.extremum(ij+iip1)) qn(ij+iip1,l)=q(ij,l)
            else
               qs(ij,l)=zqv(ij,l)
               qn(ij,l)=zqv(ij-iip1,l)
            endif
         enddo

         do ij=1,iip1
            qs(ij,l)=q(ij,l)
            qn(ij,l)=q(ij,l)
            qs(ip1jm+ij,l)=q(ip1jm+ij,l)
            qn(ip1jm+ij,l)=q(ip1jm+ij,l)
         enddo

      enddo
      endif
      RETURN
      END

      SUBROUTINE advnqz(q,qh,qb)
c
c     Auteurs:   Calcul des valeurs de q aux point v.
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use conf_gcm_m
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real q(ip1jmp1,llm),qh(ip1jmp1,llm),qb(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l
c
      real dzqw(ip1jmp1,llm+1),zqw(ip1jmp1,llm+1)
      real zqmax(ip1jmp1,llm),zqmin(ip1jmp1,llm)
      logical extremum(ip1jmp1,llm)

      integer mode
      save mode

      data mode/1/

c   calcul des pentes en u:
c   -----------------------

      if (mode.eq.0) then
         do l=1,llm
            do ij=1,ip1jmp1
               qb(ij,l)=q(ij,l)
               qh(ij,l)=q(ij,l)
            enddo
         enddo
      else
      do l = 2, llm
         do ij=1,ip1jmp1
            dzqw(ij,l)=q(ij,l-1)-q(ij,l)
            zqw(ij,l)=0.5*(q(ij,l-1)+q(ij,l))
         enddo
      enddo
      do ij=1,ip1jmp1
         dzqw(ij,1)=0.
         dzqw(ij,llm+1)=0.
      enddo
      do l=2,llm
         do ij=1,ip1jmp1
            zqw(ij,l)=zqw(ij,l)+(dzqw(ij,l+1)-dzqw(ij,l-1))/12.
         enddo
      enddo
      do l=2,llm-1
         do ij=1,ip1jmp1
            extremum(ij,l)=dzqw(ij,l)*dzqw(ij,l+1).le.0.
         enddo
      enddo

c Pas de pentes en bas et en haut
         do ij=1,ip1jmp1
            zqw(ij,2)=q(ij,1)
            zqw(ij,llm)=q(ij,llm)
            extremum(ij,1)=.true.
            extremum(ij,llm)=.true.
         enddo

c   calcul des valeurs max et min acceptees aux interfaces
      do l=2,llm
         do ij=1,ip1jmp1
            zqmax(ij,l)=max(q(ij,l-1),q(ij,l))
            zqmin(ij,l)=min(q(ij,l-1),q(ij,l))
         enddo
      enddo

      do l=2,llm
         do ij=1,ip1jmp1
            zqw(ij,l)=min(max(zqmin(ij,l),zqw(ij,l)),zqmax(ij,l))
         enddo
      enddo

      do l=2,llm-1
         do ij=1,ip1jmp1
            if(extremum(ij,l)) then
               qh(ij,l)=q(ij,l)
               qb(ij,l)=q(ij,l)
            else
               qh(ij,l)=zqw(ij,l+1)
               qb(ij,l)=zqw(ij,l)
            endif
         enddo
      enddo
c     do l=2,llm-1
c        do ij=1,ip1jmp1
c           if(extremum(ij,l)) then
c              if (.not.extremum(ij,l-1)) qh(ij,l-1)=q(ij,l)
c              if (.not.extremum(ij,l+1)) qb(ij,l+1)=q(ij,l)
c           endif
c        enddo
c     enddo

      do ij=1,ip1jmp1
         qb(ij,1)=q(ij,1)
         qh(ij,1)=q(ij,1)
         qb(ij,llm)=q(ij,llm)
         qh(ij,llm)=q(ij,llm)
      enddo

      endif

      RETURN
      END

      SUBROUTINE advnx(q,qg,qd,masse,u_m,mode)
c
c     Auteur : F. Hourdin
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
      use conf_gcm_m
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      integer mode
      real masse(ip1jmp1,llm)
      real u_m( ip1jmp1,llm )
      real q(ip1jmp1,llm),qd(ip1jmp1,llm),qg(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,j,ij,l,indu(ip1jmp1),niju,iju,ijq
      integer n0,nl(llm)
c
      real new_m,zu_m,zdq,zz
      real zsigg(ip1jmp1,llm),zsigd(ip1jmp1,llm),zsig
      real u_mq(ip1jmp1,llm)

      real zm,zq,zsigm,zsigp,zqm,zqp,zu

      logical ladvplus(ip1jmp1,llm)

      real prec
      save prec

      data prec/1.e-15/

      do l=1,llm
            do ij=iip2,ip1jm
               zdq=qd(ij,l)-qg(ij,l)
               if(abs(zdq).gt.prec) then
                  zsigd(ij,l)=(q(ij,l)-qg(ij,l))/zdq
                  zsigg(ij,l)=1.-zsigd(ij,l)
               else
                  zsigd(ij,l)=0.5
                  zsigg(ij,l)=0.5
                  qd(ij,l)=q(ij,l)
                  qg(ij,l)=q(ij,l)
               endif
            enddo
       enddo

c   calcul de la pente maximum dans la maille en valeur absolue

       do l=1,llm
       do ij=iip2,ip1jm-1
          if (u_m(ij,l).ge.0.) then
             zsigp=zsigd(ij,l)
             zsigm=zsigg(ij,l)
             zqp=qd(ij,l)
             zqm=qg(ij,l)
             zm=masse(ij,l)
             zq=q(ij,l)
          else
             zsigm=zsigd(ij+1,l)
             zsigp=zsigg(ij+1,l)
             zqm=qd(ij+1,l)
             zqp=qg(ij+1,l)
             zm=masse(ij+1,l)
             zq=q(ij+1,l)
          endif
          zu=abs(u_m(ij,l))
          ladvplus(ij,l)=zu.gt.zm
          zsig=zu/zm
          if(zsig.eq.0.) zsigp=0.1
          if (mode.eq.1) then
             if (zsig.le.zsigp) then
                 u_mq(ij,l)=u_m(ij,l)*zqp
             else if (mode.eq.1) then
                 u_mq(ij,l)=
     s           sign(zm,u_m(ij,l))*(zsigp*zqp+(zsig-zsigp)*zqm)
             endif 
          else
             if (zsig.le.zsigp) then
                 u_mq(ij,l)=u_m(ij,l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
             else
                zz=0.5*(zsig-zsigp)/zsigm
                u_mq(ij,l)=sign(zm,u_m(ij,l))*( 0.5*(zq+zqp)*zsigp
     s          +(zsig-zsigp)*(zq+zz*(zqm-zq)) )
             endif
          endif
      enddo
      enddo

      do l=1,llm
       do ij=iip1+iip1,ip1jm,iip1
          u_mq(ij,l)=u_mq(ij-iim,l)
          ladvplus(ij,l)=ladvplus(ij-iim,l)
       enddo
      enddo

c=================================================================
C   SCHEMA SEMI-LAGRAGIEN EN X DANS LES REGIONS POLAIRES
c=================================================================
c   tris des regions a traiter
      n0=0
      do l=1,llm
         nl(l)=0
         do ij=iip2,ip1jm
            if(ladvplus(ij,l)) then
               nl(l)=nl(l)+1
               u_mq(ij,l)=0.
            endif
         enddo
         n0=n0+nl(l)
      enddo

      if(n0.gt.1) then
      IF (prt_level > 9) print *,
     & 'Nombre de points pour lesquels on advect plus que le'
     &       ,'contenu de la maille : ',n0

         do l=1,llm
            if(nl(l).gt.0) then
               iju=0
c   indicage des mailles concernees par le traitement special
               do ij=iip2,ip1jm
                  if(ladvplus(ij,l).and.mod(ij,iip1).ne.0) then
                     iju=iju+1
                     indu(iju)=ij
                  endif
               enddo
               niju=iju

c  traitement des mailles
               do iju=1,niju
                  ij=indu(iju)
                  j=(ij-1)/iip1+1
                  zu_m=u_m(ij,l)
                  u_mq(ij,l)=0.
                  if(zu_m.gt.0.) then
                     ijq=ij
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)+q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m-masse(ijq,l)
                        i=mod(i-2+iim,iim)+1
                        ijq=(j-1)*iip1+i
                     enddo
c   MODIFS SPECIFIQUES DU SCHEMA
c   ajout de la maille non completement advectee
             zsig=zu_m/masse(ijq,l)
             if(zsig.le.zsigd(ijq,l)) then
                u_mq(ij,l)=u_mq(ij,l)+zu_m*(qd(ijq,l)
     s          -0.5*zsig/zsigd(ijq,l)*(qd(ijq,l)-q(ijq,l)))
             else
c               u_mq(ij,l)=u_mq(ij,l)+zu_m*q(ijq,l)
c         goto 8888
                zz=0.5*(zsig-zsigd(ijq,l))/zsigg(ijq,l)
                if(.not.(zz.gt.0..and.zz.le.0.5)) then
                     print *,'probleme2 au point ij=',ij,
     s               '  l=',l
                     print *,'zz=',zz
                     stop
                endif
                u_mq(ij,l)=u_mq(ij,l)+masse(ijq,l)*(
     s          0.5*(q(ijq,l)+qd(ijq,l))*zsigd(ijq,l)
     s        +(zsig-zsigd(ijq,l))*(q(ijq,l)+zz*(qg(ijq,l)-q(ijq,l))) )
             endif
                  else
                     ijq=ij+1
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(-zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)-q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m+masse(ijq,l)
                        i=mod(i,iim)+1
                        ijq=(j-1)*iip1+i
                     enddo
c   ajout de la maille non completement advectee
c 2eme MODIF SPECIFIQUE
             zsig=-zu_m/masse(ij+1,l)
             if(zsig.le.zsigg(ijq,l)) then
                u_mq(ij,l)=u_mq(ij,l)+zu_m*(qg(ijq,l)
     s          -0.5*zsig/zsigg(ijq,l)*(qg(ijq,l)-q(ijq,l)))
             else
c               u_mq(ij,l)=u_mq(ij,l)+zu_m*q(ijq,l)
c           goto 9999
                zz=0.5*(zsig-zsigg(ijq,l))/zsigd(ijq,l)
                if(.not.(zz.gt.0..and.zz.le.0.5)) then
                     print *,'probleme22 au point ij=',ij
     s               ,'  l=',l
                     print *,'zz=',zz
                     stop
                endif
                u_mq(ij,l)=u_mq(ij,l)-masse(ijq,l)*(
     s          0.5*(q(ijq,l)+qg(ijq,l))*zsigg(ijq,l)
     s          +(zsig-zsigg(ijq,l))*
     s           (q(ijq,l)+zz*(qd(ijq,l)-q(ijq,l))) )
             endif
c   fin de la modif
                  endif
               enddo
            endif
         enddo
      endif  ! n0.gt.0 

c   bouclage en latitude
      do l=1,llm
        do ij=iip1+iip1,ip1jm,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        enddo
      enddo

c=================================================================
c   CALCUL DE LA CONVERGENCE DES FLUX
c=================================================================

      do l=1,llm
         do ij=iip2+1,ip1jm
            new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+
     &      u_mq(ij-1,l)-u_mq(ij,l))
     &      /new_m
            masse(ij,l)=new_m
         enddo
c   Modif Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         do ij=iip1+iip1,ip1jm,iip1
            q(ij-iim,l)=q(ij,l)
            masse(ij-iim,l)=masse(ij,l)
         enddo
      enddo

      RETURN
      END
      SUBROUTINE advny(q,qs,qn,masse,v_m)
c
c     Auteur : F. Hourdin
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
      use comgeom
      use conf_gcm_m
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real masse(ip1jmp1,llm)
      real v_m( ip1jm,llm )
      real q(ip1jmp1,llm),qn(ip1jmp1,llm),qs(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l
c
      real new_m,zdq,zz
      real zsigs(ip1jmp1),zsign(ip1jmp1),zsig
      real v_mq(ip1jm,llm)
      real convpn,convps,convmpn,convmps,massen,masses
      real zm,zq,zsigm,zsigp,zqm,zqp
      real ssum
      real prec
      save prec

      data prec/1.e-15/
      do l=1,llm
            do ij=1,ip1jmp1
               zdq=qn(ij,l)-qs(ij,l)
               if(abs(zdq).gt.prec) then
                  zsign(ij)=(q(ij,l)-qs(ij,l))/zdq
                  zsigs(ij)=1.-zsign(ij)
               else
                  zsign(ij)=0.5
                  zsigs(ij)=0.5
               endif
            enddo

c   calcul de la pente maximum dans la maille en valeur absolue

       do ij=1,ip1jm
          if (v_m(ij,l).ge.0.) then
             zsigp=zsign(ij+iip1)
             zsigm=zsigs(ij+iip1)
             zqp=qn(ij+iip1,l)
             zqm=qs(ij+iip1,l)
             zm=masse(ij+iip1,l)
             zq=q(ij+iip1,l)
          else
             zsigm=zsign(ij)
             zsigp=zsigs(ij)
             zqm=qn(ij,l)
             zqp=qs(ij,l)
             zm=masse(ij,l)
             zq=q(ij,l)
          endif
          zsig=abs(v_m(ij,l))/zm
          if(zsig.eq.0.) zsigp=0.1
          if (zsig.le.zsigp) then
              v_mq(ij,l)=v_m(ij,l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
          else
              zz=0.5*(zsig-zsigp)/zsigm
              v_mq(ij,l)=sign(zm,v_m(ij,l))*( 0.5*(zq+zqp)*zsigp 
     s        +(zsig-zsigp)*(zq+zz*(zqm-zq)) )
          endif
       enddo
      enddo

      do l=1,llm
         do ij=iip2,ip1jm
            new_m=masse(ij,l)
     &      +v_m(ij,l)-v_m(ij-iip1,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+v_mq(ij,l)-v_mq(ij-iip1,l))
     &         /new_m
            masse(ij,l)=new_m
         enddo
c.-. ancienne version
         convpn=SSUM(iim,v_mq(1,l),1)
         convmpn=ssum(iim,v_m(1,l),1)
         massen=ssum(iim,masse(1,l),1)
         new_m=massen+convmpn
         q(1,l)=(q(1,l)*massen+convpn)/new_m
         do ij = 1,iip1
            q(ij,l)=q(1,l)
            masse(ij,l)=new_m*aire(ij)/apoln
         enddo

         convps=-SSUM(iim,v_mq(ip1jm-iim,l),1)
         convmps=-ssum(iim,v_m(ip1jm-iim,l),1)
         masses=ssum(iim,masse(ip1jm+1,l),1)
         new_m=masses+convmps
         q(ip1jm+1,l)=(q(ip1jm+1,l)*masses+convps)/new_m
         do ij = ip1jm+1,ip1jmp1
            q(ij,l)=q(ip1jm+1,l)
            masse(ij,l)=new_m*aire(ij)/apols
         enddo
      enddo

      RETURN
      END
      SUBROUTINE advnz(q,qh,qb,masse,w_m)
c
c     Auteurs:   F.Hourdin
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c     b designe le bas et h le haut
c     il y a une correspondance entre le b en z et le d en x
c    ********************************************************************
c
c
c   --------------------------------------------------------------------
      use dimens_m
      use paramet_m
      use comgeom
      use conf_gcm_m
      IMPLICIT NONE
c
c
c
c   Arguments:
c   ----------
      real masse(ip1jmp1,llm)
      real w_m( ip1jmp1,llm+1)
      real q(ip1jmp1,llm),qb(ip1jmp1,llm),qh(ip1jmp1,llm)

c
c      Local 
c   ---------
c
      INTEGER ij,l
c
      real new_m,zdq,zz
      real zsigh(ip1jmp1,llm),zsigb(ip1jmp1,llm),zsig
      real w_mq(ip1jmp1,llm+1)
      real zm,zq,zsigm,zsigp,zqm,zqp
      real prec
      save prec

      data prec/1.e-13/

      do l=1,llm
            do ij=1,ip1jmp1
               zdq=qb(ij,l)-qh(ij,l)
               if(abs(zdq).gt.prec) then
                  zsigb(ij,l)=(q(ij,l)-qh(ij,l))/zdq
                  zsigh(ij,l)=1.-zsigb(ij,l)
                  zsigb(ij,l)=min(max(zsigb(ij,l),0.),1.)
               else
                  zsigb(ij,l)=0.5
                  zsigh(ij,l)=0.5
               endif
            enddo
       enddo

c   calcul de la pente maximum dans la maille en valeur absolue
       do l=2,llm
       do ij=1,ip1jmp1
          if (w_m(ij,l).ge.0.) then
             zsigp=zsigb(ij,l)
             zsigm=zsigh(ij,l)
             zqp=qb(ij,l)
             zqm=qh(ij,l)
             zm=masse(ij,l)
             zq=q(ij,l)
          else
             zsigm=zsigb(ij,l-1)
             zsigp=zsigh(ij,l-1)
             zqm=qb(ij,l-1)
             zqp=qh(ij,l-1)
             zm=masse(ij,l-1)
             zq=q(ij,l-1)
          endif
          zsig=abs(w_m(ij,l))/zm
          if(zsig.eq.0.) zsigp=0.1
          if (zsig.le.zsigp) then
              w_mq(ij,l)=w_m(ij,l)*(zqp-0.5*zsig/zsigp*(zqp-zq))
          else
              zz=0.5*(zsig-zsigp)/zsigm
              w_mq(ij,l)=sign(zm,w_m(ij,l))*( 0.5*(zq+zqp)*zsigp
     s        +(zsig-zsigp)*(zq+zz*(zqm-zq)) )
          endif
      enddo
      enddo

       do ij=1,ip1jmp1
          w_mq(ij,llm+1)=0.
          w_mq(ij,1)=0.
       enddo

      do l=1,llm
         do ij=1,ip1jmp1
            new_m=masse(ij,l)+w_m(ij,l+1)-w_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+w_mq(ij,l+1)-w_mq(ij,l))
     &         /new_m
            masse(ij,l)=new_m
         enddo
      enddo

      END
