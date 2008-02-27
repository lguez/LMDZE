!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/yamada.F,v 1.1 2004/06/22 11:45:36 lmdzadmin Exp $
!
      SUBROUTINE yamada(ngrid,dt,g,rconst,plev,temp
     s   ,zlev,zlay,u,v,teta,cd,q2,km,kn,ustar
     s   ,l_mix)
      use dimens_m
      use dimphy
      IMPLICIT NONE
c.......................................................................
c.......................................................................
c
c dt : pas de temps
c g  : g
c zlev : altitude a chaque niveau (interface inferieure de la couche
c        de meme indice)
c zlay : altitude au centre de chaque couche
c u,v : vitesse au centre de chaque couche
c       (en entree : la valeur au debut du pas de temps)
c teta : temperature potentielle au centre de chaque couche
c        (en entree : la valeur au debut du pas de temps)
c cd : cdrag
c      (en entree : la valeur au debut du pas de temps)
c q2 : $q^2$ au bas de chaque couche
c      (en entree : la valeur au debut du pas de temps)
c      (en sortie : la valeur a la fin du pas de temps)
c km : diffusivite turbulente de quantite de mouvement (au bas de chaque
c      couche)
c      (en sortie : la valeur a la fin du pas de temps)
c kn : diffusivite turbulente des scalaires (au bas de chaque couche)
c      (en sortie : la valeur a la fin du pas de temps)
c
c.......................................................................
      REAL dt,g,rconst
      real plev(klon,klev+1),temp(klon,klev)
      real ustar(klon),snstable
      REAL zlev(klon,klev+1)
      REAL zlay(klon,klev)
      REAL u(klon,klev)
      REAL v(klon,klev)
      REAL teta(klon,klev)
      REAL cd(klon)
      REAL q2(klon,klev+1)
      REAL km(klon,klev+1)
      REAL kn(klon,klev+1)
      integer l_mix,ngrid


      integer nlay,nlev
      PARAMETER (nlay=klev)
      PARAMETER (nlev=klev+1)

      logical first
      save first
      data first/.true./


      integer ig,k

      real ri,zrif,zalpha,zsm
      real rif(klon,klev+1),sm(klon,klev+1),alpha(klon,klev)

      real m2(klon,klev+1),dz(klon,klev+1),zq,n2(klon,klev+1)
      real l(klon,klev+1),l0(klon)

      real sq(klon),sqz(klon),zz(klon,klev+1)
      integer iter

      real ric,rifc,b1,kap
      save ric,rifc,b1,kap
      data ric,rifc,b1,kap/0.195,0.191,16.6,0.3/

      real frif,falpha,fsm

      frif(ri)=0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))
      falpha(ri)=1.318*(0.2231-ri)/(0.2341-ri)
      fsm(ri)=1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))

      if (0.eq.1.and.first) then
      do ig=1,1000
         ri=(ig-800.)/500.
         if (ri.lt.ric) then
            zrif=frif(ri)
         else
            zrif=rifc
         endif
         if(zrif.lt.0.16) then
            zalpha=falpha(zrif)
            zsm=fsm(zrif)
         else
            zalpha=1.12
            zsm=0.085
         endif
         print*,ri,rif,zalpha,zsm
      enddo
      first=.false.
      endif

c  Correction d'un bug sauvage a verifier.
c      do k=2,nlev
      do k=2,nlay
                                                          do ig=1,ngrid
         dz(ig,k)=zlay(ig,k)-zlay(ig,k-1)
         m2(ig,k)=((u(ig,k)-u(ig,k-1))**2+(v(ig,k)-v(ig,k-1))**2)
     s             /(dz(ig,k)*dz(ig,k))
         n2(ig,k)=g*2.*(teta(ig,k)-teta(ig,k-1))
     s            /(teta(ig,k-1)+teta(ig,k))  /dz(ig,k)
         ri=n2(ig,k)/max(m2(ig,k),1.e-10)
         if (ri.lt.ric) then
            rif(ig,k)=frif(ri)
         else
            rif(ig,k)=rifc
         endif
         if(rif(ig,k).lt.0.16) then
            alpha(ig,k)=falpha(rif(ig,k))
            sm(ig,k)=fsm(rif(ig,k))
         else
            alpha(ig,k)=1.12
            sm(ig,k)=0.085
         endif
         zz(ig,k)=b1*m2(ig,k)*(1.-rif(ig,k))*sm(ig,k)
                                                          enddo
      enddo

c iterration pour determiner la longueur de melange

                                                          do ig=1,ngrid
      l0(ig)=100.
                                                          enddo
      do k=2,klev-1
                                                          do ig=1,ngrid
        l(ig,k)=l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
                                                          enddo
      enddo

      do iter=1,10
                                                          do ig=1,ngrid
         sq(ig)=1.e-10
         sqz(ig)=1.e-10
                                                          enddo
         do k=2,klev-1
                                                          do ig=1,ngrid
           q2(ig,k)=l(ig,k)**2*zz(ig,k)
           l(ig,k)=min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
     s     ,0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.e-10)))
           zq=sqrt(q2(ig,k))
           sqz(ig)=sqz(ig)+zq*zlev(ig,k)*(zlay(ig,k)-zlay(ig,k-1))
           sq(ig)=sq(ig)+zq*(zlay(ig,k)-zlay(ig,k-1))
                                                          enddo
         enddo
                                                          do ig=1,ngrid
         l0(ig)=0.2*sqz(ig)/sq(ig)
                                                          enddo
c(abd 3 5 2)         print*,'ITER=',iter,'  L0=',l0

      enddo

      do k=2,klev
                                                          do ig=1,ngrid
         l(ig,k)=min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig,k)+l0(ig))
     s     ,0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.e-10)))
         q2(ig,k)=l(ig,k)**2*zz(ig,k)
         km(ig,k)=l(ig,k)*sqrt(q2(ig,k))*sm(ig,k)
         kn(ig,k)=km(ig,k)*alpha(ig,k)
                                                          enddo
      enddo

      return
      end
