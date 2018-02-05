      subroutine dvthermcell2(ngrid,nlay,ptimestep,fm,entr,masse &
          ,fraca,larga &
          ,u,v,du,dv,ua,va)
      use dimens_m
      use dimphy
      implicit none

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
!=======================================================================


      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real fraca(ngrid,nlay+1)
      real larga(ngrid)
      real entr(ngrid,nlay)
      real u(ngrid,nlay)
      real ua(ngrid,nlay)
      real du(ngrid,nlay)
      real v(ngrid,nlay)
      real va(ngrid,nlay)
      real dv(ngrid,nlay)

      real detr(klon,klev),zf,zf2
      real wvd(klon,klev+1),wud(klon,klev+1)
      real gamma0,gamma(klon,klev+1)
      real ue(klon,klev),ve(klon,klev)
      real dua,dva
      integer iter

      integer ig,k

!   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

!   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
         ue(ig,1)=u(ig,1)
         ve(ig,1)=v(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt. &
               1.e-5*masse(ig,k)) then
!   On itère sur la valeur du coeff de freinage.
!              gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
               gamma0=masse(ig,k) &
               *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) ) &
               *0.5/larga(ig) &
               *1.
!    s         *0.5
!              gamma0=0.
               zf=0.5*(fraca(ig,k)+fraca(ig,k+1))
               zf=0.
               zf2=1./(1.-zf)
!   la première fois on multiplie le coefficient de freinage
!   par le module du vent dans la couche en dessous.
               dua=ua(ig,k-1)-u(ig,k-1)
               dva=va(ig,k-1)-v(ig,k-1)
               do iter=1,5
!   On choisit une relaxation lineaire.
                  gamma(ig,k)=gamma0
!   On choisit une relaxation quadratique.
                  gamma(ig,k)=gamma0*sqrt(dua**2+dva**2)
                  ua(ig,k)=(fm(ig,k)*ua(ig,k-1) &
                     +(zf2*entr(ig,k)+gamma(ig,k))*u(ig,k)) &
                     /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2 &
                       +gamma(ig,k))
                  va(ig,k)=(fm(ig,k)*va(ig,k-1) &
                     +(zf2*entr(ig,k)+gamma(ig,k))*v(ig,k)) &
                     /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2 &
                       +gamma(ig,k))
!                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
                  dua=ua(ig,k)-u(ig,k)
                  dva=va(ig,k)-v(ig,k)
                  ue(ig,k)=(u(ig,k)-zf*ua(ig,k))*zf2
                  ve(ig,k)=(v(ig,k)-zf*va(ig,k))*zf2
               enddo
            else
               ua(ig,k)=u(ig,k)
               va(ig,k)=v(ig,k)
               ue(ig,k)=u(ig,k)
               ve(ig,k)=v(ig,k)
               gamma(ig,k)=0.
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
            wud(ig,k)=fm(ig,k)*ue(ig,k)
            wvd(ig,k)=fm(ig,k)*ve(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wud(ig,1)=0.
         wud(ig,nlay+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            du(ig,k)=((detr(ig,k)+gamma(ig,k))*ua(ig,k) &
                     -(entr(ig,k)+gamma(ig,k))*ue(ig,k) &
                     -wud(ig,k)+wud(ig,k+1)) &
                     /masse(ig,k)
            dv(ig,k)=((detr(ig,k)+gamma(ig,k))*va(ig,k) &
                     -(entr(ig,k)+gamma(ig,k))*ve(ig,k) &
                     -wvd(ig,k)+wvd(ig,k+1)) &
                     /masse(ig,k)
         enddo
      enddo

      return
      end
