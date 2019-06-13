
      subroutine dqthermcell2(ngrid,nlay,ptimestep,fm,entr,masse,frac &
          ,q,dq,qa)
      use dimensions
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
      real entr(ngrid,nlay),frac(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),wqd(klon,klev+1)
      real qe(klon,klev),zf,zf2

      integer ig,k

!   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

!   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
         qe(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt. &
               1.e-5*masse(ig,k)) then
               zf=0.5*(frac(ig,k)+frac(ig,k+1))
               zf2=1./(1.-zf)
               qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+zf2*entr(ig,k)*q(ig,k)) &
               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2)
               qe(ig,k)=(q(ig,k)-zf*qa(ig,k))*zf2
            else
               qa(ig,k)=q(ig,k)
               qe(ig,k)=q(ig,k)
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
!             wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
            wqd(ig,k)=fm(ig,k)*qe(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(detr(ig,k)*qa(ig,k)-entr(ig,k)*qe(ig,k) &
                     -wqd(ig,k)+wqd(ig,k+1)) &
                     /masse(ig,k)
         enddo
      enddo

      return
      end
