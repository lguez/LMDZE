c
c $Header: /home/cvsroot/LMDZ4/libf/phylmd/cvltr.F,v 1.1 2005/04/15 12:36:17 lmdzadmin Exp $
c
      SUBROUTINE cvltr(pdtime,da, phi, mp,paprs,pplay,x,upd,dnd,dx)
      use dimens_m
      use dimphy
      use YOMCST
      IMPLICIT NONE 
c=====================================================================
c Objet : convection des traceurs / KE
c Auteurs: M-A Filiberti and J-Y Grandpeix
c=====================================================================
c
      include "YOECUMF.h" 
c
      REAL pdtime
      REAL, intent(in):: paprs(klon,klev+1) ! pression aux 1/2 couches (bas en haut)
      REAL pplay(klon,klev)  ! pression pour le milieu de chaque couche
      REAL x(klon,klev)        ! q de traceur (bas en haut) 
      REAL dx(klon,klev)     ! tendance de traceur  (bas en haut)
      real da(klon,klev),phi(klon,klev,klev),mp(klon,klev)
      REAL upd(klon,klev)      ! saturated updraft mass flux
      REAL dnd(klon,klev)      ! saturated downdraft mass flux
c
c--variables locales      
      real zed(klon,klev),zmd(klon,klev,klev)
      real za(klon,klev,klev)
      real zmfd(klon,klev),zmfa(klon,klev)
      real zmfp(klon,klev),zmfu(klon,klev)
      integer i,k,j 
c test conservation
c      real conserv
c =========================================
c calcul des tendances liees au downdraft
c =========================================
      zed(:,:)=0.
      zmfd(:,:)=0.
      zmfa(:,:)=0.
      zmfu(:,:)=0.
      zmfp(:,:)=0.
      zmd(:,:,:)=0.
      za(:,:,:)=0.
c entrainement
      do k=1,klev-1
        do i=1,klon
          zed(i,k)=max(0.,mp(i,k)-mp(i,k+1))
        end do
      end do
c
c calcul de la matrice d echange
c matrice de distribution de la masse entrainee en k
c
      do k=1,klev
        do i=1,klon
          zmd(i,k,k)=zed(i,k)
        end do
      end do
      do k=2,klev
        do j=k-1,1,-1
          do i=1,klon
          if(mp(i,j+1).ne.0) then
          zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1))
          endif
          end do
        end do
      end do
      do k=1,klev
        do j=1,klev-1
          do i=1,klon
          za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
          end do
        end do
      end do
c
c rajout du terme lie a l ascendance induite
c
        do j=2,klev
         do i=1,klon
          za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
         end do
        end do
C
c tendances
c            
      do k=1,klev
        do j=1,klev
          do i=1,klon
          zmfd(i,j)=zmfd(i,j)+za(i,j,k)*(x(i,k)-x(i,j))
          end do
        end do
      end do
c
c =========================================
c calcul des tendances liees aux flux satures
c =========================================
      do j=1,klev
        do i=1,klon
          zmfa(i,j)=da(i,j)*(x(i,1)-x(i,j))
        end do
      end do
      do k=1,klev
        do j=1,klev
          do i=1,klon
          zmfp(i,j)=zmfp(i,j)+phi(i,j,k)*(x(i,k)-x(i,j))
          end do
        end do
      end do
      do j=1,klev-1
        do i=1,klon
          zmfu(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x(i,j+1)-x(i,j))
        end do
      end do
      do j=2,klev
        do i=1,klon
          zmfu(i,j)=zmfu(i,j)
     .             +min(0.,upd(i,j)+dnd(i,j))*(x(i,j)-x(i,j-1))
        end do
      end do

c =========================================
c--calcul final des tendances
c =========================================
      do k=1, klev
        do i=1, klon
          dx(i,k)=(zmfd(i,k)+zmfu(i,k)
     .      +zmfa(i,k)+zmfp(i,k))*pdtime
     .      *RG/(paprs(i,k)-paprs(i,k+1))
c          print*,'dx',k,dx(i,k)
        enddo
      enddo
c
c test de conservation du traceur
c      conserv=0.
c      do k=1, klev
c        do i=1, klon
c         conserv=conserv+dx(i,k)*
c     .     (paprs(i,k)-paprs(i,k+1))/RG
C
c        enddo
c      enddo
c      print *,'conserv',conserv
     
      return
      end
