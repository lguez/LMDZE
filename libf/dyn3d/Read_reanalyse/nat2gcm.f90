
!===========================================================================
      subroutine nat2gcm(u,v,t,rh,pk,ucov,vcov,teta,q)
!===========================================================================

      use dimens_m
      use paramet_m
      use comconst
      use comvert
      use comgeom
      use q_sat_m, only: q_sat
      use guide_m
      implicit none


      real u(iip1,jjp1,llm),v(iip1,jjm,llm)
      real t(iip1,jjp1,llm),pk(iip1,jjp1,llm),rh(iip1,jjp1,llm)
      real ps(iip1,jjp1)

      real ucov(iip1,jjp1,llm),vcov(iip1,jjm,llm)
      real teta(iip1,jjp1,llm),q(iip1,jjp1,llm)

      real pres(iip1,jjp1,llm),qsat(iip1,jjp1,llm)

      real unskap

      integer i,j,l


      print*,'Entree dans nat2gcm'
!    ucov(:,:,:)=0.
!    do l=1,llm
!       ucov(:,2:jjm,l)=u(:,2:jjm,l)*cu_2d(:,2:jjm)
!    enddo
!    ucov(iip1,:,:)=ucov(1,:,:)

!    teta(:,:,:)=t(:,:,:)*cpp/pk(:,:,:)
!    teta(iip1,:,:)=teta(1,:,:)

!   calcul de ucov et de la temperature potentielle
      do l=1,llm
         do j=1,jjp1
            do i=1,iim
               ucov(i,j,l)=u(i,j,l)*cu_2d(i,j)
               teta(i,j,l)=t(i,j,l)*cpp/pk(i,j,l)
            enddo
            ucov(iip1,j,l)=ucov(1,j,l)
            teta(iip1,j,l)=teta(1,j,l)
         enddo
         do i=1,iip1
            ucov(i,1,l)=0.
            ucov(i,jjp1,l)=0.
            teta(i,1,l)=teta(1,1,l)
            teta(i,jjp1,l)=teta(1,jjp1,l)
         enddo
      enddo

!   calcul de ucov
      do l=1,llm
         do j=1,jjm
            do i=1,iim
               vcov(i,j,l)=v(i,j,l)*cv_2d(i,j)
            enddo
            vcov(iip1,j,l)=vcov(1,j,l)
         enddo
      enddo

!     call dump2d(iip1,jjp1,teta,'TETA EN BAS   ')
!     call dump2d(iip1,jjp1,teta(1,1,llm),'TETA EN HAUT   ')

!  Humidite relative -> specifique
!  -------------------------------
      if (1.eq.0) then
!   FINALEMENT ON GUIDE EN HUMIDITE RELATIVE
      print*,'calcul de unskap'
      unskap   = 1./ kappa
      print*,'calcul de pres'
      pres(:,:,:)=preff*(pk(:,:,:)/cpp)**unskap
      print*,'calcul de qsat'
      qsat = q_sat(t, pres)
      print*,'calcul de q'
!   ATTENTION : humidites relatives en %
      rh(:,:,:)=max(rh(:,:,:)*0.01,1.e-6)
      q(:,:,:)=qsat(:,:,:)*rh(:,:,:)
      print*,'calcul de q OK'

      call dump2d(iip1,jjp1,pres,'PRESSION PREMIERE COUCHE   ')
      call dump2d(iip1,jjp1,q,'HUMIDITE SPECIFIQUE COUCHE 1   ')
      endif


      return
      end
