

!===========================================================================
      subroutine reanalyse2nat(nlevnc,psi &
         ,unc,vnc,tnc,qnc,psnc,pl,u,v,t,q &
         ,ps,masse,pk)
!===========================================================================

! -----------------------------------------------------------------
!   Inversion Nord/sud de la grille + interpollation sur les niveaux
!   verticaux du modele.
! -----------------------------------------------------------------

      use dimens_m
      use paramet_m
      use comconst
      use disvert_m
      use comgeom
      use exner_hyb_m, only: exner_hyb
      use conf_guide_m

      implicit none


      integer nlevnc
      real psi(iip1,jjp1)
      real u(iip1,jjp1,llm),v(iip1,jjm,llm)
      real t(iip1,jjp1,llm),ps(iip1,jjp1),q(iip1,jjp1,llm)

      real pl(nlevnc)
      real unc(iip1,jjp1,nlevnc),vnc(iip1,jjm,nlevnc)
      real tnc(iip1,jjp1,nlevnc),psnc(iip1,jjp1)
      real qnc(iip1,jjp1,nlevnc)

      real zu(iip1,jjp1,llm),zv(iip1,jjm,llm)
      real zt(iip1,jjp1,llm),zq(iip1,jjp1,llm)

      real pext(iip1,jjp1,llm)
      real pbarx(iip1,jjp1,llm),pbary(iip1,jjm,llm)
      real plunc(iip1,jjp1,llm),plvnc(iip1,jjm,llm)
      real plsnc(iip1,jjp1,llm)

      real p(iip1,jjp1,llmp1),pk(iip1,jjp1,llm),pks(iip1,jjp1)
      real pkf(iip1,jjp1,llm)
      real masse(iip1,jjp1,llm),pls(iip1,jjp1,llm)
      real prefkap,unskap


      integer i,j,l


! -----------------------------------------------------------------
!   calcul de la pression au milieu des couches.
! -----------------------------------------------------------------

      forall (l = 1: llm + 1) p(:, :, l) = ap(l) + bp(l) * psi
      call massdair(p,masse)
      CALL exner_hyb(psi,p,pks,pk,pkf)

!    ....  Calcul de pls , pression au milieu des couches ,en Pascals
      unskap=1./kappa
      prefkap =  preff  ** kappa
!     PRINT *,' Pref kappa unskap  ',preff,kappa,unskap
      DO l = 1, llm
       DO j=1,jjp1
        DO i =1, iip1
        pls(i,j,l) = preff * ( pk(i,j,l)/cpp) ** unskap
        ENDDO
       ENDDO
       ENDDO


! -----------------------------------------------------------------
!   calcul des pressions pour les grilles u et v
! -----------------------------------------------------------------

      do l=1,llm
      do j=1,jjp1
         do i=1,iip1
            pext(i,j,l)=pls(i,j,l)*aire_2d(i,j)
         enddo
      enddo
      enddo
      call massbar(pext, pbarx, pbary )
      do l=1,llm
      do j=1,jjp1
         do i=1,iip1
            plunc(i,jjp1+1-j,l)=pbarx(i,j,l)/aireu_2d(i,j)
            plsnc(i,jjp1+1-j,l)=pls(i,j,l)
         enddo
      enddo
      enddo
      do l=1,llm
      do j=1,jjm
         do i=1,iip1
            plvnc(i,jjm+1-j,l)=pbary(i,j,l)/airev_2d(i,j)
         enddo
      enddo
      enddo

! -----------------------------------------------------------------

      if (guide_P) then
      do j=1,jjp1
         do i=1,iim
            ps(i,j)=psnc(i,jjp1+1-j)
         enddo
         ps(iip1,j)=ps(1,j)
      enddo
      endif


! -----------------------------------------------------------------
      call pres2lev(unc,zu,nlevnc,llm,pl,plunc,iip1,jjp1)
      call pres2lev(vnc,zv,nlevnc,llm,pl,plvnc,iip1,jjm )
      call pres2lev(tnc,zt,nlevnc,llm,pl,plsnc,iip1,jjp1)
      call pres2lev(qnc,zq,nlevnc,llm,pl,plsnc,iip1,jjp1)

!     call dump2d(iip1,jjp1,ps,'PS    ')
!     call dump2d(iip1,jjp1,psu,'PS    ')
!     call dump2d(iip1,jjm,psv,'PS    ')
!  Inversion Nord/Sud
      do l=1,llm
         do j=1,jjp1
            do i=1,iim
               u(i,j,l)=zu(i,jjp1+1-j,l)
               t(i,j,l)=zt(i,jjp1+1-j,l)
               q(i,j,l)=zq(i,jjp1+1-j,l)
            enddo
            u(iip1,j,l)=u(1,j,l)
            t(iip1,j,l)=t(1,j,l)
            q(iip1,j,l)=q(1,j,l)
         enddo
      enddo

      do l=1,llm
         do j=1,jjm
            do i=1,iim
               v(i,j,l)=zv(i,jjm+1-j,l)
            enddo
            v(iip1,j,l)=v(1,j,l)
         enddo
      enddo

      return
      end
