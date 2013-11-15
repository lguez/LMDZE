!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/groupeun.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
      subroutine groupeun(jjmax,llmax,q)
      use dimens_m
      use paramet_m
      use comconst
      use comgeom
      implicit none


      integer jjmax,llmax
      real q(iip1,jjmax,llmax)

      integer ngroup
      parameter (ngroup=3)

      real airen,airecn,qn
      real aires,airecs,qs

      integer i,j,l,ig,j1,j2,i0,jd

Champs 3D
      jd=jjp1-jjmax
      do l=1,llm
      j1=1+jd
      j2=2
      do ig=1,ngroup
         do j=j1-jd,j2-jd
            do i0=1,iim,2**(ngroup-ig+1)
               airen=0.
               airecn=0.
               qn=0.
               aires=0.
               airecs=0.
               qs=0.
               do i=i0,i0+2**(ngroup-ig+1)-1
                  airen=airen+aire_2d(i,j)
                  aires=aires+aire_2d(i,jjp1-j+1)
                  qn=qn+q(i,j,l)
                  qs=qs+q(i,jjp1-j+1-jd,l)
               enddo
               airecn=0.
               airecs=0.
               do i=i0,i0+2**(ngroup-ig+1)-1
                  q(i,j,l)=qn*aire_2d(i,j)/airen
                  q(i,jjp1-j+1-jd,l)=qs*aire_2d(i,jjp1-j+1)/aires
               enddo
            enddo
            q(iip1,j,l)=q(1,j,l)
            q(iip1,jjp1-j+1-jd,l)=q(1,jjp1-j+1-jd,l)
         enddo
         j1=j2+1
         j2=j2+2**ig
      enddo
      enddo

      return
      end
