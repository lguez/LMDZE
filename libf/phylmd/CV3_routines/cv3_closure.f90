
      SUBROUTINE cv3_closure(nloc,ncum,nd,icb,inb &
                            ,pbase,p,ph,tv,buoy &
                            ,sig,w0,cape,m)
            use cvparam3
            use cvthermo
      implicit none

!===================================================================
! ---  CLOSURE OF CONVECT3
!
! vectorization: S. Bony
!===================================================================


! input:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real pbase(nloc)
      real p(nloc,nd), ph(nloc,nd+1)
      real tv(nloc,nd), buoy(nloc,nd)

! input/output:
      real sig(nloc,nd), w0(nloc,nd)

! output:
      real cape(nloc)
      real m(nloc,nd)

! local variables:
      integer i, j, k, icbmax
      real deltap, fac, w, amu
      real dtmin(nloc,nd), sigold(nloc,nd)


! -------------------------------------------------------
! -- Initialization
! -------------------------------------------------------

      do k=1,nl
       do i=1,ncum
        m(i,k)=0.0
       enddo
      enddo

! -------------------------------------------------------
! -- Reset sig(i) and w0(i) for i>inb and i<icb
! -------------------------------------------------------

! update sig and w0 above LNB:

      do 100 k=1,nl-1
       do 110 i=1,ncum
        if ((inb(i).lt.(nl-1)).and.(k.ge.(inb(i)+1)))then
         sig(i,k)=beta*sig(i,k) &
                  +2.*alpha*buoy(i,inb(i))*ABS(buoy(i,inb(i)))
         sig(i,k)=AMAX1(sig(i,k),0.0)
         w0(i,k)=beta*w0(i,k)
        endif
 110   continue
 100  continue

! compute icbmax:

      icbmax=2
      do 200 i=1,ncum
        icbmax=MAX(icbmax,icb(i))
 200  continue

! update sig and w0 below cloud base:

      do 300 k=1,icbmax
       do 310 i=1,ncum
        if (k.le.icb(i))then
         sig(i,k)=beta*sig(i,k)-2.*alpha*buoy(i,icb(i))*buoy(i,icb(i))
         sig(i,k)=amax1(sig(i,k),0.0)
         w0(i,k)=beta*w0(i,k)
        endif
310    continue
300    continue

!!      if(inb.lt.(nl-1))then
!!         do 85 i=inb+1,nl-1
!!            sig(i)=beta*sig(i)+2.*alpha*buoy(inb)*
!!     1              abs(buoy(inb))
!!            sig(i)=amax1(sig(i),0.0)
!!            w0(i)=beta*w0(i)
!!   85    continue
!!      end if

!!      do 87 i=1,icb
!!         sig(i)=beta*sig(i)-2.*alpha*buoy(icb)*buoy(icb)
!!         sig(i)=amax1(sig(i),0.0)
!!         w0(i)=beta*w0(i)
!!   87 continue

! -------------------------------------------------------------
! -- Reset fractional areas of updrafts and w0 at initial time
! -- and after 10 time steps of no convection
! -------------------------------------------------------------

      do 400 k=1,nl-1
       do 410 i=1,ncum
        if (sig(i,nd).lt.1.5.or.sig(i,nd).gt.12.0)then
         sig(i,k)=0.0
         w0(i,k)=0.0
        endif
 410   continue
 400  continue

! -------------------------------------------------------------
! -- Calculate convective available potential energy (cape),
! -- vertical velocity (w), fractional area covered by
! -- undilute updraft (sig), and updraft mass flux (m)
! -------------------------------------------------------------

      do 500 i=1,ncum
       cape(i)=0.0
 500  continue

! compute dtmin (minimum buoyancy between ICB and given level k):

      do i=1,ncum
       do k=1,nl
         dtmin(i,k)=100.0
       enddo
      enddo

      do 550 i=1,ncum
       do 560 k=1,nl
         do 570 j=minorig,nl
          if ( (k.ge.(icb(i)+1)).and.(k.le.inb(i)).and. &
               (j.ge.icb(i)).and.(j.le.(k-1)) )then
           dtmin(i,k)=AMIN1(dtmin(i,k),buoy(i,j))
          endif
 570     continue
 560   continue
 550  continue

! the interval on which cape is computed starts at pbase :

      do 600 k=1,nl
       do 610 i=1,ncum

        if ((k.ge.(icb(i)+1)).and.(k.le.inb(i))) then

         deltap = MIN(pbase(i),ph(i,k-1))-MIN(pbase(i),ph(i,k))
         cape(i)=cape(i)+rrd*buoy(i,k-1)*deltap/p(i,k-1)
         cape(i)=AMAX1(0.0,cape(i))
         sigold(i,k)=sig(i,k)

!         dtmin(i,k)=100.0
!         do 97 j=icb(i),k-1 ! mauvaise vectorisation
!          dtmin(i,k)=AMIN1(dtmin(i,k),buoy(i,j))
!  97     continue

         sig(i,k)=beta*sig(i,k)+alpha*dtmin(i,k)*ABS(dtmin(i,k))
         sig(i,k)=amax1(sig(i,k),0.0)
         sig(i,k)=amin1(sig(i,k),0.01)
         fac=AMIN1(((dtcrit-dtmin(i,k))/dtcrit),1.0)
         w=(1.-beta)*fac*SQRT(cape(i))+beta*w0(i,k)
         amu=0.5*(sig(i,k)+sigold(i,k))*w
         m(i,k)=amu*0.007*p(i,k)*(ph(i,k)-ph(i,k+1))/tv(i,k)
         w0(i,k)=w
        endif

 610   continue
 600  continue

      do 700 i=1,ncum
       w0(i,icb(i))=0.5*w0(i,icb(i)+1)
       m(i,icb(i))=0.5*m(i,icb(i)+1) &
                   *(ph(i,icb(i))-ph(i,icb(i)+1)) &
                   /(ph(i,icb(i)+1)-ph(i,icb(i)+2))
       sig(i,icb(i))=sig(i,icb(i)+1)
       sig(i,icb(i)-1)=sig(i,icb(i))
 700  continue


!!      cape=0.0
!!      do 98 i=icb+1,inb
!!         deltap = min(pbase,ph(i-1))-min(pbase,ph(i))
!!         cape=cape+rrd*buoy(i-1)*deltap/p(i-1)
!!         dcape=rrd*buoy(i-1)*deltap/p(i-1)
!!         dlnp=deltap/p(i-1)
!!         cape=amax1(0.0,cape)
!!         sigold=sig(i)

!!         dtmin=100.0
!!         do 97 j=icb,i-1
!!            dtmin=amin1(dtmin,buoy(j))
!!   97    continue

!!         sig(i)=beta*sig(i)+alpha*dtmin*abs(dtmin)
!!         sig(i)=amax1(sig(i),0.0)
!!         sig(i)=amin1(sig(i),0.01)
!!         fac=amin1(((dtcrit-dtmin)/dtcrit),1.0)
!!         w=(1.-beta)*fac*sqrt(cape)+beta*w0(i)
!!         amu=0.5*(sig(i)+sigold)*w
!!         m(i)=amu*0.007*p(i)*(ph(i)-ph(i+1))/tv(i)
!!         w0(i)=w
!!   98 continue
!!      w0(icb)=0.5*w0(icb+1)
!!      m(icb)=0.5*m(icb+1)*(ph(icb)-ph(icb+1))/(ph(icb+1)-ph(icb+2))
!!      sig(icb)=sig(icb+1)
!!      sig(icb-1)=sig(icb)

       return
       end
