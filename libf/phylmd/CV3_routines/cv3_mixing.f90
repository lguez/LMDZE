
      SUBROUTINE cv3_mixing(nloc,ncum,nd,na,ntra,icb,nk,inb &
                          ,ph,t,rr,rs,u,v,tra,h,lv,qnk &
                          ,hp,tv,tvp,ep,clw,m,sig &
         ,ment,qent,uent,vent, nent, sij,elij,ments,qents,traent)
            use cvparam3
            use cvthermo
      implicit none

!---------------------------------------------------------------------
! a faire:
!     - changer rr(il,1) -> qnk(il)
!   - vectorisation de la partie normalisation des flux (do 789...)
!---------------------------------------------------------------------


! inputs:
      integer ncum, nd, na, ntra, nloc
      integer icb(nloc), inb(nloc), nk(nloc)
      real sig(nloc,nd)
      real qnk(nloc)
      real ph(nloc,nd+1)
      real t(nloc,nd), rr(nloc,nd), rs(nloc,nd)
      real u(nloc,nd), v(nloc,nd)
      real tra(nloc,nd,ntra) ! input of convect3
      real lv(nloc,na), h(nloc,na), hp(nloc,na)
      real tv(nloc,na), tvp(nloc,na), ep(nloc,na), clw(nloc,na)
      real m(nloc,na)        ! input of convect3

! outputs:
      real ment(nloc,na,na), qent(nloc,na,na)
      real uent(nloc,na,na), vent(nloc,na,na)
      real sij(nloc,na,na), elij(nloc,na,na)
      real traent(nloc,nd,nd,ntra)
      real ments(nloc,nd,nd), qents(nloc,nd,nd)
      real sigij(nloc,nd,nd)
      integer nent(nloc,nd)

! local variables:
      integer i, j, k, il, im, jm
      integer num1, num2
      real rti, bf2, anum, denom, dei, altem, cwat, stemp, qp
      real alt, smid, sjmin, sjmax, delp, delm
      real asij(nloc), smax(nloc), scrit(nloc)
      real asum(nloc,nd),bsum(nloc,nd),csum(nloc,nd)
      real wgh
      real zm(nloc,na)
      logical lwork(nloc)

!=====================================================================
! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
!=====================================================================

! ori        do 360 i=1,ncum*nlp
        do 361 j=1,nl
        do 360 i=1,ncum
          nent(i,j)=0
! in convect3, m is computed in cv3_closure
! ori          m(i,1)=0.0
 360    continue
 361    continue

! ori      do 400 k=1,nlp
! ori       do 390 j=1,nlp
      do 400 j=1,nl
       do 390 k=1,nl
          do 385 i=1,ncum
            qent(i,k,j)=rr(i,j)
            uent(i,k,j)=u(i,j)
            vent(i,k,j)=v(i,j)
            elij(i,k,j)=0.0
!ym            ment(i,k,j)=0.0
!ym            sij(i,k,j)=0.0
 385      continue
 390    continue
 400  continue

!ym
      ment(1:ncum,1:nd,1:nd)=0.0
      sij(1:ncum,1:nd,1:nd)=0.0

!      do k=1,ntra
!       do j=1,nd  ! instead nlp
!        do i=1,nd ! instead nlp
!         do il=1,ncum
!            traent(il,i,j,k)=tra(il,j,k)
!         enddo
!        enddo
!       enddo
!      enddo
      zm(:,:)=0.

!=====================================================================
! --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
! --- FRACTION (sij)
!=====================================================================

      do 750 i=minorig+1, nl

       do 710 j=minorig,nl
        do 700 il=1,ncum
         if( (i.ge.icb(il)).and.(i.le.inb(il)).and. &
            (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then

          rti=rr(il,1)-ep(il,i)*clw(il,i)
          bf2=1.+lv(il,j)*lv(il,j)*rs(il,j)/(rrv*t(il,j)*t(il,j)*cpd)
          anum=h(il,j)-hp(il,i)+(cpv-cpd)*t(il,j)*(rti-rr(il,j))
          denom=h(il,i)-hp(il,i)+(cpd-cpv)*(rr(il,i)-rti)*t(il,j)
          dei=denom
          if(abs(dei).lt.0.01)dei=0.01
          sij(il,i,j)=anum/dei
          sij(il,i,i)=1.0
          altem=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti-rs(il,j)
          altem=altem/bf2
          cwat=clw(il,j)*(1.-ep(il,j))
          stemp=sij(il,i,j)
          if((stemp.lt.0.0.or.stemp.gt.1.0.or.altem.gt.cwat) &
                       .and.j.gt.i)then
           anum=anum-lv(il,j)*(rti-rs(il,j)-cwat*bf2)
           denom=denom+lv(il,j)*(rr(il,i)-rti)
           if(abs(denom).lt.0.01)denom=0.01
           sij(il,i,j)=anum/denom
           altem=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti-rs(il,j)
           altem=altem-(bf2-1.)*cwat
          end if
         if(sij(il,i,j).gt.0.0.and.sij(il,i,j).lt.0.95)then
          qent(il,i,j)=sij(il,i,j)*rr(il,i)+(1.-sij(il,i,j))*rti
          uent(il,i,j)=sij(il,i,j)*u(il,i)+(1.-sij(il,i,j))*u(il,nk(il))
          vent(il,i,j)=sij(il,i,j)*v(il,i)+(1.-sij(il,i,j))*v(il,nk(il))
!!!!      do k=1,ntra
!!!!      traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
!!!!     :      +(1.-sij(il,i,j))*tra(il,nk(il),k)
!!!!      end do
          elij(il,i,j)=altem
          elij(il,i,j)=amax1(0.0,elij(il,i,j))
          ment(il,i,j)=m(il,i)/(1.-sij(il,i,j))
          nent(il,i)=nent(il,i)+1
         end if
         sij(il,i,j)=amax1(0.0,sij(il,i,j))
         sij(il,i,j)=amin1(1.0,sij(il,i,j))
         endif ! new
 700   continue
 710  continue

!       do k=1,ntra
!        do j=minorig,nl
!         do il=1,ncum
!          if( (i.ge.icb(il)).and.(i.le.inb(il)).and.
!     :       (j.ge.(icb(il)-1)).and.(j.le.inb(il)))then
!            traent(il,i,j,k)=sij(il,i,j)*tra(il,i,k)
!     :            +(1.-sij(il,i,j))*tra(il,nk(il),k)
!          endif
!         enddo
!        enddo
!       enddo

!
!   ***   if no air can entrain at level i assume that updraft detrains  ***
!   ***   at that level and calculate detrained air flux and properties  ***
!

!@      do 170 i=icb(il),inb(il)

      do 740 il=1,ncum
      if ((i.ge.icb(il)).and.(i.le.inb(il)).and.(nent(il,i).eq.0)) then
!@      if(nent(il,i).eq.0)then
      ment(il,i,i)=m(il,i)
      qent(il,i,i)=rr(il,nk(il))-ep(il,i)*clw(il,i)
      uent(il,i,i)=u(il,nk(il))
      vent(il,i,i)=v(il,nk(il))
      elij(il,i,i)=clw(il,i)
!MAF      sij(il,i,i)=1.0
      sij(il,i,i)=0.0
      end if
 740  continue
 750  continue

!      do j=1,ntra
!       do i=minorig+1,nl
!        do il=1,ncum
!         if (i.ge.icb(il) .and. i.le.inb(il) .and. nent(il,i).eq.0) then
!          traent(il,i,i,j)=tra(il,nk(il),j)
!         endif
!        enddo
!       enddo
!      enddo

      do 100 j=minorig,nl
      do 101 i=minorig,nl
      do 102 il=1,ncum
      if ((j.ge.(icb(il)-1)).and.(j.le.inb(il)) &
          .and.(i.ge.icb(il)).and.(i.le.inb(il)))then
       sigij(il,i,j)=sij(il,i,j)
      endif
 102  continue
 101  continue
 100  continue
!@      enddo

!@170   continue

!=====================================================================
!   ---  NORMALIZE ENTRAINED AIR MASS FLUXES
!   ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
!=====================================================================

!ym      call zilch(asum,ncum*nd)
!ym      call zilch(bsum,ncum*nd)
!ym      call zilch(csum,ncum*nd)
      call zilch(asum,nloc*nd)
      call zilch(csum,nloc*nd)
      call zilch(csum,nloc*nd)

      do il=1,ncum
       lwork(il) = .FALSE.
      enddo

      DO 789 i=minorig+1,nl

      num1=0
      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) ) num1=num1+1
      enddo
      if (num1.le.0) goto 789


      do 781 il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) ) then
        lwork(il)=(nent(il,i).ne.0)
        qp=rr(il,1)-ep(il,i)*clw(il,i)
        anum=h(il,i)-hp(il,i)-lv(il,i)*(qp-rs(il,i)) &
                 +(cpv-cpd)*t(il,i)*(qp-rr(il,i))
        denom=h(il,i)-hp(il,i)+lv(il,i)*(rr(il,i)-qp) &
                 +(cpd-cpv)*t(il,i)*(rr(il,i)-qp)
        if(abs(denom).lt.0.01)denom=0.01
        scrit(il)=anum/denom
        alt=qp-rs(il,i)+scrit(il)*(rr(il,i)-qp)
        if(scrit(il).le.0.0.or.alt.le.0.0)scrit(il)=1.0
        smax(il)=0.0
        asij(il)=0.0
       endif
781   continue

      do 175 j=nl,minorig,-1

      num2=0
      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) .and. &
            j.ge.(icb(il)-1) .and. j.le.inb(il)  &
            .and. lwork(il) ) num2=num2+1
      enddo
      if (num2.le.0) goto 175

      do 782 il=1,ncum
      if ( i.ge.icb(il) .and. i.le.inb(il) .and. &
            j.ge.(icb(il)-1) .and. j.le.inb(il)  &
            .and. lwork(il) ) then

       if(sij(il,i,j).gt.1.0e-16.and.sij(il,i,j).lt.0.95)then
        wgh=1.0
        if(j.gt.i)then
         sjmax=amax1(sij(il,i,j+1),smax(il))
         sjmax=amin1(sjmax,scrit(il))
         smax(il)=amax1(sij(il,i,j),smax(il))
         sjmin=amax1(sij(il,i,j-1),smax(il))
         sjmin=amin1(sjmin,scrit(il))
         if(sij(il,i,j).lt.(smax(il)-1.0e-16))wgh=0.0
         smid=amin1(sij(il,i,j),scrit(il))
        else
         sjmax=amax1(sij(il,i,j+1),scrit(il))
         smid=amax1(sij(il,i,j),scrit(il))
         sjmin=0.0
         if(j.gt.1)sjmin=sij(il,i,j-1)
         sjmin=amax1(sjmin,scrit(il))
        endif
        delp=abs(sjmax-smid)
        delm=abs(sjmin-smid)
        asij(il)=asij(il)+wgh*(delp+delm)
        ment(il,i,j)=ment(il,i,j)*(delp+delm)*wgh
       endif
      endif
782   continue

175   continue

      do il=1,ncum
       if (i.ge.icb(il).and.i.le.inb(il).and.lwork(il)) then
        asij(il)=amax1(1.0e-16,asij(il))
        asij(il)=1.0/asij(il)
        asum(il,i)=0.0
        bsum(il,i)=0.0
        csum(il,i)=0.0
       endif
      enddo

      do 180 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il) &
         .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         ment(il,i,j)=ment(il,i,j)*asij(il)
        endif
       enddo
180   continue

      do 190 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il) &
         .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         asum(il,i)=asum(il,i)+ment(il,i,j)
         ment(il,i,j)=ment(il,i,j)*sig(il,j)
         bsum(il,i)=bsum(il,i)+ment(il,i,j)
        endif
       enddo
190   continue

      do il=1,ncum
       if (i.ge.icb(il).and.i.le.inb(il).and.lwork(il)) then
        bsum(il,i)=amax1(bsum(il,i),1.0e-16)
        bsum(il,i)=1.0/bsum(il,i)
       endif
      enddo

      do 195 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il) &
         .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         ment(il,i,j)=ment(il,i,j)*asum(il,i)*bsum(il,i)
        endif
       enddo
195   continue

      do 197 j=minorig,nl
       do il=1,ncum
        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il) &
         .and. j.ge.(icb(il)-1) .and. j.le.inb(il) ) then
         csum(il,i)=csum(il,i)+ment(il,i,j)
        endif
       enddo
197   continue

      do il=1,ncum
       if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il) &
           .and. csum(il,i).lt.m(il,i) ) then
        nent(il,i)=0
        ment(il,i,i)=m(il,i)
        qent(il,i,i)=rr(il,1)-ep(il,i)*clw(il,i)
        uent(il,i,i)=u(il,nk(il))
        vent(il,i,i)=v(il,nk(il))
        elij(il,i,i)=clw(il,i)
!MAF        sij(il,i,i)=1.0
        sij(il,i,i)=0.0
       endif
      enddo ! il

!      do j=1,ntra
!       do il=1,ncum
!        if ( i.ge.icb(il) .and. i.le.inb(il) .and. lwork(il)
!     :     .and. csum(il,i).lt.m(il,i) ) then
!         traent(il,i,i,j)=tra(il,nk(il),j)
!        endif
!       enddo
!      enddo
789   continue
!
! MAF: renormalisation de MENT
      do jm=1,nd
        do im=1,nd
          do il=1,ncum
          zm(il,im)=zm(il,im)+(1.-sij(il,im,jm))*ment(il,im,jm)
         end do
        end do
      end do
!
      do jm=1,nd
        do im=1,nd
          do il=1,ncum
          if(zm(il,im).ne.0.) then
          ment(il,im,jm)=ment(il,im,jm)*m(il,im)/zm(il,im)
          endif
         end do
       end do
      end do
!
      do jm=1,nd
       do im=1,nd
        do 999 il=1,ncum
         qents(il,im,jm)=qent(il,im,jm)
         ments(il,im,jm)=ment(il,im,jm)
999     continue
       enddo
      enddo

      return
      end
