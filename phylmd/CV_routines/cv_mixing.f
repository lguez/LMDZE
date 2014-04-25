
      SUBROUTINE cv_mixing(nloc,ncum,nd,icb,nk,inb,inb1 &
                          ,ph,t,q,qs,u,v,h,lv,qnk &
                          ,hp,tv,tvp,ep,clw,cbmf &
                          ,m,ment,qent,uent,vent,nent,sij,elij)
            use cvthermo
            use cvparam
      implicit none


! inputs:
      integer, intent(in):: ncum, nd, nloc
      integer icb(nloc), inb(nloc), inb1(nloc), nk(nloc)
      real cbmf(nloc), qnk(nloc)
      real ph(nloc,nd+1)
      real t(nloc,nd), q(nloc,nd), qs(nloc,nd), lv(nloc,nd)
      real u(nloc,nd), v(nloc,nd), h(nloc,nd), hp(nloc,nd)
      real tv(nloc,nd), tvp(nloc,nd), ep(nloc,nd), clw(nloc,nd)

! outputs:
      integer nent(nloc,nd)
      real m(nloc,nd), ment(nloc,nd,nd), qent(nloc,nd,nd)
      real uent(nloc,nd,nd), vent(nloc,nd,nd)
      real sij(nloc,nd,nd), elij(nloc,nd,nd)

! local variables:
      integer i, j, k, ij
      integer num1, num2
      real dbo, qti, bf2, anum, denom, dei, altem, cwat, stemp
      real alt, qp1, smid, sjmin, sjmax, delp, delm
      real work(nloc), asij(nloc), smin(nloc), scrit(nloc)
      real bsum(nloc,nd)
      logical lwork(nloc)

!=====================================================================
! --- INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS
!=====================================================================
!
        do 360 i=1,ncum*nlp
          nent(i,1)=0
          m(i,1)=0.0
 360    continue
!
      do 400 k=1,nlp
       do 390 j=1,nlp
          do 385 i=1,ncum
            qent(i,k,j)=q(i,j)
            uent(i,k,j)=u(i,j)
            vent(i,k,j)=v(i,j)
            elij(i,k,j)=0.0
            ment(i,k,j)=0.0
            sij(i,k,j)=0.0
 385      continue
 390    continue
 400  continue
!
!-------------------------------------------------------------------
! --- Calculate rates of mixing,  m(i)
!-------------------------------------------------------------------
!
      call zilch(work,ncum)
!
      do 670 j=minorig+1,nl
        do 660 i=1,ncum
          if((j.ge.(icb(i)+1)).and.(j.le.inb(i)))then
             k=min(j,inb1(i))
             dbo=abs(tv(i,k+1)-tvp(i,k+1)-tv(i,k-1)+tvp(i,k-1)) &
             +entp*0.04*(ph(i,k)-ph(i,k+1))
             work(i)=work(i)+dbo
             m(i,j)=cbmf(i)*dbo
          endif
 660    continue
 670  continue
      do 690 k=minorig+1,nl
        do 680 i=1,ncum
          if((k.ge.(icb(i)+1)).and.(k.le.inb(i)))then
            m(i,k)=m(i,k)/work(i)
          endif
 680    continue
 690  continue
!
!
!=====================================================================
! --- CALCULATE ENTRAINED AIR MASS FLUX (ment), TOTAL WATER MIXING
! --- RATIO (QENT), TOTAL CONDENSED WATER (elij), AND MIXING
! --- FRACTION (sij)
!=====================================================================
!
!
       do 750 i=minorig+1,nl
         do 710 j=minorig+1,nl
           do 700 ij=1,ncum
             if((i.ge.(icb(ij)+1)).and.(j.ge.icb(ij)) &
               .and.(i.le.inb(ij)).and.(j.le.inb(ij)))then
               qti=qnk(ij)-ep(ij,i)*clw(ij,i)
               bf2=1.+lv(ij,j)*lv(ij,j)*qs(ij,j) &
               /(rrv*t(ij,j)*t(ij,j)*cpd)
               anum=h(ij,j)-hp(ij,i)+(cpv-cpd)*t(ij,j)*(qti-q(ij,j))
               denom=h(ij,i)-hp(ij,i)+(cpd-cpv)*(q(ij,i)-qti)*t(ij,j)
               dei=denom
               if(abs(dei).lt.0.01)dei=0.01
               sij(ij,i,j)=anum/dei
               sij(ij,i,i)=1.0
               altem=sij(ij,i,j)*q(ij,i)+(1.-sij(ij,i,j))*qti-qs(ij,j)
               altem=altem/bf2
               cwat=clw(ij,j)*(1.-ep(ij,j))
               stemp=sij(ij,i,j)
               if((stemp.lt.0.0.or.stemp.gt.1.0.or. &
                 altem.gt.cwat).and.j.gt.i)then
                 anum=anum-lv(ij,j)*(qti-qs(ij,j)-cwat*bf2)
                 denom=denom+lv(ij,j)*(q(ij,i)-qti)
                 if(abs(denom).lt.0.01)denom=0.01
                 sij(ij,i,j)=anum/denom
                 altem=sij(ij,i,j)*q(ij,i)+(1.-sij(ij,i,j))*qti-qs(ij,j)
                 altem=altem-(bf2-1.)*cwat
               endif
               if(sij(ij,i,j).gt.0.0.and.sij(ij,i,j).lt.0.9)then
                 qent(ij,i,j)=sij(ij,i,j)*q(ij,i) &
                              +(1.-sij(ij,i,j))*qti
                 uent(ij,i,j)=sij(ij,i,j)*u(ij,i) &
                              +(1.-sij(ij,i,j))*u(ij,nk(ij))
                 vent(ij,i,j)=sij(ij,i,j)*v(ij,i) &
                              +(1.-sij(ij,i,j))*v(ij,nk(ij))
                 elij(ij,i,j)=altem
                 elij(ij,i,j)=max(0.0,elij(ij,i,j))
                 ment(ij,i,j)=m(ij,i)/(1.-sij(ij,i,j))
                 nent(ij,i)=nent(ij,i)+1
               endif
             sij(ij,i,j)=max(0.0,sij(ij,i,j))
             sij(ij,i,j)=min(1.0,sij(ij,i,j))
             endif
  700      continue
  710    continue
!
!   ***   If no air can entrain at level i assume that updraft detrains  ***
!   ***   at that level and calculate detrained air flux and properties  ***
!
           do 740 ij=1,ncum
             if((i.ge.(icb(ij)+1)).and.(i.le.inb(ij)) &
             .and.(nent(ij,i).eq.0))then
               ment(ij,i,i)=m(ij,i)
               qent(ij,i,i)=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
               uent(ij,i,i)=u(ij,nk(ij))
               vent(ij,i,i)=v(ij,nk(ij))
               elij(ij,i,i)=clw(ij,i)
               sij(ij,i,i)=1.0
             endif
 740       continue
 750   continue
!
      do 770 i=1,ncum
        sij(i,inb(i),inb(i))=1.0
 770  continue
!
!=====================================================================
!   ---  NORMALIZE ENTRAINED AIR MASS FLUXES
!   ---  TO REPRESENT EQUAL PROBABILITIES OF MIXING
!=====================================================================
!
       call zilch(bsum,ncum*nlp)
       do 780 ij=1,ncum
         lwork(ij)=.false.
 780   continue
       do 789 i=minorig+1,nl
!
         num1=0
         do 779 ij=1,ncum
           if((i.ge.icb(ij)+1).and.(i.le.inb(ij)))num1=num1+1
 779     continue
         if(num1.le.0)go to 789
!
           do 781 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij)))then
                lwork(ij)=(nent(ij,i).ne.0)
                qp1=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
                anum=h(ij,i)-hp(ij,i)-lv(ij,i)*(qp1-qs(ij,i))
                denom=h(ij,i)-hp(ij,i)+lv(ij,i)*(q(ij,i)-qp1)
                if(abs(denom).lt.0.01)denom=0.01
                scrit(ij)=anum/denom
                alt=qp1-qs(ij,i)+scrit(ij)*(q(ij,i)-qp1)
                if(scrit(ij).lt.0.0.or.alt.lt.0.0)scrit(ij)=1.0
                asij(ij)=0.0
                smin(ij)=1.0
             endif
 781       continue
         do 783 j=minorig,nl
!
         num2=0
         do 778 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij)) &
             .and.(j.ge.icb(ij)).and.(j.le.inb(ij)) &
             .and.lwork(ij))num2=num2+1
 778     continue
         if(num2.le.0)go to 783
!
           do 782 ij=1,ncum
             if((i.ge.icb(ij)+1).and.(i.le.inb(ij)) &
             .and.(j.ge.icb(ij)).and.(j.le.inb(ij)).and.lwork(ij))then
                  if(sij(ij,i,j).gt.0.0.and.sij(ij,i,j).lt.0.9)then
                    if(j.gt.i)then
                      smid=min(sij(ij,i,j),scrit(ij))
                      sjmax=smid
                      sjmin=smid
                        if(smid.lt.smin(ij) &
                        .and.sij(ij,i,j+1).lt.smid)then
                          smin(ij)=smid
                          sjmax=min(sij(ij,i,j+1),sij(ij,i,j),scrit(ij))
                          sjmin=max(sij(ij,i,j-1),sij(ij,i,j))
                          sjmin=min(sjmin,scrit(ij))
                        endif
                    else
                      sjmax=max(sij(ij,i,j+1),scrit(ij))
                      smid=max(sij(ij,i,j),scrit(ij))
                      sjmin=0.0
                      if(j.gt.1)sjmin=sij(ij,i,j-1)
                      sjmin=max(sjmin,scrit(ij))
                    endif
                    delp=abs(sjmax-smid)
                    delm=abs(sjmin-smid)
                    asij(ij)=asij(ij)+(delp+delm) &
                                 *(ph(ij,j)-ph(ij,j+1))
                    ment(ij,i,j)=ment(ij,i,j)*(delp+delm) &
                                 *(ph(ij,j)-ph(ij,j+1))
                  endif
              endif
  782    continue
  783    continue
            do 784 ij=1,ncum
            if((i.ge.icb(ij)+1).and.(i.le.inb(ij)).and.lwork(ij))then
               asij(ij)=max(1.0e-21,asij(ij))
               asij(ij)=1.0/asij(ij)
               bsum(ij,i)=0.0
            endif
 784        continue
            do 786 j=minorig,nl+1
              do 785 ij=1,ncum
                if((i.ge.icb(ij)+1).and.(i.le.inb(ij)) &
                .and.(j.ge.icb(ij)).and.(j.le.inb(ij)) &
                .and.lwork(ij))then
                   ment(ij,i,j)=ment(ij,i,j)*asij(ij)
                   bsum(ij,i)=bsum(ij,i)+ment(ij,i,j)
                endif
 785     continue
 786     continue
             do 787 ij=1,ncum
               if((i.ge.icb(ij)+1).and.(i.le.inb(ij)) &
               .and.(bsum(ij,i).lt.1.0e-18).and.lwork(ij))then
                 nent(ij,i)=0
                 ment(ij,i,i)=m(ij,i)
                 qent(ij,i,i)=q(ij,nk(ij))-ep(ij,i)*clw(ij,i)
                 uent(ij,i,i)=u(ij,nk(ij))
                 vent(ij,i,i)=v(ij,nk(ij))
                 elij(ij,i,i)=clw(ij,i)
                 sij(ij,i,i)=1.0
               endif
  787        continue
  789  continue
!
       return
       end
