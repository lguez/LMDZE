
      SUBROUTINE cv_unsat(nloc,ncum,nd,inb,t,q,qs,gz,u,v,p,ph &
                        ,h,lv,ep,sigp,clw,m,ment,elij &
                        ,iflag,mp,qp,up,vp,wt,water,evap)
            use cvthermo
            use cvparam
      implicit none



! inputs:
      integer, intent(in):: ncum, nd, nloc
      integer inb(nloc)
      real t(nloc,nd), q(nloc,nd), qs(nloc,nd)
      real gz(nloc,nd), u(nloc,nd), v(nloc,nd)
      real p(nloc,nd), ph(nloc,nd+1), h(nloc,nd)
      real lv(nloc,nd), ep(nloc,nd), sigp(nloc,nd), clw(nloc,nd)
      real m(nloc,nd), ment(nloc,nd,nd), elij(nloc,nd,nd)

! outputs:
      integer iflag(nloc) ! also an input
      real mp(nloc,nd), qp(nloc,nd), up(nloc,nd), vp(nloc,nd)
      real water(nloc,nd), evap(nloc,nd), wt(nloc,nd)

! local variables:
      integer i,j,k,ij,num1
      integer jtt(nloc)
      real awat, coeff, qsm, afac, sigt, b6, c6, revap
      real dhdp, fac, qstm, rat
      real wdtrain(nloc)
      logical lwork(nloc)

!=====================================================================
! --- PRECIPITATING DOWNDRAFT CALCULATION
!=====================================================================
!
! Initializations:
!
         do i = 1, ncum
         do k = 1, nl+1
          wt(i,k) = omtsnow
          mp(i,k) = 0.0
          evap(i,k) = 0.0
          water(i,k) = 0.0
         enddo
         enddo

         do 420 i=1,ncum
          qp(i,1)=q(i,1)
          up(i,1)=u(i,1)
          vp(i,1)=v(i,1)
 420     continue

         do 440 k=2,nl+1
         do 430 i=1,ncum
          qp(i,k)=q(i,k-1)
          up(i,k)=u(i,k-1)
          vp(i,k)=v(i,k-1)
 430     continue
 440     continue


!   ***  Check whether ep(inb)=0, if so, skip precipitating    ***
!   ***             downdraft calculation                      ***
!
!
!   ***  Integrate liquid water equation to find condensed water   ***
!   ***                and condensed water flux                    ***
!
!
      do 890 i=1,ncum
        jtt(i)=2
        if(ep(i,inb(i)).le.0.0001)iflag(i)=2
        if(iflag(i).eq.0)then
          lwork(i)=.true.
        else
          lwork(i)=.false.
        endif
 890  continue
!
!    ***                    Begin downdraft loop                    ***
!
!
        call zilch(wdtrain,ncum)
        do 899 i=nl+1,1,-1
!
          num1=0
          do 879 ij=1,ncum
            if((i.le.inb(ij)).and.lwork(ij))num1=num1+1
 879      continue
          if(num1.le.0)go to 899
!
!
!    ***        Calculate detrained precipitation             ***
!
          do 891 ij=1,ncum
            if((i.le.inb(ij)).and.(lwork(ij)))then
            wdtrain(ij)=g*ep(ij,i)*m(ij,i)*clw(ij,i)
            endif
 891      continue
!
          if(i.gt.1)then
            do 893 j=1,i-1
              do 892 ij=1,ncum
                if((i.le.inb(ij)).and.(lwork(ij)))then
                  awat=elij(ij,j,i)-(1.-ep(ij,i))*clw(ij,i)
                  awat=max(0.0,awat)
                  wdtrain(ij)=wdtrain(ij)+g*awat*ment(ij,j,i)
                endif
 892          continue
 893      continue
          endif
!
!    ***    Find rain water and evaporation using provisional   ***
!    ***              estimates of qp(i)and qp(i-1)             ***
!
!
!  ***  Value of terminal velocity and coeffecient of evaporation for snow   ***
!
          do 894 ij=1,ncum
            if((i.le.inb(ij)).and.(lwork(ij)))then
            coeff=coeffs
            wt(ij,i)=omtsnow
!
!  ***  Value of terminal velocity and coeffecient of evaporation for rain   ***
!
            if(t(ij,i).gt.273.0)then
              coeff=coeffr
              wt(ij,i)=omtrain
            endif
            qsm=0.5*(q(ij,i)+qp(ij,i+1))
            afac=coeff*ph(ij,i)*(qs(ij,i)-qsm) &
             /(1.0e4+2.0e3*ph(ij,i)*qs(ij,i))
            afac=max(afac,0.0)
            sigt=sigp(ij,i)
            sigt=max(0.0,sigt)
            sigt=min(1.0,sigt)
            b6=100.*(ph(ij,i)-ph(ij,i+1))*sigt*afac/wt(ij,i)
            c6=(water(ij,i+1)*wt(ij,i+1)+wdtrain(ij)/sigd)/wt(ij,i)
            revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
            evap(ij,i)=sigt*afac*revap
            water(ij,i)=revap*revap
!
!    ***  Calculate precipitating downdraft mass flux under     ***
!    ***              hydrostatic approximation                 ***
!
            if(i.gt.1)then
              dhdp=(h(ij,i)-h(ij,i-1))/(p(ij,i-1)-p(ij,i))
              dhdp=max(dhdp,10.0)
              mp(ij,i)=100.*ginv*lv(ij,i)*sigd*evap(ij,i)/dhdp
              mp(ij,i)=max(mp(ij,i),0.0)
!
!   ***   Add small amount of inertia to downdraft              ***
!
              fac=20.0/(ph(ij,i-1)-ph(ij,i))
              mp(ij,i)=(fac*mp(ij,i+1)+mp(ij,i))/(1.+fac)
!
!    ***      Force mp to decrease linearly to zero                 ***
!    ***      between about 950 mb and the surface                  ***
!
              if(p(ij,i).gt.(0.949*p(ij,1)))then
                 jtt(ij)=max(jtt(ij),i)
                 mp(ij,i)=mp(ij,jtt(ij))*(p(ij,1)-p(ij,i)) &
                 /(p(ij,1)-p(ij,jtt(ij)))
              endif
            endif
!
!    ***       Find mixing ratio of precipitating downdraft     ***
!
            if(i.ne.inb(ij))then
              if(i.eq.1)then
                qstm=qs(ij,1)
              else
                qstm=qs(ij,i-1)
              endif
              if(mp(ij,i).gt.mp(ij,i+1))then
                 rat=mp(ij,i+1)/mp(ij,i)
                 qp(ij,i)=qp(ij,i+1)*rat+q(ij,i)*(1.0-rat)+100.*ginv* &
                   sigd*(ph(ij,i)-ph(ij,i+1))*(evap(ij,i)/mp(ij,i))
                 up(ij,i)=up(ij,i+1)*rat+u(ij,i)*(1.-rat)
                 vp(ij,i)=vp(ij,i+1)*rat+v(ij,i)*(1.-rat)
               else
                 if(mp(ij,i+1).gt.0.0)then
                   qp(ij,i)=(gz(ij,i+1)-gz(ij,i) &
                     +qp(ij,i+1)*(lv(ij,i+1)+t(ij,i+1) &
                     *(cl-cpd))+cpd*(t(ij,i+1)-t(ij,i))) &
                     /(lv(ij,i)+t(ij,i)*(cl-cpd))
                   up(ij,i)=up(ij,i+1)
                   vp(ij,i)=vp(ij,i+1)
                 endif
              endif
              qp(ij,i)=min(qp(ij,i),qstm)
              qp(ij,i)=max(qp(ij,i),0.0)
            endif
            endif
 894      continue
 899    continue
!
        return
        end
