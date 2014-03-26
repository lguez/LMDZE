
      SUBROUTINE cv_compress( len,nloc,ncum,nd &
         ,iflag1,nk1,icb1 &
         ,cbmf1,plcl1,tnk1,qnk1,gznk1 &
         ,t1,q1,qs1,u1,v1,gz1 &
         ,h1,lv1,cpn1,p1,ph1,tv1,tp1,tvp1,clw1 &
         ,iflag,nk,icb &
         ,cbmf,plcl,tnk,qnk,gznk &
         ,t,q,qs,u,v,gz,h,lv,cpn,p,ph,tv,tp,tvp,clw  &
         ,dph          )
            use cvparam
      implicit none


! inputs:
      integer len,ncum,nd,nloc
      integer iflag1(len),nk1(len),icb1(len)
      real cbmf1(len),plcl1(len),tnk1(len),qnk1(len),gznk1(len)
      real, intent(in):: t1(len,nd)
      real, intent(in):: q1(len,nd),qs1(len,nd),u1(len,nd),v1(len,nd)
      real gz1(len,nd),h1(len,nd),lv1(len,nd),cpn1(len,nd)
      real p1(len,nd),ph1(len,nd+1),tv1(len,nd),tp1(len,nd)
      real tvp1(len,nd),clw1(len,nd)

! outputs:
      integer iflag(nloc),nk(nloc),icb(nloc)
      real cbmf(nloc),plcl(nloc),tnk(nloc),qnk(nloc),gznk(nloc)
      real t(nloc,nd),q(nloc,nd),qs(nloc,nd),u(nloc,nd),v(nloc,nd)
      real gz(nloc,nd),h(nloc,nd),lv(nloc,nd),cpn(nloc,nd)
      real p(nloc,nd),ph(nloc,nd+1),tv(nloc,nd),tp(nloc,nd)
      real tvp(nloc,nd),clw(nloc,nd)
      real dph(nloc,nd)

! local variables:
      integer i,k,nn


      do 110 k=1,nl+1
       nn=0
      do 100 i=1,len
      if(iflag1(i).eq.0)then
        nn=nn+1
        t(nn,k)=t1(i,k)
        q(nn,k)=q1(i,k)
        qs(nn,k)=qs1(i,k)
        u(nn,k)=u1(i,k)
        v(nn,k)=v1(i,k)
        gz(nn,k)=gz1(i,k)
        h(nn,k)=h1(i,k)
        lv(nn,k)=lv1(i,k)
        cpn(nn,k)=cpn1(i,k)
        p(nn,k)=p1(i,k)
        ph(nn,k)=ph1(i,k)
        tv(nn,k)=tv1(i,k)
        tp(nn,k)=tp1(i,k)
        tvp(nn,k)=tvp1(i,k)
        clw(nn,k)=clw1(i,k)
      endif
 100    continue
 110  continue

      if (nn.ne.ncum) then
         print*,'strange! nn not equal to ncum: ',nn,ncum
         stop
      endif

      nn=0
      do 150 i=1,len
      if(iflag1(i).eq.0)then
      nn=nn+1
      cbmf(nn)=cbmf1(i)
      plcl(nn)=plcl1(i)
      tnk(nn)=tnk1(i)
      qnk(nn)=qnk1(i)
      gznk(nn)=gznk1(i)
      nk(nn)=nk1(i)
      icb(nn)=icb1(i)
      iflag(nn)=iflag1(i)
      endif
 150  continue

      do 170 k=1,nl
       do 160 i=1,ncum
        dph(i,k)=ph(i,k)-ph(i,k+1)
 160   continue
 170  continue

      return
      end
