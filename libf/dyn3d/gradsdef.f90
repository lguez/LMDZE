module gradsdef

  ! From dyn3d/gradsdef.h,v 1.1.1.1 2004/05/19 12:53:07

  implicit none

  integer, parameter:: nfmx=10,imx=200,jmx=150,lmx=100,nvarmx=1000

  real xd(imx,nfmx),yd(jmx,nfmx),zd(lmx,nfmx),dtime(nfmx)

  integer imd(imx),jmd(jmx),lmd(lmx)
  integer iid(imx),jid(jmx)
  integer ifd(imx),jfd(jmx)
  integer unit(nfmx),irec(nfmx),itime(nfmx),nld(nvarmx,nfmx)

  integer nvar(nfmx),ivar(nfmx)
  logical firsttime(nfmx)

  character(len=10) var(nvarmx,nfmx),fichier(nfmx)
  character(len=40) title(nfmx),tvar(nvarmx,nfmx)

end module gradsdef
