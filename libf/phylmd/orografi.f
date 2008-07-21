!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/orografi.F,v 1.4 2005/12/01 11:27:29 fairhead Exp $
!
      SUBROUTINE drag_noro (nlon,nlev,dtime,paprs,pplay,
     e                   pmea,pstd, psig, pgam, pthe,ppic,pval,
     e                   kgwd,kdx,ktest,
     e                   t, u, v,
     s                   pulow, pvlow, pustr, pvstr,
     s                   d_t, d_u, d_v)
c
      use dimens_m
      use dimphy
      use YOMCST
      IMPLICIT none
c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 19950201
c Objet: Frottement de la montagne Interface
c======================================================================
c Arguments:
c dtime---input-R- pas d'integration (s)
c paprs---input-R-pression pour chaque inter-couche (en Pa)
c pplay---input-R-pression pour le mileu de chaque couche (en Pa)
c t-------input-R-temperature (K)
c u-------input-R-vitesse horizontale (m/s)
c v-------input-R-vitesse horizontale (m/s)
c
c d_t-----output-R-increment de la temperature             
c d_u-----output-R-increment de la vitesse u
c d_v-----output-R-increment de la vitesse v
c======================================================================
c
c ARGUMENTS
c
      INTEGER nlon,nlev
      REAL, intent(in):: dtime
      REAL, intent(in):: paprs(klon,klev+1)
      REAL, intent(in):: pplay(klon,klev)
      REAL pmea(nlon),pstd(nlon),psig(nlon),pgam(nlon),pthe(nlon)
      REAL ppic(nlon),pval(nlon)
      REAL pulow(nlon),pvlow(nlon),pustr(nlon),pvstr(nlon)
      REAL t(nlon,nlev), u(nlon,nlev), v(nlon,nlev)
      REAL d_t(nlon,nlev), d_u(nlon,nlev), d_v(nlon,nlev)
c
      INTEGER i, k, kgwd, kdx(nlon), ktest(nlon)
c
c Variables locales:
c
      REAL zgeom(klon,klev)
      REAL pdtdt(klon,klev), pdudt(klon,klev), pdvdt(klon,klev)
      REAL pt(klon,klev), pu(klon,klev), pv(klon,klev)
      REAL papmf(klon,klev),papmh(klon,klev+1)
c
c initialiser les variables de sortie (pour securite)
c
      DO i = 1,klon
         pulow(i) = 0.0
         pvlow(i) = 0.0
         pustr(i) = 0.0
         pvstr(i) = 0.0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
         pdudt(i,k)=0.0
         pdvdt(i,k)=0.0
         pdtdt(i,k)=0.0
      ENDDO
      ENDDO
c
c preparer les variables d'entree (attention: l'ordre des niveaux 
c verticaux augmente du haut vers le bas)
c
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1) 
         pu(i,k) = u(i,klev-k+1)
         pv(i,k) = v(i,klev-k+1)
         papmf(i,k) = pplay(i,klev-k+1)
      ENDDO
      ENDDO
      DO k = 1, klev+1
      DO i = 1, klon
         papmh(i,k) = paprs(i,klev-k+2)
      ENDDO
      ENDDO
      DO i = 1, klon
         zgeom(i,klev) = RD * pt(i,klev)
     .                  * LOG(papmh(i,klev+1)/papmf(i,klev))
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         zgeom(i,k) = zgeom(i,k+1) + RD * (pt(i,k)+pt(i,k+1))/2.0
     .               * LOG(papmf(i,k+1)/papmf(i,k))
      ENDDO
      ENDDO
c
c appeler la routine principale
c
      CALL orodrag(klon,klev,kgwd,kdx,ktest,
     .            dtime,
     .            papmh, papmf, zgeom,
     .            pt, pu, pv,
     .            pmea, pstd, psig, pgam, pthe, ppic,pval,
     .            pulow,pvlow,
     .            pdudt,pdvdt,pdtdt)
C
      DO k = 1, klev
      DO i = 1, klon
         d_u(i,klev+1-k) = dtime*pdudt(i,k)
         d_v(i,klev+1-k) = dtime*pdvdt(i,k)
         d_t(i,klev+1-k) = dtime*pdtdt(i,k)
         pustr(i)        = pustr(i)
cIM BUG  .                +rg*pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))
     .                    +pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))/RG
         pvstr(i)        = pvstr(i)
cIM BUG  .                +rg*pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))
     .                    +pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))/RG
      ENDDO
      ENDDO
c
      RETURN
      END
      SUBROUTINE orodrag( nlon,nlev 
     i                 , kgwd, kdx, ktest
     r                 , ptsphy
     r                 , paphm1,papm1,pgeom1,ptm1,pum1,pvm1
     r                 , pmea, pstd, psig, pgamma, ptheta, ppic, pval
c outputs
     r                 , pulow,pvlow
     r                 , pvom,pvol,pte )

      use dimens_m
      use dimphy
      use YOMCST
            use yoegwd
      implicit none

c
c
c**** *gwdrag* - does the gravity wave parametrization.
c
c     purpose.
c     --------
c
c          this routine computes the physical tendencies of the
c     prognostic variables u,v  and t due to  vertical transports by
c     subgridscale orographically excited gravity waves
c
c**   interface.
c     ----------
c          called from *callpar*.
c
c          the routine takes its input from the long-term storage:
c          u,v,t and p at t-1.
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
c     ==== outputs ===
c
c        implicit arguments :   none
c        --------------------
c
c      implicit logical (l)
c
c     method.
c     -------
c
c     externals.
c     ----------
      integer ismin, ismax
      external ismin, ismax
c
c     reference.
c     ----------
c
c     author.
c     -------
c     m.miller + b.ritter   e.c.m.w.f.     15/06/86.
c
c     f.lott + m. miller    e.c.m.w.f.     22/11/94
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
c
      integer nlon, nlev, klevm1
      integer kgwd, jl, ilevp1, jk, ji
      real zdelp, ztemp, zforc, ztend
      real rover, zb, zc, zconb, zabsv
      real zzd1, ratio, zbet, zust,zvst, zdis
      real  pte(nlon,nlev),
     *      pvol(nlon,nlev),
     *      pvom(nlon,nlev),
     *      pulow(klon),
     *      pvlow(klon)
      real  pum1(nlon,nlev),
     *      pvm1(nlon,nlev),
     *      ptm1(nlon,nlev),
     *      pmea(nlon),pstd(nlon),psig(nlon),
     *      pgamma(nlon),ptheta(nlon),ppic(nlon),pval(nlon),
     *      pgeom1(nlon,nlev),
     *      papm1(nlon,nlev),
     *      paphm1(nlon,nlev+1)
c
      integer  kdx(nlon),ktest(nlon)
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------
      integer  isect(klon),
     *         icrit(klon),
     *         ikcrith(klon),
     *         ikenvh(klon),
     *         iknu(klon),
     *         iknu2(klon),
     *         ikcrit(klon),
     *         ikhlim(klon)
c
      real   ztau(klon,klev+1),
     $       ztauf(klon,klev+1),
     *       zstab(klon,klev+1),
     *       zvph(klon,klev+1),
     *       zrho(klon,klev+1),
     *       zri(klon,klev+1),
     *       zpsi(klon,klev+1),
     *       zzdep(klon,klev)
      real   zdudt(klon),
     *       zdvdt(klon),
     *       zdtdt(klon),
     *       zdedt(klon),
     *       zvidis(klon),
     *       znu(klon),
     *       zd1(klon),
     *       zd2(klon),
     *       zdmod(klon)
      real ztmst, zrtmst 
      real, intent(in):: ptsphy
c
c------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
 100  continue
c
c     ------------------------------------------------------------------
c
c*         1.1   computational constants
c                -----------------------
c
 110  continue
c
c     ztmst=twodt
c     if(nstep.eq.nstart) ztmst=0.5*twodt
      klevm1=klev-1
      ztmst=ptsphy
      zrtmst=1./ztmst
c     ------------------------------------------------------------------
c
 120  continue
c
c     ------------------------------------------------------------------
c
c*         1.3   check whether row contains point for printing
c                ---------------------------------------------
c
 130  continue
c
c     ------------------------------------------------------------------
c
c*         2.     precompute basic state variables.
c*                ---------- ----- ----- ----------
c*                define low level wind, project winds in plane of
c*                low level wind, determine sector in which to take
c*                the variance and set indicator for critical levels.
c
  200 continue
c
c
c
      call orosetup
     *     ( nlon, ktest 
     *     , ikcrit, ikcrith, icrit,  ikenvh,iknu,iknu2
     *     , paphm1, papm1 , pum1   , pvm1 , ptm1 , pgeom1, pstd
     *     , zrho  , zri   , zstab  , ztau , zvph , zpsi, zzdep
     *     , pulow, pvlow 
     *     , ptheta,pgamma,pmea,ppic,pval,znu  ,zd1,  zd2,  zdmod )
c
c
c
c***********************************************************
c
c
c*         3.      compute low level stresses using subcritical and
c*                 supercritical forms.computes anisotropy coefficient
c*                 as measure of orographic twodimensionality.
c
  300 continue
c
      call gwstress
     *    ( nlon  , nlev
     *    , ktest , icrit, ikenvh, iknu
     *    , zrho  , zstab, zvph  , pstd,  psig, pmea, ppic
     *    , ztau 
     *    , pgeom1,zdmod)
c
c
c*         4.      compute stress profile.
c*                 ------- ------ --------
c
  400 continue
c
c
      call gwprofil
     *       (  nlon , nlev
     *       , kgwd   , kdx , ktest
     *       , ikcrith, icrit
     *       , paphm1, zrho   , zstab ,  zvph
     *       , zri   , ztau   
     *       , zdmod , psig  , pstd)
c
c
c*         5.      compute tendencies.
c*                 -------------------
c
  500 continue
c
c  explicit solution at all levels for the gravity wave
c  implicit solution for the blocked levels

      do 510 jl=1,klon
      zvidis(jl)=0.0
      zdudt(jl)=0.0
      zdvdt(jl)=0.0
      zdtdt(jl)=0.0
  510 continue
c
      ilevp1=klev+1
c
c
      do 524 jk=1,klev
c
c
c     do 523 jl=1,kgwd
c     ji=kdx(jl)
c  Modif vectorisation 02/04/2004
      do 523 ji=1,klon
      if(ktest(ji).eq.1) then

      zdelp=paphm1(ji,jk+1)-paphm1(ji,jk)
      ztemp=-rg*(ztau(ji,jk+1)-ztau(ji,jk))/(zvph(ji,ilevp1)*zdelp)
      zdudt(ji)=(pulow(ji)*zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
      zdvdt(ji)=(pvlow(ji)*zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)
c
c controle des overshoots:
c
      zforc=sqrt(zdudt(ji)**2+zdvdt(ji)**2)+1.E-12
      ztend=sqrt(pum1(ji,jk)**2+pvm1(ji,jk)**2)/ztmst+1.E-12
      rover=0.25
      if(zforc.ge.rover*ztend)then
        zdudt(ji)=rover*ztend/zforc*zdudt(ji)
        zdvdt(ji)=rover*ztend/zforc*zdvdt(ji)
      endif
c
c fin du controle des overshoots
c
      if(jk.ge.ikenvh(ji)) then
         zb=1.0-0.18*pgamma(ji)-0.04*pgamma(ji)**2
         zc=0.48*pgamma(ji)+0.3*pgamma(ji)**2
         zconb=2.*ztmst*gkwake*psig(ji)/(4.*pstd(ji))
         zabsv=sqrt(pum1(ji,jk)**2+pvm1(ji,jk)**2)/2.
         zzd1=zb*cos(zpsi(ji,jk))**2+zc*sin(zpsi(ji,jk))**2
	     ratio=(cos(zpsi(ji,jk))**2+pgamma(ji)*sin(zpsi(ji,jk))**2)/
     *   (pgamma(ji)*cos(zpsi(ji,jk))**2+sin(zpsi(ji,jk))**2)
         zbet=max(0.,2.-1./ratio)*zconb*zzdep(ji,jk)*zzd1*zabsv
c
c simplement oppose au vent
c
         zdudt(ji)=-pum1(ji,jk)/ztmst
         zdvdt(ji)=-pvm1(ji,jk)/ztmst
c
c  projection dans la direction de l'axe principal de l'orographie
cmod     zdudt(ji)=-(pum1(ji,jk)*cos(ptheta(ji)*rpi/180.)
cmod *              +pvm1(ji,jk)*sin(ptheta(ji)*rpi/180.))
cmod *              *cos(ptheta(ji)*rpi/180.)/ztmst
cmod     zdvdt(ji)=-(pum1(ji,jk)*cos(ptheta(ji)*rpi/180.)
cmod *              +pvm1(ji,jk)*sin(ptheta(ji)*rpi/180.))
cmod *              *sin(ptheta(ji)*rpi/180.)/ztmst
         zdudt(ji)=zdudt(ji)*(zbet/(1.+zbet))
         zdvdt(ji)=zdvdt(ji)*(zbet/(1.+zbet))
      end if
      pvom(ji,jk)=zdudt(ji)
      pvol(ji,jk)=zdvdt(ji)
      zust=pum1(ji,jk)+ztmst*zdudt(ji)
      zvst=pvm1(ji,jk)+ztmst*zdvdt(ji)
      zdis=0.5*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
      zdedt(ji)=zdis/ztmst
      zvidis(ji)=zvidis(ji)+zdis*zdelp
      zdtdt(ji)=zdedt(ji)/rcpd
c     pte(ji,jk)=zdtdt(ji)
c
c  ENCORE UN TRUC POUR EVITER LES EXPLOSIONS
c
      pte(ji,jk)=0.0

      endif
  523 continue

  524 continue
c
c
      return
      end
      SUBROUTINE orosetup
     *         ( nlon   , ktest
     *         , kkcrit, kkcrith, kcrit
     *         , kkenvh, kknu  , kknu2
     *         , paphm1, papm1 , pum1   , pvm1 , ptm1  , pgeom1, pstd
     *         , prho  , pri   , pstab  , ptau , pvph  ,ppsi, pzdep
     *         , pulow , pvlow  
     *         , ptheta, pgamma, pmea, ppic, pval
     *         , pnu  ,  pd1  ,  pd2  ,pdmod  )
c
c**** *gwsetup*
c
c     purpose.
c     --------
c
c**   interface.
c     ----------
c          from *orodrag*
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
c     ==== outputs ===
c
c        implicit arguments :   none
c        --------------------
c
c     method.
c     -------
c
c
c     externals.
c     ----------
c
c
c     reference.
c     ----------
c
c        see ecmwf research department documentation of the "i.f.s."
c
c     author.
c     -------
c
c     modifications.
c     --------------
c     f.lott  for the new-gwdrag scheme november 1993
c
c-----------------------------------------------------------------------
      use dimens_m
      use dimphy
      use YOMCST
            use yoegwd
      implicit none
c


c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
      integer nlon
      integer jl, jk
      real zdelp

      integer kkcrit(nlon),kkcrith(nlon),kcrit(nlon),
     *        ktest(nlon),kkenvh(nlon)

c
      real paphm1(nlon,klev+1),papm1(nlon,klev),pum1(nlon,klev),
     *     pvm1(nlon,klev),ptm1(nlon,klev),pgeom1(nlon,klev),
     *     prho(nlon,klev+1),pri(nlon,klev+1),pstab(nlon,klev+1),
     *     ptau(nlon,klev+1),pvph(nlon,klev+1),ppsi(nlon,klev+1),
     *     pzdep(nlon,klev)
       real pulow(nlon),pvlow(nlon),ptheta(nlon),pgamma(nlon),pnu(nlon),
     *     pd1(nlon),pd2(nlon),pdmod(nlon)
      real pstd(nlon),pmea(nlon),ppic(nlon),pval(nlon)
c
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------
c
c
      integer ilevm1, ilevm2, ilevh
      real zcons1, zcons2,zcons3, zhgeo
      real zu, zphi, zvt1,zvt2, zst, zvar, zdwind, zwind
      real zstabm, zstabp, zrhom,  zrhop, alpha
      real zggeenv, zggeom1,zgvar 
      logical lo 
      logical ll1(klon,klev+1)
      integer kknu(klon),kknu2(klon),kknub(klon),kknul(klon),
     *        kentp(klon),ncount(klon)  
c
      real zhcrit(klon,klev),zvpf(klon,klev),
     *     zdp(klon,klev)
      real znorm(klon),zb(klon),zc(klon),
     *      zulow(klon),zvlow(klon),znup(klon),znum(klon)
c
c     ------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
c     print *,' entree gwsetup'
 100  continue
c
c     ------------------------------------------------------------------
c
c*         1.1   computational constants
c                -----------------------
c
 110  continue
c
      ilevm1=klev-1
      ilevm2=klev-2
      ilevh =klev/3
c
      zcons1=1./rd
cold  zcons2=g**2/cpd
      zcons2=rg**2/rcpd
cold  zcons3=1.5*api
      zcons3=1.5*rpi
c
c
c     ------------------------------------------------------------------
c
c*         2.
c                --------------
c
 200  continue
c
c     ------------------------------------------------------------------
c
c*         2.1     define low level wind, project winds in plane of
c*                 low level wind, determine sector in which to take
c*                 the variance and set indicator for critical levels.
c
c
c
      do 2001 jl=1,klon
      kknu(jl)    =klev
      kknu2(jl)   =klev
      kknub(jl)   =klev
      kknul(jl)   =klev
      pgamma(jl) =max(pgamma(jl),gtsec)
      ll1(jl,klev+1)=.false.
 2001 continue
c
c Ajouter une initialisation (L. Li, le 23fev99):
c
      do jk=klev,ilevh,-1
      do jl=1,klon
      ll1(jl,jk)= .FALSE.
      ENDDO
      ENDDO
c
c*      define top of low level flow
c       ----------------------------
      do 2002 jk=klev,ilevh,-1
      do 2003 jl=1,klon
      lo=(paphm1(jl,jk)/paphm1(jl,klev+1)).ge.gsigcr
      if(lo) then
        kkcrit(jl)=jk
      endif
      zhcrit(jl,jk)=ppic(jl)
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        kknu(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknu(jl)=ilevh
 2003 continue
 2002 continue
      do 2004 jk=klev,ilevh,-1
      do 2005 jl=1,klon
      zhcrit(jl,jk)=ppic(jl)-pval(jl)
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        kknu2(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknu2(jl)=ilevh
 2005 continue
 2004 continue
      do 2006 jk=klev,ilevh,-1
      do 2007 jl=1,klon
      zhcrit(jl,jk)=amax1(ppic(jl)-pmea(jl),pmea(jl)-pval(jl))
      zhgeo=pgeom1(jl,jk)/rg
      ll1(jl,jk)=(zhgeo.gt.zhcrit(jl,jk))
      if(ll1(jl,jk).neqv.ll1(jl,jk+1)) then
        kknub(jl)=jk
      endif
      if(.not.ll1(jl,ilevh))kknub(jl)=ilevh
 2007 continue
 2006 continue
c
      do 2010 jl=1,klon  
      kknu(jl)=min(kknu(jl),nktopg)
      kknu2(jl)=min(kknu2(jl),nktopg)
      kknub(jl)=min(kknub(jl),nktopg)
      kknul(jl)=klev
 2010 continue      
c

 210  continue
c
c
cc*     initialize various arrays
c
      do 2107 jl=1,klon
      prho(jl,klev+1)  =0.0
      pstab(jl,klev+1) =0.0
      pstab(jl,1)      =0.0
      pri(jl,klev+1)   =9999.0
      ppsi(jl,klev+1)  =0.0
      pri(jl,1)        =0.0
      pvph(jl,1)       =0.0
      pulow(jl)        =0.0
      pvlow(jl)        =0.0
      zulow(jl)        =0.0
      zvlow(jl)        =0.0
      kkcrith(jl)      =klev
      kkenvh(jl)       =klev
      kentp(jl)        =klev
      kcrit(jl)        =1
      ncount(jl)       =0
      ll1(jl,klev+1)   =.false.
 2107 continue
c
c*     define low-level flow
c      ---------------------
c
      do 223 jk=klev,2,-1
      do 222 jl=1,klon
      if(ktest(jl).eq.1) then
        zdp(jl,jk)=papm1(jl,jk)-papm1(jl,jk-1)
        prho(jl,jk)=2.*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
        pstab(jl,jk)=2.*zcons2/(ptm1(jl,jk)+ptm1(jl,jk-1))*
     *  (1.-rcpd*prho(jl,jk)*(ptm1(jl,jk)-ptm1(jl,jk-1))/zdp(jl,jk))
        pstab(jl,jk)=max(pstab(jl,jk),gssec)
      endif
  222 continue
  223 continue
c
c********************************************************************
c
c*     define blocked flow
c      -------------------
      do 2115 jk=klev,ilevh,-1
      do 2116 jl=1,klon
      if(jk.ge.kknub(jl).and.jk.le.kknul(jl)) then
        pulow(jl)=pulow(jl)+pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        pvlow(jl)=pvlow(jl)+pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
      end if
 2116 continue
 2115 continue
      do 2110 jl=1,klon
      pulow(jl)=pulow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknub(jl)))
      pvlow(jl)=pvlow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknub(jl)))
      znorm(jl)=max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
      pvph(jl,klev+1)=znorm(jl)
 2110 continue
c
c*******  setup orography axes and define plane of profiles  *******
c
      do 2112 jl=1,klon
      lo=(pulow(jl).lt.gvsec).and.(pulow(jl).ge.-gvsec)
      if(lo) then
        zu=pulow(jl)+2.*gvsec
      else
        zu=pulow(jl)
      endif
      zphi=atan(pvlow(jl)/zu)
      ppsi(jl,klev+1)=ptheta(jl)*rpi/180.-zphi
      zb(jl)=1.-0.18*pgamma(jl)-0.04*pgamma(jl)**2
      zc(jl)=0.48*pgamma(jl)+0.3*pgamma(jl)**2
      pd1(jl)=zb(jl)-(zb(jl)-zc(jl))*(sin(ppsi(jl,klev+1))**2)
      pd2(jl)=(zb(jl)-zc(jl))*sin(ppsi(jl,klev+1))*cos(ppsi(jl,klev+1))
      pdmod(jl)=sqrt(pd1(jl)**2+pd2(jl)**2)
 2112 continue
c
c  ************ define flow in plane of lowlevel stress *************
c
      do 213 jk=1,klev
      do 212 jl=1,klon
      if(ktest(jl).eq.1)  then
        zvt1       =pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk)
        zvt2       =-pvlow(jl)*pum1(jl,jk)+pulow(jl)*pvm1(jl,jk)
        zvpf(jl,jk)=(zvt1*pd1(jl)+zvt2*pd2(jl))/(znorm(jl)*pdmod(jl))
      endif
      ptau(jl,jk)  =0.0
      pzdep(jl,jk) =0.0
      ppsi(jl,jk)  =0.0
      ll1(jl,jk)   =.false.
  212 continue
  213 continue
      do 215 jk=2,klev
      do 214 jl=1,klon
      if(ktest(jl).eq.1) then
        zdp(jl,jk)=papm1(jl,jk)-papm1(jl,jk-1)
        pvph(jl,jk)=((paphm1(jl,jk)-papm1(jl,jk-1))*zvpf(jl,jk)+
     *            (papm1(jl,jk)-paphm1(jl,jk))*zvpf(jl,jk-1))
     *            /zdp(jl,jk)
        if(pvph(jl,jk).lt.gvsec) then
          pvph(jl,jk)=gvsec
          kcrit(jl)=jk
        endif
      endif
  214 continue
  215 continue
c
c
c*         2.2     brunt-vaisala frequency and density at half levels.
c
  220 continue
c
      do 2211 jk=ilevh,klev
      do 221 jl=1,klon
      if(ktest(jl).eq.1) then
      if(jk.ge.(kknub(jl)+1).and.jk.le.kknul(jl)) then
           zst=zcons2/ptm1(jl,jk)*(1.-rcpd*prho(jl,jk)*
     *                   (ptm1(jl,jk)-ptm1(jl,jk-1))/zdp(jl,jk))
           pstab(jl,klev+1)=pstab(jl,klev+1)+zst*zdp(jl,jk)
           pstab(jl,klev+1)=max(pstab(jl,klev+1),gssec)
           prho(jl,klev+1)=prho(jl,klev+1)+paphm1(jl,jk)*2.*zdp(jl,jk)
     *                   *zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
      endif
      endif
  221 continue
 2211 continue
c
      do 2212 jl=1,klon
        pstab(jl,klev+1)=pstab(jl,klev+1)/(papm1(jl,kknul(jl))
     *                                          -papm1(jl,kknub(jl)))
        prho(jl,klev+1)=prho(jl,klev+1)/(papm1(jl,kknul(jl))
     *                                          -papm1(jl,kknub(jl)))
        zvar=pstd(jl)
 2212 continue
c
c*         2.3     mean flow richardson number.
c*                 and critical height for froude layer
c
  230 continue
c
      do 232 jk=2,klev
      do 231 jl=1,klon
      if(ktest(jl).eq.1) then
        zdwind=max(abs(zvpf(jl,jk)-zvpf(jl,jk-1)),gvsec)
        pri(jl,jk)=pstab(jl,jk)*(zdp(jl,jk)
     *          /(rg*prho(jl,jk)*zdwind))**2
        pri(jl,jk)=max(pri(jl,jk),grcrit)
      endif
  231 continue
  232 continue
  
c
c
c*      define top of 'envelope' layer
c       ----------------------------

      do 233 jl=1,klon
      pnu (jl)=0.0
      znum(jl)=0.0
 233  continue
      
      do 234 jk=2,klev-1
      do 234 jl=1,klon
      
      if(ktest(jl).eq.1) then
       
      if (jk.ge.kknub(jl)) then
          
            znum(jl)=pnu(jl)
            zwind=(pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/
     *            max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
            zwind=max(sqrt(zwind**2),gvsec)
            zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
            zstabm=sqrt(max(pstab(jl,jk  ),gssec))
            zstabp=sqrt(max(pstab(jl,jk+1),gssec))
            zrhom=prho(jl,jk  )
            zrhop=prho(jl,jk+1)
            pnu(jl) = pnu(jl) + (zdelp/rg)*
     *            ((zstabp/zrhop+zstabm/zrhom)/2.)/zwind     
            if((znum(jl).le.gfrcrit).and.(pnu(jl).gt.gfrcrit)
     *                          .and.(kkenvh(jl).eq.klev))
     *      kkenvh(jl)=jk
     
      endif    

      endif
      
 234  continue
      
c  calculation of a dynamical mixing height for the breaking
c  of gravity waves:

              
      do 235 jl=1,klon
      znup(jl)=0.0
      znum(jl)=0.0
 235  continue

      do 236 jk=klev-1,2,-1
      do 236 jl=1,klon
      
      if(ktest(jl).eq.1) then

            znum(jl)=znup(jl)
            zwind=(pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/
     *            max(sqrt(pulow(jl)**2+pvlow(jl)**2),gvsec)
            zwind=max(sqrt(zwind**2),gvsec)
            zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
            zstabm=sqrt(max(pstab(jl,jk  ),gssec))
            zstabp=sqrt(max(pstab(jl,jk+1),gssec))
            zrhom=prho(jl,jk  )
            zrhop=prho(jl,jk+1)
            znup(jl) = znup(jl) + (zdelp/rg)*
     *            ((zstabp/zrhop+zstabm/zrhom)/2.)/zwind     
            if((znum(jl).le.rpi/2.).and.(znup(jl).gt.rpi/2.)
     *                          .and.(kkcrith(jl).eq.klev))
     *      kkcrith(jl)=jk
     
      endif
      
 236  continue
 
      do 237 jl=1,klon
      kkcrith(jl)=min0(kkcrith(jl),kknu2(jl))
      kkcrith(jl)=max0(kkcrith(jl),ilevh*2)
 237  continue         
c
c     directional info for flow blocking ************************* 
c
      do 251 jk=ilevh,klev    
      do 252 jl=1,klon
      if(jk.ge.kkenvh(jl)) then
      lo=(pum1(jl,jk).lt.gvsec).and.(pum1(jl,jk).ge.-gvsec)
      if(lo) then
        zu=pum1(jl,jk)+2.*gvsec
      else
        zu=pum1(jl,jk)
      endif
       zphi=atan(pvm1(jl,jk)/zu)
       ppsi(jl,jk)=ptheta(jl)*rpi/180.-zphi
      end if
 252  continue
 251  continue
c      forms the vertical 'leakiness' **************************

      alpha=3.
      
      do 254  jk=ilevh,klev
      do 253  jl=1,klon
      if(jk.ge.kkenvh(jl)) then
        zggeenv=amax1(1.,
     *          (pgeom1(jl,kkenvh(jl))+pgeom1(jl,kkenvh(jl)-1))/2.)      
        zggeom1=amax1(pgeom1(jl,jk),1.)
        zgvar=amax1(pstd(jl)*rg,1.)     
cmod    pzdep(jl,jk)=sqrt((zggeenv-zggeom1)/(zggeom1+zgvar))      
        pzdep(jl,jk)=(pgeom1(jl,kkenvh(jl)-1)-pgeom1(jl,  jk))/
     *               (pgeom1(jl,kkenvh(jl)-1)-pgeom1(jl,klev))
      end if
 253  continue
 254  continue

 260  continue

      return
      end
      SUBROUTINE gwstress
     *         (  nlon  , nlev
     *         , ktest, kcrit, kkenvh
     *         , kknu
     *         , prho  , pstab , pvph  , pstd, psig
     *         , pmea , ppic  , ptau  
     *         , pgeom1 , pdmod )
c
c**** *gwstress*
c
c     purpose.
c     --------
c
c**   interface.
c     ----------
c     call *gwstress*  from *gwdrag*
c
c        explicit arguments :
c        --------------------
c     ==== inputs ===
c     ==== outputs ===
c
c        implicit arguments :   none
c        --------------------
c
c     method.
c     -------
c
c
c     externals.
c     ----------
c
c
c     reference.
c     ----------
c
c        see ecmwf research department documentation of the "i.f.s."
c
c     author.
c     -------
c
c     modifications.
c     --------------
c     f. lott put the new gwd on ifs      22/11/93
c
c-----------------------------------------------------------------------
      use dimens_m
      use dimphy
      use YOMCST
            use yoegwd
      implicit none

c-----------------------------------------------------------------------
c
c*       0.1   arguments
c              ---------
c
      integer nlon, nlev
      integer kcrit(nlon),
     *        ktest(nlon),kkenvh(nlon),kknu(nlon)
c
      real prho(nlon,nlev+1),pstab(nlon,nlev+1),ptau(nlon,nlev+1),
     *     pvph(nlon,nlev+1),
     *     pgeom1(nlon,nlev),pstd(nlon)
c
      real psig(nlon)
      real pmea(nlon),ppic(nlon)
      real pdmod(nlon)
c
c-----------------------------------------------------------------------
c
c*       0.2   local arrays
c              ------------
      integer jl
      real zblock, zvar, zeff
      logical lo
c
c-----------------------------------------------------------------------
c
c*       0.3   functions
c              ---------
c     ------------------------------------------------------------------
c
c*         1.    initialization
c                --------------
c
 100  continue
c
c*         3.1     gravity wave stress.
c
  300 continue
c
c
      do 301 jl=1,klon
      if(ktest(jl).eq.1) then
      
c  effective mountain height above the blocked flow
  
         if(kkenvh(jl).eq.klev)then
         zblock=0.0 
         else
         zblock=(pgeom1(jl,kkenvh(jl))+pgeom1(jl,kkenvh(jl)+1))/2./rg          
         endif
      
        zvar=ppic(jl)-pmea(jl)
        zeff=amax1(0.,zvar-zblock)

        ptau(jl,klev+1)=prho(jl,klev+1)*gkdrag*psig(jl)*zeff**2
     *    /4./pstd(jl)*pvph(jl,klev+1)*pdmod(jl)*sqrt(pstab(jl,klev+1))

c  too small value of stress or  low level flow include critical level
c  or low level flow:  gravity wave stress nul.
                
        lo=(ptau(jl,klev+1).lt.gtsec).or.(kcrit(jl).ge.kknu(jl))
     *      .or.(pvph(jl,klev+1).lt.gvcrit)
c       if(lo) ptau(jl,klev+1)=0.0
      
      else
      
          ptau(jl,klev+1)=0.0
          
      endif
      
  301 continue
c
      return
      end
      SUBROUTINE GWPROFIL
     *         ( NLON, NLEV
     *         , kgwd, kdx , ktest
     *         , KKCRITH, KCRIT
     *         , PAPHM1, PRHO   , PSTAB  , PVPH , PRI , PTAU
     *         , pdmod   , psig , pvar)

C**** *GWPROFIL*
C
C     PURPOSE.
C     --------
C
C**   INTERFACE.
C     ----------
C          FROM *GWDRAG*
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C     ==== OUTPUTS ===
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD:
C     -------
C     THE STRESS PROFILE FOR GRAVITY WAVES IS COMPUTED AS FOLLOWS:
C     IT IS CONSTANT (NO GWD) AT THE LEVELS BETWEEN THE GROUND
C     AND THE TOP OF THE BLOCKED LAYER (KKENVH).
C     IT DECREASES LINEARLY WITH HEIGHTS FROM THE TOP OF THE 
C     BLOCKED LAYER TO 3*VAROR (kKNU), TO SIMULATES LEE WAVES OR 
C     NONLINEAR GRAVITY WAVE BREAKING.
C     ABOVE IT IS CONSTANT, EXCEPT WHEN THE WAVE ENCOUNTERS A CRITICAL
C     LEVEL (KCRIT) OR WHEN IT BREAKS.
C     
C
C
C     EXTERNALS.
C     ----------
C
C
C     REFERENCE.
C     ----------
C
C        SEE ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S."
C
C     AUTHOR.
C     -------
C
C     MODIFICATIONS.
C     --------------
C     PASSAGE OF THE NEW GWDRAG TO I.F.S. (F. LOTT, 22/11/93)
C-----------------------------------------------------------------------
      use dimens_m
      use dimphy
      use YOMCST
            use yoegwd
      implicit none
C

C


C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      integer nlon,nlev
      INTEGER KKCRITH(NLON),KCRIT(NLON)
     *       ,kdx(nlon) , ktest(nlon)

C
      REAL PAPHM1(NLON,NLEV+1), PSTAB(NLON,NLEV+1),
     *     PRHO  (NLON,NLEV+1), PVPH (NLON,NLEV+1),
     *     PRI   (NLON,NLEV+1), PTAU(NLON,NLEV+1)
     
      REAL pdmod (NLON) , psig(NLON),
     *     pvar(NLON)
     
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
      integer ilevh, ji, kgwd, jl, jk
      real zsqr, zalfa, zriw, zdel, zb, zalpha,zdz2n
      real zdelp, zdelpt 
      REAL ZDZ2 (KLON,KLEV) , ZNORM(KLON) , zoro(KLON)
      REAL ZTAU (KLON,KLEV+1)
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
c      print *,' entree gwprofil' 
 100  CONTINUE
C
C
C*    COMPUTATIONAL CONSTANTS.
C     ------------- ----------
C
      ilevh=KLEV/3
C
c     DO 400 ji=1,kgwd
c     jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      DO 400 jl=1,klon
      if (ktest(jl).eq.1) then
      Zoro(JL)=Psig(JL)*Pdmod(JL)/4./max(pvar(jl),1.0)
      ZTAU(JL,KLEV+1)=PTAU(JL,KLEV+1)
      endif
  400 CONTINUE
  
C
      DO 430 JK=KLEV,2,-1
C
C
C*         4.1    CONSTANT WAVE STRESS UNTIL TOP OF THE
C                 BLOCKING LAYER.
  410 CONTINUE
C
c     DO 411 ji=1,kgwd
c     jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      do 411 jl=1,klon
      if (ktest(jl).eq.1) then
           IF(JK.GT.KKCRITH(JL)) THEN
           PTAU(JL,JK)=ZTAU(JL,KLEV+1)
C          ENDIF
C          IF(JK.EQ.KKCRITH(JL)) THEN
           ELSE                    
           PTAU(JL,JK)=GRAHILO*ZTAU(JL,KLEV+1)
           ENDIF
      endif
 411  CONTINUE             
C
C*         4.15   CONSTANT SHEAR STRESS UNTIL THE TOP OF THE
C                 LOW LEVEL FLOW LAYER.
 415  CONTINUE
C        
C
C*         4.2    WAVE DISPLACEMENT AT NEXT LEVEL.
C
  420 CONTINUE
C
c     DO 421 ji=1,kgwd
c     jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      do 421 jl=1,klon
      if(ktest(jl).eq.1) then
      IF(JK.LT.KKCRITH(JL)) THEN
      ZNORM(JL)=gkdrag*PRHO(JL,JK)*SQRT(PSTAB(JL,JK))*PVPH(JL,JK)
     *                                                    *zoro(jl)
      ZDZ2(JL,JK)=PTAU(JL,JK+1)/max(ZNORM(JL),gssec)
      ENDIF
      endif
  421 CONTINUE
C
C*         4.3    WAVE RICHARDSON NUMBER, NEW WAVE DISPLACEMENT
C*                AND STRESS:  BREAKING EVALUATION AND CRITICAL 
C                 LEVEL
C
                          
c     DO 431 ji=1,kgwd
c     jl=Kdx(ji)
c  Modif vectorisation 02/04/2004
      do 431 jl=1,klon
      if(ktest(jl).eq.1) then

          IF(JK.LT.KKCRITH(JL)) THEN
          IF((PTAU(JL,JK+1).LT.GTSEC).OR.(JK.LE.KCRIT(JL))) THEN
            PTAU(JL,JK)=0.0
          ELSE
               ZSQR=SQRT(PRI(JL,JK))
               ZALFA=SQRT(PSTAB(JL,JK)*ZDZ2(JL,JK))/PVPH(JL,JK)
               ZRIW=PRI(JL,JK)*(1.-ZALFA)/(1+ZALFA*ZSQR)**2
               IF(ZRIW.LT.GRCRIT) THEN
                 ZDEL=4./ZSQR/GRCRIT+1./GRCRIT**2+4./GRCRIT
                 ZB=1./GRCRIT+2./ZSQR
                 ZALPHA=0.5*(-ZB+SQRT(ZDEL))
                 ZDZ2N=(PVPH(JL,JK)*ZALPHA)**2/PSTAB(JL,JK)
                 PTAU(JL,JK)=ZNORM(JL)*ZDZ2N
               ELSE
                 PTAU(JL,JK)=ZNORM(JL)*ZDZ2(JL,JK)
               ENDIF
            PTAU(JL,JK)=MIN(PTAU(JL,JK),PTAU(JL,JK+1))
          ENDIF
          ENDIF
      endif
  431 CONTINUE
  
  430 CONTINUE
  440 CONTINUE
  
C  REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

c     DO 530 ji=1,kgwd
c     jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      do 530 jl=1,klon
      if(ktest(jl).eq.1) then
      ZTAU(JL,KKCRITH(JL))=PTAU(JL,KKCRITH(JL))
      ZTAU(JL,NSTRA)=PTAU(JL,NSTRA)
      endif
 530  CONTINUE      

      DO 531 JK=1,KLEV
      
c     DO 532 ji=1,kgwd
c     jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      do 532 jl=1,klon
      if(ktest(jl).eq.1) then

                
         IF(JK.GT.KKCRITH(JL))THEN

          ZDELP=PAPHM1(JL,JK)-PAPHM1(JL,KLEV+1    )
          ZDELPT=PAPHM1(JL,KKCRITH(JL))-PAPHM1(JL,KLEV+1    )
          PTAU(JL,JK)=ZTAU(JL,KLEV+1    ) +
     .                (ZTAU(JL,KKCRITH(JL))-ZTAU(JL,KLEV+1    ) )*
     .                ZDELP/ZDELPT
     
        ENDIF
            
      endif
 532  CONTINUE    
 
C  REORGANISATION IN THE STRATOSPHERE

c     DO 533 ji=1,kgwd
c     jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      do 533 jl=1,klon
      if(ktest(jl).eq.1) then


         IF(JK.LT.NSTRA)THEN

          ZDELP =PAPHM1(JL,NSTRA)
          ZDELPT=PAPHM1(JL,JK)
          PTAU(JL,JK)=ZTAU(JL,NSTRA)*ZDELPT/ZDELP 

        ENDIF

      endif
 533  CONTINUE

C REORGANISATION IN THE TROPOSPHERE

c      DO 534 ji=1,kgwd
c      jl=kdx(ji)
c  Modif vectorisation 02/04/2004
      do 534 jl=1,klon
      if(ktest(jl).eq.1) then


         IF(JK.LT.KKCRITH(JL).AND.JK.GT.NSTRA)THEN

           ZDELP=PAPHM1(JL,JK)-PAPHM1(JL,KKCRITH(JL))
           ZDELPT=PAPHM1(JL,NSTRA)-PAPHM1(JL,KKCRITH(JL))
           PTAU(JL,JK)=ZTAU(JL,KKCRITH(JL)) +
     *                 (ZTAU(JL,NSTRA)-ZTAU(JL,KKCRITH(JL)))*ZDELP
     .                                                     /ZDELPT

       ENDIF
      endif
 534   CONTINUE

 
 531  CONTINUE        


      RETURN
      END
      SUBROUTINE lift_noro (nlon,nlev,dtime,paprs,pplay,      
     e                   plat,pmea,pstd, ppic,
     e                   ktest,
     e                   t, u, v,
     s                   pulow, pvlow, pustr, pvstr,
     s                   d_t, d_u, d_v)
c
      use dimens_m
      use dimphy
      use YOMCST
      IMPLICIT none
c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 19950201
c Objet: Frottement de la montagne Interface
c======================================================================
c Arguments:
c dtime---input-R- pas d'integration (s)
c paprs---input-R-pression pour chaque inter-couche (en Pa)
c pplay---input-R-pression pour le mileu de chaque couche (en Pa)
c t-------input-R-temperature (K)
c u-------input-R-vitesse horizontale (m/s)
c v-------input-R-vitesse horizontale (m/s)
c
c d_t-----output-R-increment de la temperature
c d_u-----output-R-increment de la vitesse u
c d_v-----output-R-increment de la vitesse v
c======================================================================
c
c ARGUMENTS
c
      INTEGER nlon,nlev
      REAL, intent(in):: dtime
      REAL, intent(in):: paprs(klon,klev+1)
      REAL, intent(in):: pplay(klon,klev)
      REAL, intent(in):: plat(nlon)
      real pmea(nlon)
      REAL pstd(nlon)
      REAL ppic(nlon)
      REAL pulow(nlon),pvlow(nlon),pustr(nlon),pvstr(nlon)
      REAL t(nlon,nlev), u(nlon,nlev), v(nlon,nlev)
      REAL d_t(nlon,nlev), d_u(nlon,nlev), d_v(nlon,nlev)
c
      INTEGER i, k, ktest(nlon)
c
c Variables locales:
c
      REAL zgeom(klon,klev)
      REAL pdtdt(klon,klev), pdudt(klon,klev), pdvdt(klon,klev)
      REAL pt(klon,klev), pu(klon,klev), pv(klon,klev)
      REAL papmf(klon,klev),papmh(klon,klev+1)
c
c initialiser les variables de sortie (pour securite)
c
      DO i = 1,klon
         pulow(i) = 0.0
         pvlow(i) = 0.0
         pustr(i) = 0.0
         pvstr(i) = 0.0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
         pdudt(i,k)=0.0
         pdvdt(i,k)=0.0
         pdtdt(i,k)=0.0
      ENDDO
      ENDDO
c
c preparer les variables d'entree (attention: l'ordre des niveaux 
c verticaux augmente du haut vers le bas)
c
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1) 
         pu(i,k) = u(i,klev-k+1)
         pv(i,k) = v(i,klev-k+1)
         papmf(i,k) = pplay(i,klev-k+1)
      ENDDO
      ENDDO
      DO k = 1, klev+1
      DO i = 1, klon
         papmh(i,k) = paprs(i,klev-k+2)
      ENDDO
      ENDDO
      DO i = 1, klon
         zgeom(i,klev) = RD * pt(i,klev)
     .                  * LOG(papmh(i,klev+1)/papmf(i,klev))
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         zgeom(i,k) = zgeom(i,k+1) + RD * (pt(i,k)+pt(i,k+1))/2.0
     .               * LOG(papmf(i,k+1)/papmf(i,k))
      ENDDO
      ENDDO
c
c appeler la routine principale
c
      CALL OROLIFT(klon,klev,ktest,
     .            dtime,
     .            papmh, zgeom,
     .            pt, pu, pv,
     .            plat,pmea, pstd, ppic,
     .            pulow,pvlow,
     .            pdudt,pdvdt,pdtdt)
C
      DO k = 1, klev
      DO i = 1, klon
         d_u(i,klev+1-k) = dtime*pdudt(i,k)
         d_v(i,klev+1-k) = dtime*pdvdt(i,k)
         d_t(i,klev+1-k) = dtime*pdtdt(i,k)
         pustr(i)        = pustr(i)
cIM BUG .                 +RG*pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))
     .                    +pdudt(i,k)*(papmh(i,k+1)-papmh(i,k))/RG
         pvstr(i)        = pvstr(i)
cIM BUG .                 +RG*pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))
     .                    +pdvdt(i,k)*(papmh(i,k+1)-papmh(i,k))/RG
      ENDDO
      ENDDO
c
      RETURN
      END
      SUBROUTINE OROLIFT( NLON,NLEV
     I                 , KTEST
     R                 , PTSPHY
     R                 , PAPHM1,PGEOM1,PTM1,PUM1,PVM1
     R                 , PLAT
     R                 , PMEA, PVAROR, ppic
C OUTPUTS
     R                 , PULOW,PVLOW
     R                 , PVOM,PVOL,PTE )

C
C**** *OROLIFT: SIMULATE THE GEOSTROPHIC LIFT.
C
C     PURPOSE.
C     --------
C
C**   INTERFACE.
C     ----------
C          CALLED FROM *lift_noro
C     ----------
C
C     AUTHOR.
C     -------
C     F.LOTT  LMD 22/11/95
C
      use dimens_m
      use dimphy
      use YOMCST
            use yoegwd
      implicit none
C
C
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
C
      integer nlon, nlev
      REAL  PTE(NLON,NLEV),
     *      PVOL(NLON,NLEV),
     *      PVOM(NLON,NLEV),
     *      PULOW(NLON),
     *      PVLOW(NLON)
      REAL  PUM1(NLON,NLEV),
     *      PVM1(NLON,NLEV),
     *      PTM1(NLON,NLEV)
      real, intent(in):: PLAT(NLON)
      real PMEA(NLON),
     *      PVAROR(NLON),
     *      ppic(NLON),
     *      PGEOM1(NLON,NLEV),
     *      PAPHM1(NLON,NLEV+1)
C
      INTEGER  KTEST(NLON)
      real, intent(in):: ptsphy
C-----------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
      logical lifthigh
      integer klevm1, jl, ilevh, jk
      real zcons1, ztmst, zrtmst,zpi, zhgeo
      real zdelp, zslow, zsqua, zscav, zbet
      INTEGER  
     *         IKNUB(klon),
     *         IKNUL(klon)
      LOGICAL LL1(KLON,KLEV+1)
C
      REAL   ZTAU(KLON,KLEV+1),
     *       ZTAV(KLON,KLEV+1),
     *       ZRHO(KLON,KLEV+1)
      REAL   ZDUDT(KLON),
     *       ZDVDT(KLON)
      REAL ZHCRIT(KLON,KLEV)
C-----------------------------------------------------------------------
C
C*         1.1  INITIALIZATIONS
C               ---------------

      LIFTHIGH=.FALSE.

      IF(NLON.NE.KLON.OR.NLEV.NE.KLEV)STOP
      ZCONS1=1./RD
      KLEVM1=KLEV-1
      ZTMST=PTSPHY
      ZRTMST=1./ZTMST
      ZPI=ACOS(-1.)
C
      DO 1001 JL=1,klon
      ZRHO(JL,KLEV+1)  =0.0
      PULOW(JL)        =0.0
      PVLOW(JL)        =0.0
      iknub(JL)   =klev
      iknul(JL)   =klev
      ilevh=klev/3
      ll1(jl,klev+1)=.false.
      DO 1000 JK=1,KLEV
      PVOM(JL,JK)=0.0
      PVOL(JL,JK)=0.0
      PTE (JL,JK)=0.0
 1000 CONTINUE
 1001 CONTINUE

C
C*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
C*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
C*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
C
C
C
      DO 2006 JK=KLEV,1,-1
      DO 2007 JL=1,klon
      IF(KTEST(JL).EQ.1) THEN
      ZHCRIT(JL,JK)=amax1(Ppic(JL)-pmea(JL),100.)
      ZHGEO=PGEOM1(JL,JK)/RG
      ll1(JL,JK)=(ZHGEO.GT.ZHCRIT(JL,JK))
      IF(ll1(JL,JK).neqv.ll1(JL,JK+1)) THEN
        iknub(JL)=JK
      ENDIF
      ENDIF
 2007 CONTINUE
 2006 CONTINUE
C
      do 2010 jl=1,klon
      IF(KTEST(JL).EQ.1) THEN
      iknub(jl)=max(iknub(jl),klev/2)
      iknul(jl)=max(iknul(jl),2*klev/3)
      if(iknub(jl).gt.nktopg) iknub(jl)=nktopg
      if(iknub(jl).eq.nktopg) iknul(jl)=klev
      if(iknub(jl).eq.iknul(jl)) iknub(jl)=iknul(jl)-1
      ENDIF
 2010 continue

C     do 2011 jl=1,klon
C     IF(KTEST(JL).EQ.1) THEN
C       print *,' iknul= ',iknul(jl),'  iknub=',iknub(jl)
C     ENDIF
C2011 continue

C     PRINT *,'  DANS OROLIFT: 2010'

      DO 223 JK=KLEV,2,-1
      DO 222 JL=1,klon
        ZRHO(JL,JK)=2.*PAPHM1(JL,JK)*ZCONS1/(PTM1(JL,JK)+PTM1(JL,JK-1))
  222 CONTINUE
  223 CONTINUE
C     PRINT *,'  DANS OROLIFT: 223'

C********************************************************************
C
C*     DEFINE LOW LEVEL FLOW
C      -------------------
      DO 2115 JK=klev,1,-1
      DO 2116 JL=1,klon
      IF(KTEST(JL).EQ.1) THEN
      if(jk.ge.iknub(jl).and.jk.le.iknul(jl)) then
        pulow(JL)=pulow(JL)+PUM1(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
        pvlow(JL)=pvlow(JL)+PVM1(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
        zrho(JL,klev+1)=zrho(JL,klev+1)
     *                 +zrho(JL,JK)*(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
      end if
      ENDIF
 2116 CONTINUE
 2115 CONTINUE
      DO 2110 JL=1,klon
      IF(KTEST(JL).EQ.1) THEN
      pulow(JL)=pulow(JL)/(PAPHM1(JL,iknul(jl)+1)-PAPHM1(JL,iknub(jl)))
      pvlow(JL)=pvlow(JL)/(PAPHM1(JL,iknul(jl)+1)-PAPHM1(JL,iknub(jl)))
      zrho(JL,klev+1)=zrho(JL,klev+1)
     *               /(PAPHM1(JL,iknul(jl)+1)-PAPHM1(JL,iknub(jl)))
      ENDIF
 2110 CONTINUE


200   CONTINUE

C***********************************************************
C
C*         3.      COMPUTE MOUNTAIN LIFT
C
  300 CONTINUE
C
      DO 301 JL=1,klon
      IF(KTEST(JL).EQ.1) THEN
       ZTAU(JL,KLEV+1)= - GKLIFT*ZRHO(JL,KLEV+1)*2.*ROMEGA*
C    *                 (2*PVAROR(JL)+PMEA(JL))*
     *                 2*PVAROR(JL)*
     *                 SIN(ZPI/180.*PLAT(JL))*PVLOW(JL)
       ZTAV(JL,KLEV+1)=   GKLIFT*ZRHO(JL,KLEV+1)*2.*ROMEGA*
C    *                 (2*PVAROR(JL)+PMEA(JL))*
     *                 2*PVAROR(JL)*
     *                 SIN(ZPI/180.*PLAT(JL))*PULOW(JL)
      ELSE
       ZTAU(JL,KLEV+1)=0.0
       ZTAV(JL,KLEV+1)=0.0
      ENDIF
301   CONTINUE

C
C*         4.      COMPUTE LIFT PROFILE         
C*                 --------------------   
C

  400 CONTINUE

      DO 401 JK=1,KLEV
      DO 401 JL=1,klon
      IF(KTEST(JL).EQ.1) THEN
      ZTAU(JL,JK)=ZTAU(JL,KLEV+1)*PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1)
      ZTAV(JL,JK)=ZTAV(JL,KLEV+1)*PAPHM1(JL,JK)/PAPHM1(JL,KLEV+1)
      ELSE
      ZTAU(JL,JK)=0.0
      ZTAV(JL,JK)=0.0
      ENDIF
401   CONTINUE
C
C
C*         5.      COMPUTE TENDENCIES.
C*                 -------------------
      IF(LIFTHIGH)THEN
C
  500 CONTINUE
C     PRINT *,'  DANS OROLIFT: 500'
C
C  EXPLICIT SOLUTION AT ALL LEVELS
C
      DO 524 JK=1,klev
      DO 523 JL=1,KLON
      IF(KTEST(JL).EQ.1) THEN
      ZDELP=PAPHM1(JL,JK+1)-PAPHM1(JL,JK)
      ZDUDT(JL)=-RG*(ZTAU(JL,JK+1)-ZTAU(JL,JK))/ZDELP
      ZDVDT(JL)=-RG*(ZTAV(JL,JK+1)-ZTAV(JL,JK))/ZDELP
      ENDIF  
  523 CONTINUE
  524 CONTINUE
C
C  PROJECT PERPENDICULARLY TO U NOT TO DESTROY ENERGY
C
      DO 530 JK=1,klev
      DO 530 JL=1,KLON
      IF(KTEST(JL).EQ.1) THEN

        ZSLOW=SQRT(PULOW(JL)**2+PVLOW(JL)**2)
        ZSQUA=AMAX1(SQRT(PUM1(JL,JK)**2+PVM1(JL,JK)**2),GVSEC)
        ZSCAV=-ZDUDT(JL)*PVM1(JL,JK)+ZDVDT(JL)*PUM1(JL,JK)
        IF(ZSQUA.GT.GVSEC)THEN
          PVOM(JL,JK)=-ZSCAV*PVM1(JL,JK)/ZSQUA**2
          PVOL(JL,JK)= ZSCAV*PUM1(JL,JK)/ZSQUA**2
        ELSE
          PVOM(JL,JK)=0.0
          PVOL(JL,JK)=0.0      
        ENDIF  
        ZSQUA=SQRT(PUM1(JL,JK)**2+PUM1(JL,JK)**2)               
        IF(ZSQUA.LT.ZSLOW)THEN
          PVOM(JL,JK)=ZSQUA/ZSLOW*PVOM(JL,JK)
          PVOL(JL,JK)=ZSQUA/ZSLOW*PVOL(JL,JK)
        ENDIF 

      ENDIF  
530   CONTINUE
C
C  6.  LOW LEVEL LIFT, SEMI IMPLICIT:
C  ----------------------------------

      ELSE

        DO 601 JL=1,KLON
        IF(KTEST(JL).EQ.1) THEN
          DO JK=KLEV,IKNUB(JL),-1
          ZBET=GKLIFT*2.*ROMEGA*SIN(ZPI/180.*PLAT(JL))*ztmst*
     *        (PGEOM1(JL,IKNUB(JL)-1)-PGEOM1(JL,  JK))/
     *        (PGEOM1(JL,IKNUB(JL)-1)-PGEOM1(JL,KLEV))
          ZDUDT(JL)=-PUM1(JL,JK)/ztmst/(1+ZBET**2)
          ZDVDT(JL)=-PVM1(JL,JK)/ztmst/(1+ZBET**2)
          PVOM(JL,JK)= ZBET**2*ZDUDT(JL) - ZBET   *ZDVDT(JL)
          PVOL(JL,JK)= ZBET   *ZDUDT(JL) + ZBET**2*ZDVDT(JL)    
          ENDDO
        ENDIF
 601    CONTINUE

      ENDIF

      RETURN
      END
      SUBROUTINE SUGWD(NLON,NLEV,paprs,pplay)
C
C**** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG
C
C     PURPOSE.
C     --------
C           INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
C           GRAVITY WAVE DRAG PARAMETRIZATION.
C
C**   INTERFACE.
C     ----------
C        CALL *SUGWD* FROM *SUPHEC*
C              -----        ------
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C        PSIG        : VERTICAL COORDINATE TABLE
C        NLEV        : NUMBER OF MODEL LEVELS
C
C        IMPLICIT ARGUMENTS :
C        --------------------
C        COMMON YOEGWD
C
C     METHOD.
C     -------
C        SEE DOCUMENTATION
C
C     EXTERNALS.
C     ----------
C        NONE
C
C     REFERENCE.
C     ----------
C        ECMWF Research Department documentation of the IFS
C
C     AUTHOR.
C     -------
C        MARTIN MILLER             *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 90-01-01
C     ------------------------------------------------------------------
            use yoegwd
      implicit none
C
C     -----------------------------------------------------------------
C      ----------------------------------------------------------------
C
      integer nlon,nlev, jk
      REAL, intent(in):: paprs(nlon,nlev+1)
      REAL, intent(in):: pplay(nlon,nlev)
      real zpr,zstra,zsigt,zpm1r
C
C*       1.    SET THE VALUES OF THE PARAMETERS
C              --------------------------------
C
 100  CONTINUE
C
      PRINT *,' DANS SUGWD NLEV=',NLEV
      GHMAX=10000.
C
      ZPR=100000.
      ZSTRA=0.1 
      ZSIGT=0.94
cold  ZPR=80000.
cold  ZSIGT=0.85
C
      DO 110 JK=1,NLEV
      ZPM1R=pplay(nlon/2,jk)/paprs(nlon/2,1) 
      IF(ZPM1R.GE.ZSIGT)THEN
         nktopg=JK
      ENDIF
      ZPM1R=pplay(nlon/2,jk)/paprs(nlon/2,1) 
      IF(ZPM1R.GE.ZSTRA)THEN
         NSTRA=JK
      ENDIF
  110 CONTINUE
c
c  inversion car dans orodrag on compte les niveaux a l'envers
      nktopg=nlev-nktopg+1
      nstra=nlev-nstra
      print *,' DANS SUGWD nktopg=', nktopg
      print *,' DANS SUGWD nstra=', nstra
C
      GSIGCR=0.80
C
      GKDRAG=0.2 
      GRAHILO=1.    
      GRCRIT=0.01
      GFRCRIT=1.0
      GKWAKE=0.50 
C
      GKLIFT=0.50  
      GVCRIT =0.0
C
C
C      ----------------------------------------------------------------
C
C*       2.    SET VALUES OF SECURITY PARAMETERS
C              ---------------------------------
C
 200  CONTINUE
C
      GVSEC=0.10
      GSSEC=1.E-12
C
      GTSEC=1.E-07
C
C      ----------------------------------------------------------------
C
      RETURN
      END
