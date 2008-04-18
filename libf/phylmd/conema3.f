!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/conema3.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE conema3 (dtime,paprs,pplay,t,q,u,v,tra,ntra,
     .             work1,work2,d_t,d_q,d_u,d_v,d_tra,
     .             rain, snow, kbas, ktop,
     .             upwd,dnwd,dnwdbis,bas,top,Ma,cape,tvp,rflag,
     .             pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr,
     .             qcond_incld)

      use dimens_m
      use dimphy
      use YOMCST
      use conema3_m
      use yoethf
      use fcttre
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: schema de convection de Emanuel (1991) interface
c Mai 1998: Interface modifiee pour implementation dans LMDZ
c======================================================================
c Arguments:
c dtime---input-R-pas d'integration (s)
c paprs---input-R-pression inter-couches (Pa)
c pplay---input-R-pression au milieu des couches (Pa)
c t-------input-R-temperature (K)
c q-------input-R-humidite specifique (kg/kg)
c u-------input-R-vitesse du vent zonal (m/s)
c v-------input-R-vitesse duvent meridien (m/s)
c tra-----input-R-tableau de rapport de melange des traceurs
c work*: input et output: deux variables de travail,
c                            on peut les mettre a 0 au debut
c
C d_t-----output-R-increment de la temperature
c d_q-----output-R-increment de la vapeur d'eau
c d_u-----output-R-increment de la vitesse zonale
c d_v-----output-R-increment de la vitesse meridienne
c d_tra---output-R-increment du contenu en traceurs
c rain----output-R-la pluie (mm/s)
c snow----output-R-la neige (mm/s)
c kbas----output-R-bas du nuage (integer)
c ktop----output-R-haut du nuage (integer)
c upwd----output-R-saturated updraft mass flux (kg/m**2/s)
c dnwd----output-R-saturated downdraft mass flux (kg/m**2/s)
c dnwdbis-output-R-unsaturated downdraft mass flux (kg/m**2/s)
c bas-----output-R-bas du nuage (real)
c top-----output-R-haut du nuage (real)
c Ma------output-R-flux ascendant non dilue (kg/m**2/s)
c cape----output-R-CAPE
c tvp-----output-R-virtual temperature of the lifted parcel
c rflag---output-R-flag sur le fonctionnement de convect
c pbase---output-R-pression a la base du nuage (Pa)
c bbase---output-R-buoyancy a la base du nuage (K)
c dtvpdt1-output-R-derivative of parcel virtual temp wrt T1 
c dtvpdq1-output-R-derivative of parcel virtual temp wrt Q1 
c dplcldt-output-R-derivative of the PCP pressure wrt T1
c dplcldr-output-R-derivative of the PCP pressure wrt Q1
c======================================================================
c
      INTEGER i, l,m,itra
      INTEGER ntra,ntrac !number of tracers; if no tracer transport
                         ! is needed, set ntra = 1 (or 0)
      PARAMETER (ntrac=nqmx-2)
      REAL dtime
c
      REAL d_t2(klon,klev), d_q2(klon,klev) ! sbl
      REAL d_u2(klon,klev), d_v2(klon,klev) ! sbl
      REAL em_d_t2(klev), em_d_q2(klev)     ! sbl   
      REAL em_d_u2(klev), em_d_v2(klev)     ! sbl   
c 
      REAL, intent(in):: paprs(klon,klev+1)
      real, intent(in):: pplay(klon,klev)
      REAL t(klon,klev), q(klon,klev), d_t(klon,klev), d_q(klon,klev)
      REAL u(klon,klev), v(klon,klev), tra(klon,klev,ntra)
      REAL d_u(klon,klev), d_v(klon,klev), d_tra(klon,klev,ntra)
      REAL work1(klon,klev), work2(klon,klev)
      REAL upwd(klon,klev), dnwd(klon,klev), dnwdbis(klon,klev)
      REAL rain(klon)
      REAL snow(klon)
      REAL cape(klon), tvp(klon,klev), rflag(klon)
      REAL pbase(klon), bbase(klon)
      REAL dtvpdt1(klon,klev), dtvpdq1(klon,klev)
      REAL dplcldt(klon), dplcldr(klon)
      INTEGER kbas(klon), ktop(klon)

      REAL wd(klon)
      REAL qcond_incld(klon,klev)
c
      REAL em_t(klev)
      REAL em_q(klev)
      REAL em_qs(klev)
      REAL em_u(klev), em_v(klev), em_tra(klev,ntrac)
      REAL em_ph(klev+1), em_p(klev)
      REAL em_work1(klev), em_work2(klev)
      REAL em_precip, em_d_t(klev), em_d_q(klev)
      REAL em_d_u(klev), em_d_v(klev), em_d_tra(klev,ntrac)
      REAL em_upwd(klev), em_dnwd(klev), em_dnwdbis(klev)
      REAL em_dtvpdt1(klev), em_dtvpdq1(klev)
      REAL em_dplcldt, em_dplcldr
      SAVE em_t,em_q, em_qs, em_ph, em_p, em_work1, em_work2
      SAVE em_u,em_v, em_tra
      SAVE em_d_u,em_d_v, em_d_tra
      SAVE em_precip, em_d_t, em_d_q, em_upwd, em_dnwd, em_dnwdbis
      INTEGER em_bas, em_top
      SAVE em_bas, em_top

      REAL em_wd
      REAL em_qcond(klev)
      REAL em_qcondc(klev)
c
      REAL zx_t, zx_qs, zdelta, zcor
      INTEGER iflag
      REAL sigsum
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     VARIABLES A SORTIR
cccccccccccccccccccccccccccccccccccccccccccccccccc
 
      REAL emmip(klev) !variation de flux ascnon dilue i et i+1
      SAVE emmip
      real emMke(klev)
      save emMke
      real top
      real bas
      real emMa(klev)
      save emMa
      real Ma(klon,klev)
      real Ment(klev,klev)
      real Qent(klev,klev)
      real TPS(klev),TLS(klev)
      real SIJ(klev,klev)
      real em_CAPE, em_TVP(klev)
      real em_pbase, em_bbase
      integer iw,j,k,ix,iy

c -- sb: pour schema nuages:

       integer iflagcon
       integer em_ifc(klev)
     
       real em_pradj
       real em_cldf(klev), em_cldq(klev)
       real em_ftadj(klev), em_fradj(klev)

       integer ifc(klon,klev)
       real pradj(klon)
       real cldf(klon,klev), cldq(klon,klev)
       real ftadj(klon,klev), fqadj(klon,klev)

c sb --

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 
      qcond_incld(:,:) = 0.
c
c$$$      print*,'debut conema'

      DO 999 i = 1, klon
      DO l = 1, klev+1
         em_ph(l) = paprs(i,l) / 100.0
      ENDDO
c
      DO l = 1, klev
         em_p(l) = pplay(i,l) / 100.0
         em_t(l) = t(i,l)
         em_q(l) = q(i,l)
         em_u(l) = u(i,l)
         em_v(l) = v(i,l)
         do itra = 1, ntra
          em_tra(l,itra) = tra(i,l,itra)
         enddo
c$$$      print*,'em_t',em_t
c$$$      print*,'em_q',em_q
c$$$      print*,'em_qs',em_qs
c$$$      print*,'em_u',em_u
c$$$      print*,'em_v',em_v
c$$$      print*,'em_tra',em_tra
c$$$      print*,'em_p',em_p

 
c
         zx_t = em_t(l)
         zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
         zx_qs= r2es * FOEEW(zx_t,zdelta)/em_p(l)/100.0
         zx_qs=MIN(0.5,zx_qs)
c$$$       print*,'zx_qs',zx_qs
         zcor=1./(1.-retv*zx_qs) 
         zx_qs=zx_qs*zcor
         em_qs(l) = zx_qs
c$$$      print*,'em_qs',em_qs
c
         em_work1(l) = work1(i,l)
         em_work2(l) = work2(i,l)
         emMke(l)=0
c        emMa(l)=0
c        Ma(i,l)=0
     
         em_dtvpdt1(l) = 0.
         em_dtvpdq1(l) = 0.
         dtvpdt1(i,l) = 0.
         dtvpdq1(i,l) = 0.
      ENDDO
c
      em_dplcldt = 0.
      em_dplcldr = 0.
      rain(i) = 0.0
      snow(i) = 0.0
      kbas(i) = 1
      ktop(i) = 1
c ajout SB:
      bas = 1
      top = 1
 
 
c sb3d      write(*,1792) (em_work1(m),m=1,klev)
1792  format('sig avant convect ',/,10(1X,E13.5))
c
c sb d      write(*,1793) (em_work2(m),m=1,klev)
1793  format('w avant convect ',/,10(1X,E13.5))
 
c$$$      print*,'avant convect' 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 

c     print*,'avant convect i=',i
      CALL convect3(dtime,epmax,ok_adj_ema,
     .              em_t, em_q, em_qs,em_u ,em_v ,
     .              em_tra, em_p, em_ph,
     .              klev, klev+1, klev-1,ntra, dtime, iflag,
     .              em_d_t, em_d_q,em_d_u,em_d_v,
     .              em_d_tra, em_precip,
     .              em_bas, em_top,em_upwd, em_dnwd, em_dnwdbis,
     .              em_work1, em_work2,emmip,emMke,emMa,Ment,
     .  Qent,TPS,TLS,SIJ,em_CAPE,em_TVP,em_pbase,em_bbase,
     .  em_dtvpdt1,em_dtvpdq1,em_dplcldt,em_dplcldr, ! sbl
     .  em_d_t2,em_d_q2,em_d_u2,em_d_v2,em_wd,em_qcond,em_qcondc)!sbl
c     print*,'apres convect '
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c -- sb: Appel schema statistique de nuages couple a la convection
c (Bony et Emanuel 2001):

c -- creer cvthermo.h qui contiendra les cstes thermo de LMDZ:

        iflagcon = 3
c       CALL cv_thermo(iflagcon)

c -- appel schema de nuages:

c       CALL CLOUDS_SUB_LS(klev,em_q,em_qs,em_t
c    i          ,em_p,em_ph,dtime,em_qcondc
c    o          ,em_cldf,em_cldq,em_pradj,em_ftadj,em_fradj,em_ifc)

        do k = 1, klev 
         cldf(i,k)  = em_cldf(k)  ! cloud fraction (0-1)
         cldq(i,k)  = em_cldq(k)  ! in-cloud water content (kg/kg)
         ftadj(i,k) = em_ftadj(k) ! (dT/dt)_{LS adj} (K/s)
         fqadj(i,k) = em_fradj(k) ! (dq/dt)_{LS adj} (kg/kg/s)
         ifc(i,k)   = em_ifc(k)   ! flag convergence clouds_gno (1 ou 2)
        enddo
        pradj(i) = em_pradj       ! precip from LS supersat adj (mm/day)

c sb --
c
c SB:
      if (iflag.ne.1 .and. iflag.ne.4) then
         em_CAPE = 0.
      do l = 1, klev
         em_upwd(l) = 0.
         em_dnwd(l) = 0.
         em_dnwdbis(l) = 0.
         emMa(l) = 0.
         em_TVP(l) = 0.
      enddo
      endif
c fin SB
c
c  If sig has been set to zero, then set Ma to zero
c
      sigsum = 0.
      do k = 1,klev
        sigsum = sigsum + em_work1(k)
      enddo
      if (sigsum .eq. 0.0) then
        do k = 1,klev
          emMa(k) = 0.
        enddo
      endif
c
c sb3d       print*,'i, iflag=',i,iflag
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       SORTIE DES ICB ET INB
c       en fait inb et icb correspondent au niveau ou se trouve
c       le nuage,le numero d'interface
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c modif SB:
      if (iflag.EQ.1 .or. iflag.EQ.4) then
       top=em_top
       bas=em_bas
       kbas(i) = em_bas
       ktop(i) = em_top
      endif
 
      pbase(i) = em_pbase
      bbase(i) = em_bbase
      rain(i) = em_precip/ 86400.0
      snow(i) = 0.0
      cape(i) = em_CAPE
      wd(i) = em_wd
      rflag(i) = float(iflag)
c SB      kbas(i) = em_bas
c SB      ktop(i) = em_top
      dplcldt(i) = em_dplcldt
      dplcldr(i) = em_dplcldr
      DO l = 1, klev
         d_t2(i,l) = dtime * em_d_t2(l) 
         d_q2(i,l) = dtime * em_d_q2(l)
         d_u2(i,l) = dtime * em_d_u2(l)
         d_v2(i,l) = dtime * em_d_v2(l)

         d_t(i,l) = dtime * em_d_t(l) 
         d_q(i,l) = dtime * em_d_q(l)
         d_u(i,l) = dtime * em_d_u(l)
         d_v(i,l) = dtime * em_d_v(l)
         do itra = 1, ntra
         d_tra(i,l,itra) = dtime * em_d_tra(l,itra)
         enddo
         upwd(i,l) = em_upwd(l)
         dnwd(i,l) = em_dnwd(l)
         dnwdbis(i,l) = em_dnwdbis(l)
         work1(i,l) = em_work1(l)
         work2(i,l) = em_work2(l)
         Ma(i,l)=emMa(l)
         tvp(i,l)=em_TVP(l)
         dtvpdt1(i,l) = em_dtvpdt1(l)
         dtvpdq1(i,l) = em_dtvpdq1(l)

         if (iflag_clw.eq.0) then
            qcond_incld(i,l) = em_qcondc(l)
         else if (iflag_clw.eq.1) then
            qcond_incld(i,l) = em_qcond(l)
         endif
      ENDDO
  999 CONTINUE

c   On calcule une eau liquide diagnostique en fonction de la 
c  precip.
      if ( iflag_clw.eq.2 ) then
      do l=1,klev
         do i=1,klon
            if (ktop(i)-kbas(i).gt.0.and.
     s         l.ge.kbas(i).and.l.le.ktop(i)) then
               qcond_incld(i,l)=rain(i)*8.e4
     s         /(pplay(i,kbas(i))-pplay(i,ktop(i)))
c    s         **2
            else
               qcond_incld(i,l)=0.
            endif
         enddo
         print*,'l=',l,',   qcond_incld=',qcond_incld(1,l)
      enddo
      endif
 

      RETURN
      END

