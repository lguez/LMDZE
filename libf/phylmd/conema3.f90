SUBROUTINE conema3 (dtime,paprs,pplay,t,q,u,v,tra,ntra, &
     work1,work2,d_t,d_q,d_u,d_v,d_tra, &
     rain, snow, kbas, ktop, &
     upwd,dnwd,dnwdbis,bas,top,Ma,cape,tvp,rflag, &
     pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr, &
     qcond_incld)

  ! From LMDZ4/libf/phylmd/conema3.F,v 1.1.1.1 2004/05/19 12:53:09

  use dimens_m
  use dimphy
  use SUPHEC_M
  use conema3_m
  use yoethf_m
  use fcttre
  IMPLICIT none
  !======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: schema de convection de Emanuel (1991) interface
  ! Mai 1998: Interface modifiee pour implementation dans LMDZ
  !======================================================================
  ! Arguments:
  ! dtime---input-R-pas d'integration (s)
  ! paprs---input-R-pression inter-couches (Pa)
  ! pplay---input-R-pression au milieu des couches (Pa)
  ! t-------input-R-temperature (K)
  ! q-------input-R-humidite specifique (kg/kg)
  ! u-------input-R-vitesse du vent zonal (m/s)
  ! v-------input-R-vitesse duvent meridien (m/s)
  ! tra-----input-R-tableau de rapport de melange des traceurs
  ! work*: input et output: deux variables de travail,
  !                            on peut les mettre a 0 au debut
  !
  ! d_t-----output-R-increment de la temperature
  ! d_q-----output-R-increment de la vapeur d'eau
  ! d_u-----output-R-increment de la vitesse zonale
  ! d_v-----output-R-increment de la vitesse meridienne
  ! d_tra---output-R-increment du contenu en traceurs
  ! rain----output-R-la pluie (mm/s)
  ! snow----output-R-la neige (mm/s)
  ! kbas----output-R-bas du nuage (integer)
  ! ktop----output-R-haut du nuage (integer)
  ! upwd----output-R-saturated updraft mass flux (kg/m**2/s)
  ! dnwd----output-R-saturated downdraft mass flux (kg/m**2/s)
  ! dnwdbis-output-R-unsaturated downdraft mass flux (kg/m**2/s)
  ! bas-----output-R-bas du nuage (real)
  ! top-----output-R-haut du nuage (real)
  ! Ma------output-R-flux ascendant non dilue (kg/m**2/s)
  ! cape----output-R-CAPE
  ! tvp-----output-R-virtual temperature of the lifted parcel
  ! rflag---output-R-flag sur le fonctionnement de convect
  ! pbase---output-R-pression a la base du nuage (Pa)
  ! bbase---output-R-buoyancy a la base du nuage (K)
  ! dtvpdt1-output-R-derivative of parcel virtual temp wrt T1 
  ! dtvpdq1-output-R-derivative of parcel virtual temp wrt Q1 
  ! dplcldt-output-R-derivative of the PCP pressure wrt T1
  ! dplcldr-output-R-derivative of the PCP pressure wrt Q1
  !======================================================================
  !
  INTEGER i, l,m,itra
  INTEGER ntra,ntrac !number of tracers; if no tracer transport
  ! is needed, set ntra = 1 (or 0)
  PARAMETER (ntrac=nqmx-2)
  REAL, intent(in):: dtime
  !
  REAL d_t2(klon,klev), d_q2(klon,klev) ! sbl
  REAL d_u2(klon,klev), d_v2(klon,klev) ! sbl
  REAL em_d_t2(klev), em_d_q2(klev)     ! sbl   
  REAL em_d_u2(klev), em_d_v2(klev)     ! sbl   
  ! 
  REAL, intent(in):: paprs(klon,klev+1)
  real, intent(in):: pplay(klon,klev)
  REAL, intent(in):: t(klon,klev)
  real q(klon,klev), d_t(klon,klev), d_q(klon,klev)
  REAL u(klon,klev), v(klon,klev)
  real, intent(in):: tra(klon,klev,ntra)
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
  !
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
  !
  REAL zx_t, zx_qs, zdelta, zcor
  INTEGER iflag
  REAL sigsum
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     VARIABLES A SORTIR
  !ccccccccccccccccccccccccccccccccccccccccccccccccc

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

  ! -- sb: pour schema nuages:

  integer iflagcon
  integer em_ifc(klev)

  real em_pradj
  real em_cldf(klev), em_cldq(klev)
  real em_ftadj(klev), em_fradj(klev)

  integer ifc(klon,klev)
  real pradj(klon)
  real cldf(klon,klev), cldq(klon,klev)
  real ftadj(klon,klev), fqadj(klon,klev)

  ! sb --

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !

  qcond_incld(:,:) = 0.
  !
  !$$$      print*,'debut conema'

  DO i = 1, klon
     DO l = 1, klev+1
        em_ph(l) = paprs(i,l) / 100.0
     ENDDO
     !
     DO l = 1, klev
        em_p(l) = pplay(i,l) / 100.0
        em_t(l) = t(i,l)
        em_q(l) = q(i,l)
        em_u(l) = u(i,l)
        em_v(l) = v(i,l)
        do itra = 1, ntra
           em_tra(l,itra) = tra(i,l,itra)
        enddo
        !$$$      print*,'em_t',em_t
        !$$$      print*,'em_q',em_q
        !$$$      print*,'em_qs',em_qs
        !$$$      print*,'em_u',em_u
        !$$$      print*,'em_v',em_v
        !$$$      print*,'em_tra',em_tra
        !$$$      print*,'em_p',em_p


        !
        zx_t = em_t(l)
        zdelta=MAX(0.,SIGN(1.,rtt-zx_t))
        zx_qs= r2es * FOEEW(zx_t,zdelta)/em_p(l)/100.0
        zx_qs=MIN(0.5,zx_qs)
        !$$$       print*,'zx_qs',zx_qs
        zcor=1./(1.-retv*zx_qs) 
        zx_qs=zx_qs*zcor
        em_qs(l) = zx_qs
        !$$$      print*,'em_qs',em_qs
        !
        em_work1(l) = work1(i,l)
        em_work2(l) = work2(i,l)
        emMke(l)=0
        !        emMa(l)=0
        !        Ma(i,l)=0

        em_dtvpdt1(l) = 0.
        em_dtvpdq1(l) = 0.
        dtvpdt1(i,l) = 0.
        dtvpdq1(i,l) = 0.
     ENDDO
     !
     em_dplcldt = 0.
     em_dplcldr = 0.
     rain(i) = 0.0
     snow(i) = 0.0
     kbas(i) = 1
     ktop(i) = 1
     ! ajout SB:
     bas = 1
     top = 1


     ! sb3d      write(*,1792) (em_work1(m),m=1,klev)
1792 format('sig avant convect ',/,10(1X,E13.5))
     !
     ! sb d      write(*,1793) (em_work2(m),m=1,klev)
1793 format('w avant convect ',/,10(1X,E13.5))

     !$$$      print*,'avant convect' 
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! 

     !     print*,'avant convect i=',i
     CALL convect3(dtime,epmax,ok_adj_ema, &
          em_t, em_q, em_qs,em_u ,em_v , &
          em_tra, em_p, em_ph, &
          klev, klev+1, klev-1,ntra, dtime, iflag, &
          em_d_t, em_d_q,em_d_u,em_d_v, &
          em_d_tra, em_precip, &
          em_bas, em_top,em_upwd, em_dnwd, em_dnwdbis, &
          em_work1, em_work2,emmip,emMke,emMa,Ment, &
          Qent,TPS,TLS,SIJ,em_CAPE,em_TVP,em_pbase,em_bbase, &
          em_dtvpdt1,em_dtvpdq1,em_dplcldt,em_dplcldr, &
          em_d_t2,em_d_q2,em_d_u2,em_d_v2,em_wd,em_qcond,em_qcondc)!sbl
     !     print*,'apres convect '
     !
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !
     ! -- sb: Appel schema statistique de nuages couple a la convection
     ! (Bony et Emanuel 2001):

     ! -- creer cvthermo.h qui contiendra les cstes thermo de LMDZ:

     iflagcon = 3
     !       CALL cv_thermo(iflagcon)

     ! -- appel schema de nuages:

     do k = 1, klev 
        cldf(i,k)  = em_cldf(k)  ! cloud fraction (0-1)
        cldq(i,k)  = em_cldq(k)  ! in-cloud water content (kg/kg)
        ftadj(i,k) = em_ftadj(k) ! (dT/dt)_{LS adj} (K/s)
        fqadj(i,k) = em_fradj(k) ! (dq/dt)_{LS adj} (kg/kg/s)
        ifc(i,k)   = em_ifc(k)   ! flag convergence clouds_gno (1 ou 2)
     enddo
     pradj(i) = em_pradj       ! precip from LS supersat adj (mm/day)

     ! sb --
     !
     ! SB:
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
     ! fin SB
     !
     !  If sig has been set to zero, then set Ma to zero
     !
     sigsum = 0.
     do k = 1,klev
        sigsum = sigsum + em_work1(k)
     enddo
     if (sigsum .eq. 0.0) then
        do k = 1,klev
           emMa(k) = 0.
        enddo
     endif
     !
     ! sb3d       print*,'i, iflag=',i,iflag
     ! 
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !
     !       SORTIE DES ICB ET INB
     !       en fait inb et icb correspondent au niveau ou se trouve
     !       le nuage,le numero d'interface
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! modif SB:
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
     ! SB      kbas(i) = em_bas
     ! SB      ktop(i) = em_top
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
  end DO

  !   On calcule une eau liquide diagnostique en fonction de la 
  !  precip.
  if ( iflag_clw.eq.2 ) then
     do l=1,klev
        do i=1,klon
           if (ktop(i)-kbas(i).gt.0.and. &
                l.ge.kbas(i).and.l.le.ktop(i)) then
              qcond_incld(i,l)=rain(i)*8.e4 &
                   /(pplay(i,kbas(i))-pplay(i,ktop(i)))
              !    s         **2
           else
              qcond_incld(i,l)=0.
           endif
        enddo
        print*,'l=',l,',   qcond_incld=',qcond_incld(1,l)
     enddo
  endif

END SUBROUTINE conema3
