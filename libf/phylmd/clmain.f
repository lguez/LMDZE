      SUBROUTINE clmain(dtime,itap,date0,pctsrf,pctsrf_new,
     .                  t,q,u,v,
     .                  jour, rmu0, co2_ppm,
     .                  ok_veget, ocean, npas, nexca, ts,
     .                  soil_model,cdmmax, cdhmax,
     .                  ksta, ksta_ter, ok_kzmin, ftsoil,qsol,
     .                  paprs,pplay,snow,qsurf,evap,albe,alblw,
     .                  fluxlat,
     .                  rain_f, snow_f, solsw, sollw, sollwdown, fder,
     .                  rlon, rlat, cufi, cvfi, rugos,
     .                  debut, lafin, agesno,rugoro,
     .                  d_t,d_q,d_u,d_v,d_ts,
     .                  flux_t,flux_q,flux_u,flux_v,cdragh,cdragm,
     .                  q2,
     .                  dflux_t,dflux_q,
     .                  zcoefh,zu1,zv1, t2m, q2m, u10m, v10m,
cIM cf. AM : pbl
     .                  pblh,capCL,oliqCL,cteiCL,pblT,
     .                  therm,trmb1,trmb2,trmb3,plcl,
     .                  fqcalving,ffonte, run_off_lic_0,
cIM "slab" ocean
     .                  flux_o, flux_g, tslab, seaice)

!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/clmain.F,v 1.6 2005/11/16 14:47:19 lmdzadmin Exp $
!
c
c
cAA REM:
cAA-----
cAA Tout ce qui a trait au traceurs est dans phytrac maintenant
cAA pour l'instant le calcul de la couche limite pour les traceurs
cAA se fait avec cltrac et ne tient pas compte de la differentiation
cAA des sous-fraction de sol.
cAA REM bis :
cAA----------
cAA Pour pouvoir extraire les coefficient d'echanges et le vent 
cAA dans la premiere couche, 3 champs supplementaires ont ete crees
cAA zcoefh,zu1 et zv1. Pour l'instant nous avons moyenne les valeurs
cAA de ces trois champs sur les 4 subsurfaces du modele. Dans l'avenir 
cAA si les informations des subsurfaces doivent etre prises en compte
cAA il faudra sortir ces memes champs en leur ajoutant une dimension, 
cAA c'est a dire nbsrf (nbre de subsurface).
      USE ioipsl
      USE interface_surf
      use dimens_m
      use indicesol
      use dimphy
      use dimsoil
      use temps
      use iniprint
      use YOMCST
      use yoethf
      use fcttre
      use conf_phys_m
      use gath_cpl, only: gath2cpl
      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c Objet: interface de "couche limite" (diffusion verticale)
c Arguments:
c dtime----input-R- interval du temps (secondes)
c itap-----input-I- numero du pas de temps
c date0----input-R- jour initial
c t--------input-R- temperature (K)
c q--------input-R- vapeur d'eau (kg/kg)
c u--------input-R- vitesse u
c v--------input-R- vitesse v
c ts-------input-R- temperature du sol (en Kelvin)
c paprs----input-R- pression a intercouche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c radsol---input-R- flux radiatif net (positif vers le sol) en W/m**2
c rlat-----input-R- latitude en degree
c rugos----input-R- longeur de rugosite (en m)
c cufi-----input-R- resolution des mailles en x (m)
c cvfi-----input-R- resolution des mailles en y (m)
c
c d_t------output-R- le changement pour "t"
c d_q------output-R- le changement pour "q"
c d_u------output-R- le changement pour "u"
c d_v------output-R- le changement pour "v"
c d_ts-----output-R- le changement pour "ts"
c flux_t---output-R- flux de chaleur sensible (CpT) J/m**2/s (W/m**2)
c                    (orientation positive vers le bas)
c flux_q---output-R- flux de vapeur d'eau (kg/m**2/s)
c flux_u---output-R- tension du vent X: (kg m/s)/(m**2 s) ou Pascal
c flux_v---output-R- tension du vent Y: (kg m/s)/(m**2 s) ou Pascal
c dflux_t derive du flux sensible
c dflux_q derive du flux latent
cIM "slab" ocean
c flux_g---output-R-  flux glace (pour OCEAN='slab  ')
c flux_o---output-R-  flux ocean (pour OCEAN='slab  ')
c tslab-in/output-R temperature du slab ocean (en Kelvin) ! uniqmnt pour slab
c seaice---output-R-  glace de mer (kg/m2) (pour OCEAN='slab  ')
ccc
c ffonte----Flux thermique utilise pour fondre la neige
c fqcalving-Flux d'eau "perdue" par la surface et necessaire pour limiter la
c           hauteur de neige, en kg/m2/s
cAA on rajoute en output yu1 et yv1 qui sont les vents dans 
cAA la premiere couche
cAA ces 4 variables sont maintenant traites dans phytrac
c itr--------input-I- nombre de traceurs
c tr---------input-R- q. de traceurs
c flux_surf--input-R- flux de traceurs a la surface
c d_tr-------output-R tendance de traceurs
cIM cf. AM : PBL
c trmb1-------deep_cape
c trmb2--------inhibition 
c trmb3-------Point Omega
c Cape(klon)-------Cape du thermique
c EauLiq(klon)-------Eau liqu integr du thermique
c ctei(klon)-------Critere d'instab d'entrainmt des nuages de CL
c lcl------- Niveau de condensation
c pblh------- HCL
c pblT------- T au nveau HCL
c======================================================================
c$$$ PB ajout pour soil
c
      REAL, intent(in):: dtime
      real date0
      integer, intent(in):: itap
      REAL t(klon,klev), q(klon,klev)
      REAL u(klon,klev), v(klon,klev)
cIM 230604 BAD  REAL radsol(klon) ??? 
      REAL, intent(in):: paprs(klon,klev+1)
      real, intent(in):: pplay(klon,klev)
      REAL, intent(in):: rlon(klon), rlat(klon)
      real cufi(klon), cvfi(klon)
      REAL d_t(klon, klev), d_q(klon, klev)
      REAL d_u(klon, klev), d_v(klon, klev)
      REAL flux_t(klon,klev, nbsrf), flux_q(klon,klev, nbsrf)
      REAL dflux_t(klon), dflux_q(klon)
cIM "slab" ocean
      REAL flux_o(klon), flux_g(klon)
      REAL y_flux_o(klon), y_flux_g(klon)
      REAL tslab(klon), ytslab(klon)
      REAL seaice(klon), y_seaice(klon)
cIM cf JLD
      REAL y_fqcalving(klon), y_ffonte(klon)
      REAL fqcalving(klon,nbsrf), ffonte(klon,nbsrf)
      REAL run_off_lic_0(klon), y_run_off_lic_0(klon)

      REAL flux_u(klon,klev, nbsrf), flux_v(klon,klev, nbsrf)
      REAL rugmer(klon), agesno(klon,nbsrf),rugoro(klon)
      REAL cdragh(klon), cdragm(klon)
      integer jour            ! jour de l'annee en cours
      real rmu0(klon)         ! cosinus de l'angle solaire zenithal
      REAL co2_ppm            ! taux CO2 atmosphere
      LOGICAL, intent(in):: debut
      logical, intent(in):: lafin
      logical ok_veget
      character(len=*), intent(IN):: ocean
      integer npas, nexca
c
      REAL pctsrf(klon,nbsrf)
      REAL ts(klon,nbsrf)
      REAL d_ts(klon,nbsrf)
      REAL snow(klon,nbsrf)
      REAL qsurf(klon,nbsrf)
      REAL evap(klon,nbsrf)
      REAL albe(klon,nbsrf)
      REAL alblw(klon,nbsrf)
c$$$ PB
      REAL fluxlat(klon,nbsrf)
C
      real rain_f(klon), snow_f(klon)
      REAL fder(klon)
cIM cf. JLD   REAL sollw(klon), solsw(klon), sollwdown(klon)
      REAL sollw(klon,nbsrf), solsw(klon,nbsrf), sollwdown(klon)
      REAL rugos(klon,nbsrf)
C la nouvelle repartition des surfaces sortie de l'interface
      REAL pctsrf_new(klon,nbsrf)
cAA
      REAL zcoefh(klon,klev)
      REAL zu1(klon)
      REAL zv1(klon)
cAA
c$$$ PB ajout pour soil
      LOGICAL, intent(in):: soil_model
cIM ajout seuils cdrm, cdrh
      REAL cdmmax, cdhmax
cIM: 261103
      REAL ksta, ksta_ter
      LOGICAL ok_kzmin
cIM: 261103
      REAL ftsoil(klon,nsoilmx,nbsrf)
      REAL ytsoil(klon,nsoilmx)
      REAL qsol(klon)
c======================================================================
      EXTERNAL clqh, clvent, coefkz, calbeta, cltrac
c======================================================================
      REAL yts(klon), yrugos(klon), ypct(klon), yz0_new(klon)
      REAL yalb(klon)
      REAL yalblw(klon)
      REAL yu1(klon), yv1(klon)
      real ysnow(klon), yqsurf(klon), yagesno(klon), yqsol(klon)
      real yrain_f(klon), ysnow_f(klon)
      real ysollw(klon), ysolsw(klon), ysollwdown(klon)
      real yfder(klon), ytaux(klon), ytauy(klon)
      REAL yrugm(klon), yrads(klon),yrugoro(klon)
c$$$ PB
      REAL yfluxlat(klon)
C
      REAL y_d_ts(klon)
      REAL y_d_t(klon, klev), y_d_q(klon, klev)
      REAL y_d_u(klon, klev), y_d_v(klon, klev)
      REAL y_flux_t(klon,klev), y_flux_q(klon,klev)
      REAL y_flux_u(klon,klev), y_flux_v(klon,klev)
      REAL y_dflux_t(klon), y_dflux_q(klon)
      REAL ycoefh(klon,klev), ycoefm(klon,klev)
      REAL yu(klon,klev), yv(klon,klev)
      REAL yt(klon,klev), yq(klon,klev)
      REAL ypaprs(klon,klev+1), ypplay(klon,klev), ydelp(klon,klev)
c
      LOGICAL ok_nonloc
      PARAMETER (ok_nonloc=.FALSE.)
      REAL ycoefm0(klon,klev), ycoefh0(klon,klev)

cIM 081204 hcl_Anne ? BEG
      real yzlay(klon,klev),yzlev(klon,klev+1),yteta(klon,klev)
      real ykmm(klon,klev+1),ykmn(klon,klev+1)
      real ykmq(klon,klev+1)
      real yq2(klon,klev+1),q2(klon,klev+1,nbsrf)
      real q2diag(klon,klev+1)
cIM 081204   real yustar(klon),y_cd_m(klon),y_cd_h(klon)
cIM 081204 hcl_Anne ? END
c
      REAL u1lay(klon), v1lay(klon)
      REAL delp(klon,klev)
      INTEGER i, k, nsrf 
cAA   INTEGER it
      INTEGER ni(klon), knon, j
c Introduction d'une variable "pourcentage potentiel" pour tenir compte
c des eventuelles apparitions et/ou disparitions de la glace de mer
      REAL pctsrf_pot(klon,nbsrf)
      
c======================================================================
      REAL zx_alf1, zx_alf2 !valeur ambiante par extrapola.
c======================================================================
c
c maf pour sorties IOISPL en cas de debugagage
c
      CHARACTER*80 cldebug
      SAVE cldebug
      CHARACTER*8 cl_surf(nbsrf)
      SAVE cl_surf
      INTEGER nhoridbg, nidbg
      SAVE nhoridbg, nidbg
      INTEGER ndexbg(iim*(jjm+1))
      REAL zx_lon(iim,jjm+1), zx_lat(iim,jjm+1), zjulian
      REAL tabindx(klon)
      REAL debugtab(iim,jjm+1)
      LOGICAL first_appel
      SAVE first_appel
      DATA first_appel/.true./
      LOGICAL debugindex
      SAVE debugindex
      DATA debugindex/.false./
      integer idayref
      REAL t2m(klon,nbsrf), q2m(klon,nbsrf)
      REAL u10m(klon,nbsrf), v10m(klon,nbsrf)
c
      REAL yt2m(klon), yq2m(klon), yu10m(klon)
      REAL yustar(klon)
c -- LOOP
       REAL yu10mx(klon)
       REAL yu10my(klon)
       REAL ywindsp(klon)
c -- LOOP
c
      REAL yt10m(klon), yq10m(klon)
cIM cf. AM : pbl, hbtm2 (Comme les autres diagnostics on cumule ds physic ce qui 
c   permet de sortir les grdeurs par sous surface)
      REAL pblh(klon,nbsrf)
      REAL plcl(klon,nbsrf)
      REAL capCL(klon,nbsrf)
      REAL oliqCL(klon,nbsrf)
      REAL cteiCL(klon,nbsrf)
      REAL pblT(klon,nbsrf)
      REAL therm(klon,nbsrf)
      REAL trmb1(klon,nbsrf)
      REAL trmb2(klon,nbsrf)
      REAL trmb3(klon,nbsrf)
      REAL ypblh(klon)
      REAL ylcl(klon)
      REAL ycapCL(klon)
      REAL yoliqCL(klon)
      REAL ycteiCL(klon)
      REAL ypblT(klon)
      REAL ytherm(klon)
      REAL ytrmb1(klon)
      REAL ytrmb2(klon)
      REAL ytrmb3(klon)
      REAL y_cd_h(klon), y_cd_m(klon)
c     REAL ygamt(klon,2:klev) ! contre-gradient pour temperature
c     REAL ygamq(klon,2:klev) ! contre-gradient pour humidite
      REAL uzon(klon), vmer(klon)
      REAL tair1(klon), qair1(klon), tairsol(klon)
      REAL psfce(klon), patm(klon)
c
      REAL qairsol(klon), zgeo1(klon)
      REAL rugo1(klon)
c
      LOGICAL zxli ! utiliser un jeu de fonctions simples
      PARAMETER (zxli=.FALSE.)
c
      REAL zt, zqs, zdelta, zcor
      REAL t_coup
      PARAMETER(t_coup=273.15)
C
      character (len = 20) :: modname = 'clmain'
      LOGICAL check
      PARAMETER (check=.false.)


c initialisation Anne
      ytherm(:) = 0.
C
      if (check) THEN
          write(*,*) modname,'  klon=',klon
CC        call flush(6)
      endif
      IF (debugindex .and. first_appel) THEN
          first_appel=.false.
!
! initialisation sorties netcdf
!
          idayref = day_ini
          CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
          CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,zx_lon)
          DO i = 1, iim
            zx_lon(i,1) = rlon(i+1)
            zx_lon(i,jjm+1) = rlon(i+1)
          ENDDO
          CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,zx_lat)
          cldebug='sous_index'
          CALL histbeg_totreg(cldebug, iim,zx_lon(:,1),jjm+1,
     $        zx_lat(1,:),1,iim,1,jjm
     $        +1, itau_phy,zjulian,dtime,nhoridbg,nidbg) 
! no vertical axis
          cl_surf(1)='ter'
          cl_surf(2)='lic'
          cl_surf(3)='oce'
          cl_surf(4)='sic'
          DO nsrf=1,nbsrf
            CALL histdef(nidbg, cl_surf(nsrf),cl_surf(nsrf), "-",iim,
     $          jjm+1,nhoridbg, 1, 1, 1, -99, 32, "inst", dtime,dtime) 
          END DO
          CALL histend(nidbg)
          CALL histsync(nidbg)
      ENDIF 
          
      DO k = 1, klev   ! epaisseur de couche
      DO i = 1, klon
         delp(i,k) = paprs(i,k)-paprs(i,k+1)
      ENDDO
      ENDDO
      DO i = 1, klon  ! vent de la premiere couche
         zx_alf1 = 1.0
         zx_alf2 = 1.0 - zx_alf1
         u1lay(i) = u(i,1)*zx_alf1 + u(i,2)*zx_alf2
         v1lay(i) = v(i,1)*zx_alf1 + v(i,2)*zx_alf2
      ENDDO
c
c initialisation:
c
      DO i = 1, klon
         rugmer(i) = 0.0
         cdragh(i) = 0.0
         cdragm(i) = 0.0
         dflux_t(i) = 0.0
         dflux_q(i) = 0.0
         zu1(i) = 0.0
         zv1(i) = 0.0
      ENDDO
      ypct = 0.0
      yts = 0.0
      ysnow = 0.0
      yqsurf = 0.0
      yalb = 0.0
      yalblw = 0.0
      yrain_f = 0.0
      ysnow_f = 0.0
      yfder = 0.0
      ytaux = 0.0
      ytauy = 0.0
      ysolsw = 0.0
      ysollw = 0.0
      ysollwdown = 0.0
      yrugos = 0.0
      yu1 = 0.0
      yv1 = 0.0
      yrads = 0.0
      ypaprs = 0.0
      ypplay = 0.0
      ydelp = 0.0
      yu = 0.0
      yv = 0.0
      yt = 0.0
      yq = 0.0
      pctsrf_new = 0.0
      y_flux_u = 0.0
      y_flux_v = 0.0
C$$ PB
      y_dflux_t = 0.0
      y_dflux_q = 0.0
      ytsoil = 999999.
      yrugoro = 0.
c -- LOOP
      yu10mx = 0.0
      yu10my = 0.0
      ywindsp = 0.0
c -- LOOP
      DO nsrf = 1, nbsrf
      DO i = 1, klon
         d_ts(i,nsrf) = 0.0
      ENDDO
      END DO
C§§§ PB
      yfluxlat=0.
      flux_t = 0.
      flux_q = 0.
      flux_u = 0.
      flux_v = 0.
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
c$$$         flux_t(i,k) = 0.0
c$$$         flux_q(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
c$$$         flux_u(i,k) = 0.0
c$$$         flux_v(i,k) = 0.0
         zcoefh(i,k) = 0.0
      ENDDO
      ENDDO
cAA      IF (itr.GE.1) THEN
cAA      DO it = 1, itr
cAA      DO k = 1, klev
cAA      DO i = 1, klon
cAA         d_tr(i,k,it) = 0.0
cAA      ENDDO
cAA      ENDDO
cAA      ENDDO
cAA      ENDIF

c
c Boucler sur toutes les sous-fractions du sol:
c
C Initialisation des "pourcentages potentiels". On considere ici qu'on 
C peut avoir potentiellementdela glace sur tout le domaine oceanique 
C (a affiner)

      pctsrf_pot = pctsrf
      pctsrf_pot(:,is_oce) = 1. - zmasq(:)
      pctsrf_pot(:,is_sic) = 1. - zmasq(:)

      DO 99999 nsrf = 1, nbsrf

c chercher les indices:
      DO j = 1, klon
         ni(j) = 0
      ENDDO
      knon = 0
      DO i = 1, klon

C pour determiner le domaine a traiter on utilise les surfaces "potentielles"
C  
      IF (pctsrf_pot(i,nsrf).GT.epsfra) THEN
         knon = knon + 1
         ni(knon) = i
      ENDIF
      ENDDO
c
      if (check) THEN
          write(*,*)'CLMAIN, nsrf, knon =',nsrf, knon
CC        call flush(6)
      endif
c
c variables pour avoir une sortie IOIPSL des INDEX
c
      IF (debugindex) THEN 
          tabindx(:)=0.
c          tabindx(1:knon)=(/FLOAT(i),i=1:knon/)
          DO i=1,knon
            tabindx(1:knon)=FLOAT(i)
          END DO 
          debugtab(:,:)=0.
          ndexbg(:)=0
          CALL gath2cpl(tabindx,debugtab,klon,knon,iim,jjm,ni)
          CALL histwrite(nidbg,cl_surf(nsrf),itap,debugtab,iim*(jjm+1)
     $        ,ndexbg)
      ENDIF 
      IF (knon.EQ.0) GOTO 99999
      DO j = 1, knon
      i = ni(j)
        ypct(j) = pctsrf(i,nsrf)
        yts(j) = ts(i,nsrf)
cIM "slab" ocean
c        PRINT *, 'tslab = ', i, tslab(i)
        ytslab(i) = tslab(i)
c
        ysnow(j) = snow(i,nsrf)
        yqsurf(j) = qsurf(i,nsrf)
        yalb(j) = albe(i,nsrf)
        yalblw(j) = alblw(i,nsrf)
        yrain_f(j) = rain_f(i)
        ysnow_f(j) = snow_f(i)
        yagesno(j) = agesno(i,nsrf)
        yfder(j) = fder(i)
        ytaux(j) = flux_u(i,1,nsrf)
        ytauy(j) = flux_v(i,1,nsrf)
        ysolsw(j) = solsw(i,nsrf)
        ysollw(j) = sollw(i,nsrf)
        ysollwdown(j) = sollwdown(i)
        yrugos(j) = rugos(i,nsrf)
        yrugoro(j) = rugoro(i)
        yu1(j) = u1lay(i)
        yv1(j) = v1lay(i)
        yrads(j) =  ysolsw(j)+ ysollw(j)
        ypaprs(j,klev+1) = paprs(i,klev+1)
        y_run_off_lic_0(j) = run_off_lic_0(i)
c -- LOOP
       yu10mx(j) = u10m(i,nsrf)
       yu10my(j) = v10m(i,nsrf)
       ywindsp(j) = SQRT(yu10mx(j)*yu10mx(j) + yu10my(j)*yu10my(j) )
c -- LOOP
      END DO
C
C     IF bucket model for continent, copy soil water content
      IF ( nsrf .eq. is_ter .and. .not. ok_veget ) THEN 
          DO j = 1, knon
            i = ni(j)
            yqsol(j) = qsol(i)
          END DO
      ELSE 
          yqsol(:)=0.
      ENDIF 
c$$$ PB ajour pour soil
      DO k = 1, nsoilmx
        DO j = 1, knon
          i = ni(j)
          ytsoil(j,k) = ftsoil(i,k,nsrf)
        END DO  
      END DO 
      DO k = 1, klev
      DO j = 1, knon
      i = ni(j)
        ypaprs(j,k) = paprs(i,k)
        ypplay(j,k) = pplay(i,k)
        ydelp(j,k) = delp(i,k)
        yu(j,k) = u(i,k)
        yv(j,k) = v(i,k)
        yt(j,k) = t(i,k)
        yq(j,k) = q(i,k)
      ENDDO
      ENDDO
c
c
c calculer Cdrag et les coefficients d'echange
      CALL coefkz(nsrf, knon, ypaprs, ypplay,
cIM 261103
     .     ksta, ksta_ter,
cIM 261103
     .            yts, yrugos, yu, yv, yt, yq,
     .            yqsurf, 
     .            ycoefm, ycoefh)
cIM 081204 BEG
cCR test
      if (iflag_pbl.eq.1) then
cIM 081204 END
        CALL coefkz2(nsrf, knon, ypaprs, ypplay,yt,
     .                  ycoefm0, ycoefh0)
        DO k = 1, klev
        DO i = 1, knon
           ycoefm(i,k) = MAX(ycoefm(i,k),ycoefm0(i,k))
           ycoefh(i,k) = MAX(ycoefh(i,k),ycoefh0(i,k))
        ENDDO
        ENDDO
      endif
c
cIM cf JLD : on seuille ycoefm et ycoefh
      if (nsrf.eq.is_oce) then
         do j=1,knon
c           ycoefm(j,1)=min(ycoefm(j,1),1.1E-3)
            ycoefm(j,1)=min(ycoefm(j,1),cdmmax)
c           ycoefh(j,1)=min(ycoefh(j,1),1.1E-3)
            ycoefh(j,1)=min(ycoefh(j,1),cdhmax)
         enddo
      endif

c
cIM: 261103
      if (ok_kzmin) THEN
cIM cf FH: 201103 BEG
c   Calcul d'une diffusion minimale pour les conditions tres stables.
      call coefkzmin(knon,ypaprs,ypplay,yu,yv,yt,yq,ycoefm
     .   ,ycoefm0,ycoefh0)
c      call dump2d(iim,jjm-1,ycoefm(2:klon-1,2), 'KZ         ')
c      call dump2d(iim,jjm-1,ycoefm0(2:klon-1,2),'KZMIN      ')
 
       if ( 1.eq.1 ) then
       DO k = 1, klev
       DO i = 1, knon
          ycoefm(i,k) = MAX(ycoefm(i,k),ycoefm0(i,k))
          ycoefh(i,k) = MAX(ycoefh(i,k),ycoefh0(i,k))
       ENDDO
       ENDDO
       endif
cIM cf FH: 201103 END
      endif !ok_kzmin
cIM: 261103


      IF (iflag_pbl.ge.3) then

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c MELLOR ET YAMADA adapte a Mars Richard Fournier et Frederic Hourdin
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         yzlay(1:knon,1)=
     .   RD*yt(1:knon,1)/(0.5*(ypaprs(1:knon,1)+ypplay(1:knon,1)))
     .   *(ypaprs(1:knon,1)-ypplay(1:knon,1))/RG
         do k=2,klev
            yzlay(1:knon,k)=
     .      yzlay(1:knon,k-1)+RD*0.5*(yt(1:knon,k-1)+yt(1:knon,k))
     .      /ypaprs(1:knon,k)*(ypplay(1:knon,k-1)-ypplay(1:knon,k))/RG
         enddo
         do k=1,klev
            yteta(1:knon,k)=
     .      yt(1:knon,k)*(ypaprs(1:knon,1)/ypplay(1:knon,k))**rkappa
     .      *(1.+0.61*yq(1:knon,k))
         enddo
         yzlev(1:knon,1)=0.
         yzlev(1:knon,klev+1)=2.*yzlay(1:knon,klev)-yzlay(1:knon,klev-1)
         do k=2,klev
            yzlev(1:knon,k)=0.5*(yzlay(1:knon,k)+yzlay(1:knon,k-1))
         enddo
         DO k = 1, klev+1
            DO j = 1, knon
               i = ni(j)
               yq2(j,k)=q2(i,k,nsrf)
            enddo
         enddo


c   Bug introduit volontairement pour converger avec les resultats
c  du papier sur les thermiques.
         if (1.eq.1) then
         y_cd_m(1:knon) = ycoefm(1:knon,1)
         y_cd_h(1:knon) = ycoefh(1:knon,1)
         else
         y_cd_h(1:knon) = ycoefm(1:knon,1)
         y_cd_m(1:knon) = ycoefh(1:knon,1)
         endif
         call ustarhb(knon,yu,yv,y_cd_m, yustar)

        if (prt_level > 9) THEN
          print *,'USTAR = ',yustar
        ENDIF

c   iflag_pbl peut etre utilise comme longuer de melange

         if (iflag_pbl.ge.11) then
            call vdif_kcay(knon,dtime,rg,rd,ypaprs,yt
     s      ,yzlev,yzlay,yu,yv,yteta
     s      ,y_cd_m,yq2,q2diag,ykmm,ykmn,yustar,
     s      iflag_pbl)
         else
            call yamada4(knon,dtime,rg,rd,ypaprs,yt
     s      ,yzlev,yzlay,yu,yv,yteta
     s      ,y_cd_m,yq2,ykmm,ykmn,ykmq,yustar,
     s      iflag_pbl)
         endif

         ycoefm(1:knon,1)=y_cd_m(1:knon)
         ycoefh(1:knon,1)=y_cd_h(1:knon)
         ycoefm(1:knon,2:klev)=ykmm(1:knon,2:klev)
         ycoefh(1:knon,2:klev)=ykmn(1:knon,2:klev)


      ENDIF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c calculer la diffusion des vitesses "u" et "v"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CALL clvent(knon,dtime,yu1,yv1,ycoefm,yt,yu,ypaprs,ypplay,ydelp,
     s            y_d_u,y_flux_u)
      CALL clvent(knon,dtime,yu1,yv1,ycoefm,yt,yv,ypaprs,ypplay,ydelp,
     s            y_d_v,y_flux_v)

c pour le couplage
      ytaux = y_flux_u(:,1)
      ytauy = y_flux_v(:,1)

c FH modif sur le cdrag temperature
c$$$PB : déplace dans clcdrag
c$$$      do i=1,knon
c$$$         ycoefh(i,1)=ycoefm(i,1)*0.8
c$$$      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c calculer la diffusion de "q" et de "h"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL clqh(dtime, itap, date0,jour, debut,lafin,
     e          rlon, rlat, cufi, cvfi,
     e          knon, nsrf, ni, pctsrf,
     e          soil_model, ytsoil,yqsol,
     e          ok_veget, ocean, npas, nexca,
     e          rmu0, co2_ppm, yrugos, yrugoro,
     e          yu1, yv1, ycoefh,
     e          yt,yq,yts,ypaprs,ypplay,
     e          ydelp,yrads,yalb, yalblw, ysnow, yqsurf, 
     e          yrain_f, ysnow_f, yfder, ytaux, ytauy,
c -- LOOP
     e          ywindsp, 
c -- LOOP
c$$$     e          ysollw, ysolsw,
     e          ysollw, ysollwdown, ysolsw,yfluxlat,
     s          pctsrf_new, yagesno,
     s          y_d_t, y_d_q, y_d_ts, yz0_new,
     s          y_flux_t, y_flux_q, y_dflux_t, y_dflux_q,
     s          y_fqcalving,y_ffonte,y_run_off_lic_0,
cIM "slab" ocean
     s          y_flux_o, y_flux_g, ytslab, y_seaice)
c
c calculer la longueur de rugosite sur ocean
      yrugm=0.
      IF (nsrf.EQ.is_oce) THEN
      DO j = 1, knon
         yrugm(j) = 0.018*ycoefm(j,1) * (yu1(j)**2+yv1(j)**2)/RG 
     $      +  0.11*14e-6 / sqrt(ycoefm(j,1) * (yu1(j)**2+yv1(j)**2))
         yrugm(j) = MAX(1.5e-05,yrugm(j))
      ENDDO
      ENDIF
      DO j = 1, knon
         y_dflux_t(j) = y_dflux_t(j) * ypct(j)
         y_dflux_q(j) = y_dflux_q(j) * ypct(j)
         yu1(j) = yu1(j) *  ypct(j)
         yv1(j) = yv1(j) *  ypct(j)
      ENDDO
c
      DO k = 1, klev
        DO j = 1, knon
          i = ni(j)
          ycoefh(j,k) = ycoefh(j,k) * ypct(j)
          ycoefm(j,k) = ycoefm(j,k) * ypct(j)
          y_d_t(j,k) = y_d_t(j,k) * ypct(j)
          y_d_q(j,k) = y_d_q(j,k) * ypct(j)
C§§§ PB
          flux_t(i,k,nsrf) = y_flux_t(j,k)
          flux_q(i,k,nsrf) = y_flux_q(j,k)
          flux_u(i,k,nsrf) = y_flux_u(j,k)
          flux_v(i,k,nsrf) = y_flux_v(j,k)
c$$$ PB        y_flux_t(j,k) = y_flux_t(j,k) * ypct(j)
c$$$ PB        y_flux_q(j,k) = y_flux_q(j,k) * ypct(j)
          y_d_u(j,k) = y_d_u(j,k) * ypct(j)
          y_d_v(j,k) = y_d_v(j,k) * ypct(j)
c$$$ PB        y_flux_u(j,k) = y_flux_u(j,k) * ypct(j)
c$$$ PB        y_flux_v(j,k) = y_flux_v(j,k) * ypct(j)
        ENDDO
      ENDDO


      evap(:,nsrf) = - flux_q(:,1,nsrf)
c
      albe(:, nsrf) = 0.
      alblw(:, nsrf) = 0.
      snow(:, nsrf) = 0.
      qsurf(:, nsrf) = 0.
      rugos(:, nsrf) = 0.
      fluxlat(:,nsrf) = 0.
      DO j = 1, knon
         i = ni(j)
         d_ts(i,nsrf) = y_d_ts(j)
         albe(i,nsrf) = yalb(j)
         alblw(i,nsrf) = yalblw(j)
         snow(i,nsrf) = ysnow(j)
         qsurf(i,nsrf) = yqsurf(j)
         rugos(i,nsrf) = yz0_new(j)
         fluxlat(i,nsrf) = yfluxlat(j)
c$$$ pb         rugmer(i) = yrugm(j)
         IF (nsrf .EQ. is_oce) then 
           rugmer(i) = yrugm(j)
           rugos(i,nsrf) = yrugm(j)
         endif	
cIM cf JLD ??
         agesno(i,nsrf) = yagesno(j)
         fqcalving(i,nsrf) = y_fqcalving(j)        
         ffonte(i,nsrf) = y_ffonte(j)        
         cdragh(i) = cdragh(i) + ycoefh(j,1)
         cdragm(i) = cdragm(i) + ycoefm(j,1)
         dflux_t(i) = dflux_t(i) + y_dflux_t(j)
         dflux_q(i) = dflux_q(i) + y_dflux_q(j)
         zu1(i) = zu1(i) + yu1(j)
         zv1(i) = zv1(i) + yv1(j)
      END DO
      IF ( nsrf .eq. is_ter ) THEN 
      DO j = 1, knon
         i = ni(j)
         qsol(i) = yqsol(j)
      END DO
      END IF 
      IF ( nsrf .eq. is_lic ) THEN 
        DO j = 1, knon
          i = ni(j)
          run_off_lic_0(i) = y_run_off_lic_0(j)
        END DO
      END IF 
c$$$ PB ajout pour soil
      ftsoil(:,:,nsrf) = 0.
      DO k = 1, nsoilmx
        DO j = 1, knon
          i = ni(j)
          ftsoil(i, k, nsrf) = ytsoil(j,k)
        END DO 
      END DO 
c
      DO j = 1, knon
      i = ni(j)
      DO k = 1, klev
         d_t(i,k) = d_t(i,k) + y_d_t(j,k)
         d_q(i,k) = d_q(i,k) + y_d_q(j,k)
c$$$ PB        flux_t(i,k) = flux_t(i,k) + y_flux_t(j,k)
c$$$         flux_q(i,k) = flux_q(i,k) + y_flux_q(j,k)
         d_u(i,k) = d_u(i,k) + y_d_u(j,k)
         d_v(i,k) = d_v(i,k) + y_d_v(j,k)
c$$$  PB       flux_u(i,k) = flux_u(i,k) + y_flux_u(j,k)
c$$$         flux_v(i,k) = flux_v(i,k) + y_flux_v(j,k)
         zcoefh(i,k) = zcoefh(i,k) + ycoefh(j,k)
      ENDDO
      ENDDO
c
c
ccc diagnostic t,q a 2m et u, v a 10m
c
      DO j=1, knon
        i = ni(j)
        uzon(j) = yu(j,1) + y_d_u(j,1)
        vmer(j) = yv(j,1) + y_d_v(j,1)
        tair1(j) = yt(j,1) + y_d_t(j,1)
        qair1(j) = yq(j,1) + y_d_q(j,1)
        zgeo1(j) = RD * tair1(j) / (0.5*(ypaprs(j,1)+ypplay(j,1)))
     &                   * (ypaprs(j,1)-ypplay(j,1))
        tairsol(j) = yts(j) + y_d_ts(j)
        rugo1(j) = yrugos(j)
        IF(nsrf.EQ.is_oce) THEN
         rugo1(j) = rugos(i,nsrf)
        ENDIF
        psfce(j)=ypaprs(j,1)
        patm(j)=ypplay(j,1)
c
        qairsol(j) = yqsurf(j)
      ENDDO
c
      CALL stdlevvar(klon, knon, nsrf, zxli,
     &               uzon, vmer, tair1, qair1, zgeo1,
     &               tairsol, qairsol, rugo1, psfce, patm,
cIM  &               yt2m, yq2m, yu10m)
     &               yt2m, yq2m, yt10m, yq10m, yu10m, yustar)
cIM 081204 END
c
c
      DO j=1, knon
       i = ni(j)
       t2m(i,nsrf)=yt2m(j)

c
       q2m(i,nsrf)=yq2m(j)
c
c u10m, v10m : composantes du vent a 10m sans spirale de Ekman
       u10m(i,nsrf)=(yu10m(j) * uzon(j))/sqrt(uzon(j)**2+vmer(j)**2)
       v10m(i,nsrf)=(yu10m(j) * vmer(j))/sqrt(uzon(j)**2+vmer(j)**2)
c
      ENDDO
c
cIM cf AM : pbl, HBTM
      DO i = 1, knon
         y_cd_h(i) = ycoefh(i,1)
         y_cd_m(i) = ycoefm(i,1)
      ENDDO
c     print*,'appel hbtm2'
      CALL HBTM(knon, ypaprs, ypplay,
     .          yt2m,yt10m,yq2m,yq10m,yustar,
     .          y_flux_t,y_flux_q,yu,yv,yt,yq,
     .          ypblh,ycapCL,yoliqCL,ycteiCL,ypblT,
     .          ytherm,ytrmb1,ytrmb2,ytrmb3,ylcl)
c     print*,'fin hbtm2'
c
      DO j=1, knon
       i = ni(j)
       pblh(i,nsrf)   = ypblh(j)
       plcl(i,nsrf)   = ylcl(j)
       capCL(i,nsrf)  = ycapCL(j)
       oliqCL(i,nsrf) = yoliqCL(j)
       cteiCL(i,nsrf) = ycteiCL(j)
       pblT(i,nsrf)   = ypblT(j)
       therm(i,nsrf)  = ytherm(j)
       trmb1(i,nsrf)  = ytrmb1(j)
       trmb2(i,nsrf)  = ytrmb2(j)
       trmb3(i,nsrf)  = ytrmb3(j)
      ENDDO
c

      do j=1,knon
         do k=1,klev+1
         i=ni(j)
         q2(i,k,nsrf)=yq2(j,k)
         enddo
      enddo
cIM "slab" ocean 
       IF (nsrf.EQ.is_oce) THEN
        DO j = 1, knon
c on projette sur la grille globale
         i = ni(j)
         IF(pctsrf_new(i,is_oce).GT.epsfra) THEN
          flux_o(i) = y_flux_o(j)
         ELSE
          flux_o(i) = 0.
         ENDIF
        ENDDO
       ENDIF
c
       IF (nsrf.EQ.is_sic) THEN
        DO j = 1, knon
         i = ni(j)
cIM 230604 on pondere lorsque l'on fait le bilan au sol :  flux_g(i) = y_flux_g(j)*ypct(j)
         IF(pctsrf_new(i,is_sic).GT.epsfra) THEN
          flux_g(i) = y_flux_g(j)
         ELSE
          flux_g(i) = 0.
         ENDIF
        ENDDO
       ENDIF !nsrf.EQ.is_sic
c
      IF(OCEAN.EQ.'slab  ') THEN
       IF(nsrf.EQ.is_oce) then
        tslab(1:klon) = ytslab(1:klon)
        seaice(1:klon) = y_seaice(1:klon)
       ENDIF !nsrf
      ENDIF !OCEAN
99999 CONTINUE
C
C On utilise les nouvelles surfaces
C A rajouter: conservation de l'albedo
C
      rugos(:,is_oce) = rugmer
      pctsrf = pctsrf_new

      RETURN
      END
