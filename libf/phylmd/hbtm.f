      SUBROUTINE HBTM(knon, paprs, pplay,
     .                t2m,t10m,q2m,q10m,ustar,
     .                flux_t,flux_q,u,v,t,q,
     .                pblh,cape,EauLiq,ctei,pblT,
     .                therm,trmb1,trmb2,trmb3,plcl)
      use dimens_m
      use dimphy
      use YOMCST
      use yoethf
      use fcttre
        IMPLICIT none

c***************************************************************
c*                                                             *
c* HBTM2   D'apres Holstag&Boville et Troen&Mahrt              *
c*                 JAS 47              BLM                     *
c* Algorithme These Anne Mathieu                               *
c* Critere d'Entrainement Peter Duynkerke (JAS 50)             *
c* written by  : Anne MATHIEU & Alain LAHELLEC, 22/11/99       *
c* features : implem. exces Mathieu                            *
c***************************************************************
c* mods : decembre 99 passage th a niveau plus bas. voir fixer *
c* la prise du th a z/Lambda = -.2 (max Ray)                   *
c* Autre algo : entrainement ~ Theta+v =cste mais comment=>The?*
c* on peut fixer q a .7qsat(cf non adiab)=>T2 et The2          *
c* voir aussi //KE pblh = niveau The_e ou l = env.             *
c***************************************************************
c* fin therm a la HBTM passage a forme Mathieu 12/09/2001      *
c***************************************************************
c*
c
c
cAM Fev 2003
c Adaptation a LMDZ version couplee
c
c Pour le moment on fait passer en argument les grdeurs de surface : 
c flux, t,q2m, t,q10m, on va utiliser systematiquement les grdeurs a 2m ms 
c on garde la possibilite de changer si besoin est (jusqu'a present la 
c forme de HB avec le 1er niveau modele etait conservee)
c
c
c
c
c
      REAL RLvCp, REPS
c Arguments:
c
      INTEGER knon ! nombre de points a calculer
cAM
      REAL t2m(klon), t10m(klon) ! temperature a 2 et 10m
      REAL q2m(klon), q10m(klon) ! q a 2 et 10m
      REAL ustar(klon)
      REAL paprs(klon,klev+1) ! pression a inter-couche (Pa)
      REAL pplay(klon,klev)   ! pression au milieu de couche (Pa)
      REAL flux_t(klon,klev), flux_q(klon,klev)     ! Flux 
      REAL u(klon,klev) ! vitesse U (m/s)
      REAL v(klon,klev) ! vitesse V (m/s)
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! vapeur d'eau (kg/kg)
cAM      REAL cd_h(klon) ! coefficient de friction au sol pour chaleur
cAM      REAL cd_m(klon) ! coefficient de friction au sol pour vitesse
c
      INTEGER isommet
      PARAMETER (isommet=klev) ! limite max sommet pbl
      REAL vk
      PARAMETER (vk=0.35)     ! Von Karman => passer a .41 ! cf U.Olgstrom
      REAL ricr
      PARAMETER (ricr=0.4)
      REAL fak
      PARAMETER (fak=8.5)     ! b calcul du Prandtl et de dTetas
      REAL fakn
      PARAMETER (fakn=7.2)    ! a
      REAL onet
      PARAMETER (onet=1.0/3.0)
      REAL t_coup
      PARAMETER(t_coup=273.15)
      REAL zkmin
      PARAMETER (zkmin=0.01)
      REAL betam
      PARAMETER (betam=15.0)  ! pour Phim / h dans la S.L stable
      REAL betah
      PARAMETER (betah=15.0)
      REAL betas
      PARAMETER (betas=5.0)   ! Phit dans la S.L. stable (mais 2 formes / z/OBL<>1
      REAL sffrac
      PARAMETER (sffrac=0.1)  ! S.L. = z/h < .1
      REAL binm
      PARAMETER (binm=betam*sffrac)
      REAL binh
      PARAMETER (binh=betah*sffrac)
      REAL ccon
      PARAMETER (ccon=fak*sffrac*vk)
c
      REAL q_star,t_star
      REAL b1,b2,b212,b2sr     ! Lambert correlations T' q' avec T* q*
      PARAMETER (b1=70.,b2=20.)
c
      REAL z(klon,klev)
cAM      REAL pcfm(klon,klev), pcfh(klon,klev)
cAM
      REAL zref
      PARAMETER (zref=2.)    ! Niveau de ref a 2m peut eventuellement 
c                              etre choisi a 10m
cMA
c
      INTEGER i, k, j
      REAL zxt
cAM      REAL zxt, zxq, zxu, zxv, zxmod, taux, tauy
cAM      REAL zx_alf1, zx_alf2 ! parametres pour extrapolation
      REAL khfs(klon)       ! surface kinematic heat flux [mK/s]
      REAL kqfs(klon)       ! sfc kinematic constituent flux [m/s]
      REAL heatv(klon)      ! surface virtual heat flux
      REAL rhino(klon,klev) ! bulk Richardon no. mais en Theta_v
      LOGICAL unstbl(klon)  ! pts w/unstbl pbl (positive virtual ht flx)
      LOGICAL stblev(klon)  ! stable pbl with levels within pbl
      LOGICAL unslev(klon)  ! unstbl pbl with levels within pbl
      LOGICAL unssrf(klon)  ! unstb pbl w/lvls within srf pbl lyr
      LOGICAL unsout(klon)  ! unstb pbl w/lvls in outer pbl lyr
      LOGICAL check(klon)   ! True=>chk if Richardson no.>critcal
      LOGICAL omegafl(klon) ! flag de prolongerment cape pour pt Omega
      REAL pblh(klon)
      REAL pblT(klon)
      REAL plcl(klon)
cAM      REAL cgh(klon,2:klev) ! counter-gradient term for heat [K/m]
cAM      REAL cgq(klon,2:klev) ! counter-gradient term for constituents
cAM      REAL cgs(klon,2:klev) ! counter-gradient star (cg/flux)
      REAL obklen(klon)     ! Monin-Obukhov lengh
cAM      REAL ztvd, ztvu, 
      REAL zdu2
      REAL therm(klon)      ! thermal virtual temperature excess
      REAL trmb1(klon),trmb2(klon),trmb3(klon)
C  Algorithme thermique
      REAL s(klon,klev)     ! [P/Po]^Kappa milieux couches
      REAL Th_th(klon)      ! potential temperature of thermal
      REAL The_th(klon)     ! equivalent potential temperature of thermal
      REAL qT_th(klon)      ! total water  of thermal
      REAL Tbef(klon)       ! T thermique niveau precedent
      REAL qsatbef(klon)
      LOGICAL Zsat(klon)    ! le thermique est sature
      REAL Cape(klon)       ! Cape du thermique
      REAL Kape(klon)       ! Cape locale
      REAL EauLiq(klon)     ! Eau liqu integr du thermique
      REAL ctei(klon)       ! Critere d'instab d'entrainmt des nuages de CL
      REAL the1,the2,aa,bb,zthvd,zthvu,xintpos,qqsat
cIM 091204 BEG
      REAL a1,a2,a3
cIM 091204 END
      REAL xhis,rnum,denom,th1,th2,thv1,thv2,ql2
      REAL dqsat_dt,qsat2,qT1,q2,t1,t2,xnull,delt_the
      REAL delt_qt,delt_2,quadsat,spblh,reduc
c
      REAL phiminv(klon)    ! inverse phi function for momentum
      REAL phihinv(klon)    ! inverse phi function for heat
      REAL wm(klon)         ! turbulent velocity scale for momentum
      REAL fak1(klon)       ! k*ustar*pblh
      REAL fak2(klon)       ! k*wm*pblh
      REAL fak3(klon)       ! fakn*wstr/wm
      REAL pblk(klon)       ! level eddy diffusivity for momentum
      REAL pr(klon)         ! Prandtl number for eddy diffusivities
      REAL zl(klon)         ! zmzp / Obukhov length
      REAL zh(klon)         ! zmzp / pblh
      REAL zzh(klon)        ! (1-(zmzp/pblh))**2
      REAL wstr(klon)       ! w*, convective velocity scale
      REAL zm(klon)         ! current level height
      REAL zp(klon)         ! current level height + one level up
      REAL zcor, zdelta, zcvm5
cAM      REAL zxqs
      REAL fac, pblmin, zmzp, term
c



! initialisations (Anne)
      th_th(:) = 0.
      q_star = 0
      t_star = 0


      b212=sqrt(b1*b2)
      b2sr=sqrt(b2)
c
C ============================================================
C     Fonctions thermo implicites
C ============================================================
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Tetens : pression partielle de vap d'eau e_sat(T)
c =================================================
C++ e_sat(T) = r2*exp( r3*(T-Tf)/(T-r4) ) id a r2*FOEWE
C++ avec :
C++ Tf = 273.16 K  (Temp de fusion de la glace)
C++ r2 = 611.14 Pa
C++ r3 = 17.269 (liquide) 21.875 (solide) adim
C++ r4 = 35.86             7.66           Kelvin
C++  q_sat = eps*e_sat/(p-(1-eps)*e_sat)
C++ derivée :
C++ =========
C++                   r3*(Tf-r4)*q_sat(T,p)
C++ d_qsat_dT = --------------------------------
C++             (T-r4)^2*( 1-(1-eps)*e_sat(T)/p )
c++ pour zcvm5=Lv, c'est FOEDE
c++ Rq :(1.-REPS)*esarg/Parg id a RETV*Qsat
C     ------------------------------------------------------------------
c
c Initialisation
      RLvCp = RLVTT/RCPD
      REPS  = RD/RV

c
c      DO i = 1, klon
c         pcfh(i,1) = cd_h(i)
c         pcfm(i,1) = cd_m(i)
c      ENDDO
c      DO k = 2, klev
c      DO i = 1, klon
c         pcfh(i,k) = zkmin
c         pcfm(i,k) = zkmin
c         cgs(i,k) = 0.0
c         cgh(i,k) = 0.0
c         cgq(i,k) = 0.0
c      ENDDO
c      ENDDO
c
c Calculer les hauteurs de chaque couche
c (geopotentielle Int_dp/ro = Int_[Rd.T.dp/p] z = geop/g)
c  pourquoi ne pas utiliser Phi/RG ?
      DO i = 1, knon
         z(i,1) = RD * t(i,1) / (0.5*(paprs(i,1)+pplay(i,1)))
     .               * (paprs(i,1)-pplay(i,1)) / RG
         s(i,1) = (pplay(i,1)/paprs(i,1))**RKappa
      ENDDO
c                                 s(k) = [pplay(k)/ps]^kappa
c    + + + + + + + + + pplay  <-> s(k)   t  dp=pplay(k-1)-pplay(k)
c
c    -----------------  paprs <-> sig(k)
c
c    + + + + + + + + + pplay  <-> s(k-1)
c
c
c    + + + + + + + + + pplay  <-> s(1)   t  dp=paprs-pplay   z(1)
c
c    -----------------  paprs <-> sig(1)
c
      DO k = 2, klev
      DO i = 1, knon
         z(i,k) = z(i,k-1)
     .              + RD * 0.5*(t(i,k-1)+t(i,k)) / paprs(i,k)
     .                   * (pplay(i,k-1)-pplay(i,k)) / RG
         s(i,k) = (pplay(i,k)/paprs(i,1))**RKappa
      ENDDO
      ENDDO
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++  Determination des grandeurs de surface  +++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO i = 1, knon
cAM         IF (thermcep) THEN
cAM           zdelta=MAX(0.,SIGN(1.,RTT-tsol(i)))
c           zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
c           zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q(i,1))
cAM           zxqs= r2es * FOEEW(tsol(i),zdelta)/paprs(i,1)
cAM           zxqs=MIN(0.5,zxqs)
cAM           zcor=1./(1.-retv*zxqs)
cAM           zxqs=zxqs*zcor
cAM         ELSE
cAM           IF (tsol(i).LT.t_coup) THEN
cAM              zxqs = qsats(tsol(i)) / paprs(i,1)
cAM           ELSE
cAM              zxqs = qsatl(tsol(i)) / paprs(i,1)
cAM           ENDIF
cAM         ENDIF
c niveau de reference bulk; mais ici, c,a pourrait etre le niveau de ref du thermique
cAM        zx_alf1 = 1.0
cAM        zx_alf2 = 1.0 - zx_alf1
cAM        zxt = (t(i,1)+z(i,1)*RG/RCPD/(1.+RVTMP2*q(i,1)))
cAM     .        *(1.+RETV*q(i,1))*zx_alf1
cAM     .      + (t(i,2)+z(i,2)*RG/RCPD/(1.+RVTMP2*q(i,2)))
cAM     .        *(1.+RETV*q(i,2))*zx_alf2
cAM        zxu = u(i,1)*zx_alf1+u(i,2)*zx_alf2
cAM        zxv = v(i,1)*zx_alf1+v(i,2)*zx_alf2
cAM        zxq = q(i,1)*zx_alf1+q(i,2)*zx_alf2
cAM      
cAMAM           zxu = u10m(i)
cAMAM           zxv = v10m(i)
cAMAM           zxmod = 1.0+SQRT(zxu**2+zxv**2)
cAM Niveau de ref choisi a 2m
        zxt = t2m(i)

c ***************************************************
c attention, il doit s'agir de <w'theta'>
c   ;Calcul de tcls virtuel et de w'theta'virtuel
c   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
c   tcls=tcls*(1+.608*qcls)
c
c   ;Pour avoir w'theta',
c   ; il faut diviser par ro.Cp
c   Cp=Cpd*(1+0.84*qcls)
c   fcs=fcs/(ro_surf*Cp)
c   ;On transforme w'theta' en w'thetav'
c   Lv=(2.501-0.00237*(tcls-273.15))*1.E6
c   xle=xle/(ro_surf*Lv)
c   fcsv=fcs+.608*xle*tcls
c ***************************************************
cAM        khfs(i) = (tsol(i)*(1.+RETV*q(i,1))-zxt) *zxmod*cd_h(i)
cAM        kqfs(i) = (zxqs-zxq) *zxmod*cd_h(i) * beta(i)
cAM
cdif khfs est deja w't'_v / heatv(i) = khfs(i) + RETV*zxt*kqfs(i)
cAM calcule de Ro = paprs(i,1)/Rd zxt
cAM convention >0 vers le bas ds lmdz 
        khfs(i) = - flux_t(i,1)*zxt*Rd / (RCPD*paprs(i,1))
        kqfs(i) = - flux_q(i,1)*zxt*Rd / (paprs(i,1))
cAM   verifier que khfs et kqfs sont bien de la forme w'l'
        heatv(i) = khfs(i) + 0.608*zxt*kqfs(i)
c a comparer aussi aux sorties de clqh : flux_T/RoCp et flux_q/RoLv
cAM        heatv(i) = khfs(i)
cAM ustar est en entree
cAM        taux = zxu *zxmod*cd_m(i)
cAM        tauy = zxv *zxmod*cd_m(i)
cAM        ustar(i) = SQRT(taux**2+tauy**2)
cAM        ustar(i) = MAX(SQRT(ustar(i)),0.01)
c Theta et qT du thermique sans exces (interpolin vers surf)
c chgt de niveau du thermique (jeudi 30/12/1999)
c (interpolation lineaire avant integration phi_h)
cAM        qT_th(i) = zxqs*beta(i) + 4./z(i,1)*(q(i,1)-zxqs*beta(i))
cAM        qT_th(i) = max(qT_th(i),q(i,1))
        qT_th(i) = q2m(i)
cn The_th restera la Theta du thermique sans exces jusqu'a 2eme calcul
cn reste a regler convention P) pour Theta
c        The_th(i) = tsol(i) + 4./z(i,1)*(t(i,1)-tsol(i))
c     -                      + RLvCp*qT_th(i)
cAM        Th_th(i) = tsol(i) + 4./z(i,1)*(t(i,1)-tsol(i))
        Th_th(i) = t2m(i)
      ENDDO
c
      DO i = 1, knon
         rhino(i,1) = 0.0   ! Global Richardson
         check(i) = .TRUE.
         pblh(i) = z(i,1)   ! on initialise pblh a l'altitude du 1er niveau
         plcl(i) = 6000.
c Lambda = -u*^3 / (alpha.g.kvon.<w'Theta'v>
         obklen(i) = -t(i,1)*ustar(i)**3/(RG*vk*heatv(i))
         trmb1(i)   = 0.
         trmb2(i)   = 0.
         trmb3(i) = 0.
      ENDDO

C
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PBL height calculation:
C Search for level of pbl. Scan upward until the Richardson number between
C the first level and the current level exceeds the "critical" value.
C (bonne idee Nu de separer le Ric et l'exces de temp du thermique)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      fac = 100.0
      DO k = 2, isommet
      DO i = 1, knon
      IF (check(i)) THEN
! pourquoi / niveau 1 (au lieu du sol) et le terme en u*^2 ?
ctest     zdu2 = (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
         zdu2 = u(i,k)**2+v(i,k)**2
         zdu2 = max(zdu2,1.0e-20)
c Theta_v environnement
         zthvd=t(i,k)/s(i,k)*(1.+RETV*q(i,k))
c
c therm Theta_v sans exces (avec hypothese fausse de H&B, sinon,
c passer par Theta_e et virpot)
c         zthvu=t(i,1)/s(i,1)*(1.+RETV*q(i,1))
cAM         zthvu = Th_th(i)*(1.+RETV*q(i,1))
         zthvu = Th_th(i)*(1.+RETV*qT_th(i))
c  Le Ri par Theta_v
cAM         rhino(i,k) = (z(i,k)-z(i,1))*RG*(zthvd-zthvu)
cAM     .               /(zdu2*0.5*(zthvd+zthvu))
cAM On a nveau de ref a 2m ???
         rhino(i,k) = (z(i,k)-zref)*RG*(zthvd-zthvu)
     .               /(zdu2*0.5*(zthvd+zthvu))
c
         IF (rhino(i,k).GE.ricr) THEN
           pblh(i) = z(i,k-1) + (z(i,k-1)-z(i,k)) *
     .              (ricr-rhino(i,k-1))/(rhino(i,k-1)-rhino(i,k))
c test04
           pblh(i) = pblh(i) + 100.
           pblT(i) = t(i,k-1) + (t(i,k)-t(i,k-1)) *
     .              (pblh(i)-z(i,k-1))/(z(i,k)-z(i,k-1))
           check(i) = .FALSE.
         ENDIF
      ENDIF
      ENDDO
      ENDDO

C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
      DO i = 1, knon
        if (check(i)) pblh(i) = z(i,isommet)
      ENDDO
C
C Improve estimate of pbl height for the unstable points.
C Find unstable points (sensible heat flux is upward):
C
      DO i = 1, knon
      IF (heatv(i) .GT. 0.) THEN
        unstbl(i) = .TRUE.
        check(i) = .TRUE.
      ELSE
        unstbl(i) = .FALSE.
        check(i) = .FALSE.
      ENDIF
      ENDDO
C
C For the unstable case, compute velocity scale and the
C convective temperature excess:
C
      DO i = 1, knon
      IF (check(i)) THEN
        phiminv(i) = (1.-binm*pblh(i)/obklen(i))**onet
c ***************************************************
c Wm ? et W* ? c'est la formule pour z/h < .1
c   ;Calcul de w* ;;
c   ;;;;;;;;;;;;;;;;
c   w_star=((g/tcls)*fcsv*z(ind))^(1/3.) [ou prendre la premiere approx de h)
c   ;; CALCUL DE wm ;;
c   ;;;;;;;;;;;;;;;;;;
c   ; Ici on considerera que l'on est dans la couche de surf jusqu'a 100m
c   ; On prend svt couche de surface=0.1*h mais on ne connait pas h
c   ;;;;;;;;;;;Dans la couche de surface
c   if (z(ind) le 20) then begin
c   Phim=(1.-15.*(z(ind)/L))^(-1/3.)
c   wm=u_star/Phim
c   ;;;;;;;;;;;En dehors de la couche de surface
c   endif else if (z(ind) gt 20) then begin
c   wm=(u_star^3+c1*w_star^3)^(1/3.)
c   endif
c ***************************************************
        wm(i)= ustar(i)*phiminv(i)
c======================================================================
cvaleurs de Dominique Lambert de la campagne SEMAPHORE :
c <T'^2> = 100.T*^2; <q'^2> = 20.q*^2 a 10m
c <Tv'^2> = (1+1.2q).100.T* + 1.2Tv.sqrt(20*100).T*.q* + (.608*Tv)^2*20.q*^2;
c et dTetavS = sqrt(<Tv'^2>) ainsi calculee.
c avec : T*=<w'T'>_s/w* et q*=<w'q'>/w*
c !!! on peut donc utiliser w* pour les fluctuations <-> Lambert
c(leur corellation pourrait dependre de beta par ex)
c  if fcsv(i,j) gt 0 then begin
c    dTetavs=b1*(1.+2.*.608*q_10(i,j))*(fcs(i,j)/wm(i,j))^2+$
c    (.608*Thetav_10(i,j))^2*b2*(xle(i,j)/wm(i,j))^2+$
c    2.*.608*thetav_10(i,j)*sqrt(b1*b2)*(xle(i,j)/wm(i,j))*(fcs(i,j)/wm(i,j))
c    dqs=b2*(xle(i,j)/wm(i,j))^2
c    theta_s(i,j)=thetav_10(i,j)+sqrt(dTetavs)
c    q_s(i,j)=q_10(i,j)+sqrt(dqs)
c  endif else begin
c    Theta_s(i,j)=thetav_10(i,j)
c    q_s(i,j)=q_10(i,j)
c  endelse
c======================================================================
c
cHBTM        therm(i) = heatv(i)*fak/wm(i)
c forme Mathieu :
        q_star = kqfs(i)/wm(i)
        t_star = khfs(i)/wm(i)
cIM 091204 BEG
        IF(1.EQ.0) THEN
        IF(t_star.LT.0..OR.q_star.LT.0.) THEN
          print*,'i t_star q_star khfs kqfs wm',i,t_star,q_star,
     $    khfs(i),kqfs(i),wm(i)
        ENDIF
        ENDIF
cIM 091204 END
cAM Nveau cde ref 2m =>
cAM        therm(i) = sqrt( b1*(1.+2.*RETV*q(i,1))*t_star**2
cAM     +             + (RETV*T(i,1))**2*b2*q_star**2
cAM     +             + 2.*RETV*T(i,1)*b212*q_star*t_star
cAM     +                 )
cIM 091204 BEG
        a1=b1*(1.+2.*RETV*qT_th(i))*t_star**2
        a2=(RETV*Th_th(i))**2*b2*q_star**2
        a3=2.*RETV*Th_th(i)*b212*q_star*t_star
        aa=a1+a2+a3
        IF(1.EQ.0) THEN
        IF (aa.LT.0.) THEN 
         print*,'i a1 a2 a3 aa',i,a1,a2,a3,aa
         print*,'i qT_th Th_th t_star q_star RETV b1 b2 b212',
     $   i,qT_th(i),Th_th(i),t_star,q_star,RETV,b1,b2,b212
        ENDIF
        ENDIF
cIM 091204 END
        therm(i) = sqrt( b1*(1.+2.*RETV*qT_th(i))*t_star**2
     +             + (RETV*Th_th(i))**2*b2*q_star**2
cIM 101204  +             + 2.*RETV*Th_th(i)*b212*q_star*t_star
     +             + max(0.,2.*RETV*Th_th(i)*b212*q_star*t_star)
     +                 )
c
c Theta et qT du thermique (forme H&B) avec exces
c (attention, on ajoute therm(i) qui est virtuelle ...)
c pourquoi pas sqrt(b1)*t_star ?
c        dqs = b2sr*kqfs(i)/wm(i)
        qT_th(i) = qT_th(i)  + b2sr*q_star
cnew on differre le calcul de Theta_e
c        The_th(i) = The_th(i) + therm(i) + RLvCp*qT_th(i)
c ou:    The_th(i) = The_th(i) + sqrt(b1)*khfs(i)/wm(i) + RLvCp*qT_th(i)
        rhino(i,1) = 0.0
      ENDIF
      ENDDO
C
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ++ Improve pblh estimate for unstable conditions using the +++++++
C ++          convective temperature excess :                +++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO k = 2, isommet
      DO i = 1, knon
      IF (check(i)) THEN
ctest     zdu2 = (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
         zdu2 = u(i,k)**2+v(i,k)**2
         zdu2 = max(zdu2,1.0e-20)
c Theta_v environnement
         zthvd=t(i,k)/s(i,k)*(1.+RETV*q(i,k))
c
c et therm Theta_v (avec hypothese de constance de H&B,
c         zthvu=(t(i,1)+therm(i))/s(i,1)*(1.+RETV*q(i,1))
         zthvu = Th_th(i)*(1.+RETV*qT_th(i)) + therm(i)

c
c  Le Ri par Theta_v
cAM Niveau de ref 2m
cAM         rhino(i,k) = (z(i,k)-z(i,1))*RG*(zthvd-zthvu)
cAM     .               /(zdu2*0.5*(zthvd+zthvu))
         rhino(i,k) = (z(i,k)-zref)*RG*(zthvd-zthvu)
     .               /(zdu2*0.5*(zthvd+zthvu))
c
c
         IF (rhino(i,k).GE.ricr) THEN
           pblh(i) = z(i,k-1) + (z(i,k-1)-z(i,k)) *
     .              (ricr-rhino(i,k-1))/(rhino(i,k-1)-rhino(i,k))
c test04
           pblh(i) = pblh(i) + 100.
           pblT(i) = t(i,k-1) + (t(i,k)-t(i,k-1)) *
     .              (pblh(i)-z(i,k-1))/(z(i,k)-z(i,k-1))
           check(i) = .FALSE.
cIM 170305 BEG
      IF(1.EQ.0) THEN
c debug print -120;34       -34-        58 et    0;26 wamp
      if (i.eq.950.or.i.eq.192.or.i.eq.624.or.i.eq.118) then
            print*,' i,Th_th,Therm,qT :',i,Th_th(i),therm(i),qT_th(i)
            q_star = kqfs(i)/wm(i)
            t_star = khfs(i)/wm(i)
            print*,'q* t*, b1,b2,b212 ',q_star,t_star
     -            , b1*(1.+2.*RETV*qT_th(i))*t_star**2
     -            , (RETV*Th_th(i))**2*b2*q_star**2
     -            , 2.*RETV*Th_th(i)*b212*q_star*t_star
            print*,'zdu2 ,100.*ustar(i)**2',zdu2 ,fac*ustar(i)**2
      endif
      ENDIF !(1.EQ.0) THEN
cIM 170305 END
c             q_star = kqfs(i)/wm(i)
c             t_star = khfs(i)/wm(i)
c             trmb1(i) = b1*(1.+2.*RETV*q(i,1))*t_star**2
c             trmb2(i) = (RETV*T(i,1))**2*b2*q_star**2
c Omega now   trmb3(i) = 2.*RETV*T(i,1)*b212*q_star*t_star
         ENDIF
      ENDIF
      ENDDO
      ENDDO
C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
      DO i = 1, knon
        if (check(i)) pblh(i) = z(i,isommet)
      ENDDO
C
C PBL height must be greater than some minimum mechanical mixing depth
C Several investigators have proposed minimum mechanical mixing depth
C relationships as a function of the local friction velocity, u*.  We
C make use of a linear relationship of the form h = c u* where c=700.
C The scaling arguments that give rise to this relationship most often
C represent the coefficient c as some constant over the local coriolis
C parameter.  Here we make use of the experimental results of Koracin
C and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
C where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
C latitude value for f so that c = 0.07/f = 700.
C
      DO i = 1, knon
        pblmin  = 700.0*ustar(i)
        pblh(i) = MAX(pblh(i),pblmin)
c par exemple :
        pblT(i) = t(i,2) + (t(i,3)-t(i,2)) *
     .              (pblh(i)-z(i,2))/(z(i,3)-z(i,2))
      ENDDO

C ********************************************************************
C  pblh is now available; do preparation for diffusivity calculation :
C ********************************************************************
      DO i = 1, knon
        check(i) = .TRUE.
        Zsat(i)   = .FALSE.
c omegafl utilise pour prolongement CAPE
        omegafl(i) = .FALSE.
        Cape(i)   = 0.
        Kape(i)   = 0.
        EauLiq(i) = 0.
        CTEI(i)   = 0.
        pblk(i) = 0.0
        fak1(i) = ustar(i)*pblh(i)*vk
C
C Do additional preparation for unstable cases only, set temperature
C and moisture perturbations depending on stability.
C *** Rq: les formule sont prises dans leur forme CS ***
        IF (unstbl(i)) THEN
cAM Niveau de ref du thermique
cAM          zxt=(t(i,1)-z(i,1)*0.5*RG/RCPD/(1.+RVTMP2*q(i,1)))
cAM     .         *(1.+RETV*q(i,1))
          zxt=(Th_th(i)-zref*0.5*RG/RCPD/(1.+RVTMP2*qT_th(i)))
     .         *(1.+RETV*qT_th(i))
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1. - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*RG*pblh(i)/zxt)**onet
          fak3(i)    = fakn*wstr(i)/wm(i)
        ENDIF
c Computes Theta_e for thermal (all cases : to be modified)
c   attention ajout therm(i) = virtuelle
        The_th(i) = Th_th(i) + therm(i) + RLvCp*qT_th(i)
c ou:    The_th(i) = Th_th(i) + sqrt(b1)*khfs(i)/wm(i) + RLvCp*qT_th(i)
      ENDDO

C Main level loop to compute the diffusivities and
C counter-gradient terms:
C
      DO 1000 k = 2, isommet
C
C Find levels within boundary layer:
C
        DO i = 1, knon
          unslev(i) = .FALSE.
          stblev(i) = .FALSE.
          zm(i) = z(i,k-1)
          zp(i) = z(i,k)
          IF (zkmin.EQ.0.0 .AND. zp(i).GT.pblh(i)) zp(i) = pblh(i)
          IF (zm(i) .LT. pblh(i)) THEN
            zmzp = 0.5*(zm(i) + zp(i))
C debug
c          if (i.EQ.1864) then
c             print*,'i,pblh(1864),obklen(1864)',i,pblh(i),obklen(i)
c          endif

            zh(i) = zmzp/pblh(i)
            zl(i) = zmzp/obklen(i)
            zzh(i) = 0.
            IF (zh(i).LE.1.0) zzh(i) = (1. - zh(i))**2
C
C stblev for points zm < plbh and stable and neutral
C unslev for points zm < plbh and unstable
C
            IF (unstbl(i)) THEN
              unslev(i) = .TRUE.
            ELSE
              stblev(i) = .TRUE.
            ENDIF
          ENDIF
        ENDDO
c        print*,'fin calcul niveaux'
C
C Stable and neutral points; set diffusivities; counter-gradient
C terms zero for stable case:
C
        DO i = 1, knon
          IF (stblev(i)) THEN
            IF (zl(i).LE.1.) THEN
              pblk(i) = fak1(i)*zh(i)*zzh(i)/(1. + betas*zl(i))
            ELSE
              pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas + zl(i))
            ENDIF
c            pcfm(i,k) = pblk(i)
c            pcfh(i,k) = pcfm(i,k)
          ENDIF
        ENDDO
C
C unssrf, unstable within surface layer of pbl
C unsout, unstable within outer   layer of pbl
C
        DO i = 1, knon
          unssrf(i) = .FALSE.
          unsout(i) = .FALSE.
          IF (unslev(i)) THEN
            IF (zh(i).lt.sffrac) THEN
              unssrf(i) = .TRUE.
            ELSE
              unsout(i) = .TRUE.
            ENDIF
          ENDIF
        ENDDO
C
C Unstable for surface layer; counter-gradient terms zero
C
        DO i = 1, knon
          IF (unssrf(i)) THEN
            term = (1. - betam*zl(i))**onet
            pblk(i) = fak1(i)*zh(i)*zzh(i)*term
            pr(i) = term/sqrt(1. - betah*zl(i))
          ENDIF
        ENDDO
c        print*,'fin counter-gradient terms zero'
C
C Unstable for outer layer; counter-gradient terms non-zero:
C
        DO i = 1, knon
          IF (unsout(i)) THEN
            pblk(i) = fak2(i)*zh(i)*zzh(i)
c            cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
c            cgh(i,k) = khfs(i)*cgs(i,k)
            pr(i) = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
c            cgq(i,k) = kqfs(i)*cgs(i,k)
          ENDIF
        ENDDO
c        print*,'fin counter-gradient terms non zero'
C
C For all unstable layers, compute diffusivities and ctrgrad ter m
C
c        DO i = 1, knon
c        IF (unslev(i)) THEN
c            pcfm(i,k) = pblk(i)
c            pcfh(i,k) = pblk(i)/pr(i)
c etc cf original
c        ENDIF
c        ENDDO
C
C For all layers, compute integral info and CTEI
C
        DO i = 1, knon
        if (check(i).or.omegafl(i)) then
          if (.not.Zsat(i)) then
c            Th2 = The_th(i) - RLvCp*qT_th(i)
            Th2 = Th_th(i)
            T2 = Th2*s(i,k)
c thermodyn functions
            zdelta=MAX(0.,SIGN(1.,RTT-T2))
            qqsat= r2es * FOEEW(T2,zdelta)/pplay(i,k)
            qqsat=MIN(0.5,qqsat)
            zcor=1./(1.-retv*qqsat)
            qqsat=qqsat*zcor
c
            if (qqsat.lt.qT_th(i)) then
c on calcule lcl
              if (k.eq.2) then
                plcl(i) = z(i,k)
              else
                plcl(i) =  z(i,k-1) + (z(i,k-1)-z(i,k)) *
     .                 (qT_th(i)-qsatbef(i))/(qsatbef(i)-qqsat)
              endif
              Zsat(i) = .true.
              Tbef(i) = T2
            endif
c
          endif
          qsatbef(i) = qqsat
camn ???? cette ligne a deja ete faite normalement ?
        endif
c            print*,'hbtm2 i,k=',i,k
        ENDDO
 1000 continue           ! end of level loop
cIM 170305 BEG
        IF(1.EQ.0) THEN
            print*,'hbtm2  ok'
        ENDIF !(1.EQ.0) THEN
cIM 170305 END
      RETURN
      END
