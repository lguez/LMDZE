!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/conflx.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      SUBROUTINE conflx (dtime,pres_h,pres_f,
     e                   t, q, con_t, con_q, pqhfl, w,
     s                   d_t, d_q, rain, snow,
     s                   pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
     s                   kcbot, kctop, kdtop, pmflxr, pmflxs)
c
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      use fcttre
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19941014
c Objet: Schema flux de masse pour la convection 
c        (schema de Tiedtke avec qqs modifications mineures)
c Dec.97: Prise en compte des modifications introduites par
c         Olivier Boucher et Alexandre Armengaud pour melange
c         et lessivage des traceurs passifs.
c======================================================================
c Entree:
      REAL, intent(in):: dtime            ! pas d'integration (s)
      REAL, intent(in):: pres_h(klon,klev+1) ! pression half-level (Pa)
      REAL, intent(in):: pres_f(klon,klev)! pression full-level (Pa)
      REAL, intent(in):: t(klon,klev)     ! temperature (K)
      REAL q(klon,klev)     ! humidite specifique (g/g)
      REAL w(klon,klev)     ! vitesse verticale (Pa/s)
      REAL con_t(klon,klev) ! convergence de temperature (K/s)
      REAL con_q(klon,klev) ! convergence de l'eau vapeur (g/g/s)
      REAL pqhfl(klon)      ! evaporation (negative vers haut) mm/s
c Sortie:
      REAL d_t(klon,klev)   ! incrementation de temperature
      REAL d_q(klon,klev)   ! incrementation d'humidite
      REAL pmfu(klon,klev)  ! flux masse (kg/m2/s) panache ascendant
      REAL pmfd(klon,klev)  ! flux masse (kg/m2/s) panache descendant
      REAL pen_u(klon,klev)
      REAL pen_d(klon,klev)
      REAL pde_u(klon,klev)
      REAL pde_d(klon,klev)
      REAL rain(klon)       ! pluie (mm/s)
      REAL snow(klon)       ! neige (mm/s)
      REAL pmflxr(klon,klev+1)
      REAL pmflxs(klon,klev+1)
      INTEGER kcbot(klon)  ! niveau du bas de la convection
      INTEGER kctop(klon)  ! niveau du haut de la convection
      INTEGER kdtop(klon)  ! niveau du haut des downdrafts
c Local:
      REAL pt(klon,klev)
      REAL pq(klon,klev)
      REAL pqs(klon,klev)
      REAL pvervel(klon,klev)
      LOGICAL land(klon)
c
      REAL d_t_bis(klon,klev)
      REAL d_q_bis(klon,klev)
      REAL paprs(klon,klev+1)
      REAL paprsf(klon,klev)
      REAL zgeom(klon,klev)
      REAL zcvgq(klon,klev)
      REAL zcvgt(klon,klev)
cAA
      REAL zmfu(klon,klev) 
      REAL zmfd(klon,klev)
      REAL zen_u(klon,klev)
      REAL zen_d(klon,klev)
      REAL zde_u(klon,klev)
      REAL zde_d(klon,klev)
      REAL zmflxr(klon,klev+1)
      REAL zmflxs(klon,klev+1)
cAA

c
      INTEGER i, k
      REAL zdelta, zqsat
c
c
c initialiser les variables de sortie (pour securite)
      DO i = 1, klon
         rain(i) = 0.0
         snow(i) = 0.0
         kcbot(i) = 0
         kctop(i) = 0
         kdtop(i) = 0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
         pmfu(i,k) = 0.0
         pmfd(i,k) = 0.0
         pen_u(i,k) = 0.0
         pde_u(i,k) = 0.0
         pen_d(i,k) = 0.0
         pde_d(i,k) = 0.0
         zmfu(i,k) = 0.0
         zmfd(i,k) = 0.0
         zen_u(i,k) = 0.0
         zde_u(i,k) = 0.0
         zen_d(i,k) = 0.0
         zde_d(i,k) = 0.0
      ENDDO
      ENDDO
      DO k = 1, klev+1
      DO i = 1, klon
         zmflxr(i,k) = 0.0
         zmflxs(i,k) = 0.0
      ENDDO
      ENDDO
c
c calculer la nature du sol (pour l'instant, ocean partout)
      DO i = 1, klon
         land(i) = .FALSE.
      ENDDO
c
c preparer les variables d'entree (attention: l'ordre des niveaux 
c verticaux augmente du haut vers le bas)
      DO k = 1, klev
      DO i = 1, klon
         pt(i,k) = t(i,klev-k+1)
         pq(i,k) = q(i,klev-k+1)
         paprsf(i,k) = pres_f(i,klev-k+1)
         paprs(i,k) = pres_h(i,klev+1-k+1)
         pvervel(i,k) = w(i,klev+1-k)
         zcvgt(i,k) = con_t(i,klev-k+1)
         zcvgq(i,k) = con_q(i,klev-k+1)
c
         zdelta=MAX(0.,SIGN(1.,RTT-pt(i,k)))
         zqsat=R2ES*FOEEW ( pt(i,k), zdelta ) / paprsf(i,k)
         zqsat=MIN(0.5,zqsat)
         zqsat=zqsat/(1.-RETV  *zqsat)
         pqs(i,k) = zqsat
      ENDDO
      ENDDO
      DO i = 1, klon
         paprs(i,klev+1) = pres_h(i,1)
         zgeom(i,klev) = RD * pt(i,klev)
     .                   / (0.5*(paprs(i,klev+1)+paprsf(i,klev)))
     .                   * (paprs(i,klev+1)-paprsf(i,klev))
      ENDDO
      DO k = klev-1, 1, -1
      DO i = 1, klon
         zgeom(i,k) = zgeom(i,k+1)
     .              + RD * 0.5*(pt(i,k+1)+pt(i,k)) / paprs(i,k+1)
     .                   * (paprsf(i,k+1)-paprsf(i,k))
      ENDDO
      ENDDO
c
c appeler la routine principale
c
      CALL flxmain(dtime, pt, pq, pqs, pqhfl,
     .             paprsf, paprs, zgeom, land, zcvgt, zcvgq, pvervel,
     .             rain, snow, kcbot, kctop, kdtop,
     .             zmfu, zmfd, zen_u, zde_u, zen_d, zde_d,
     .             d_t_bis, d_q_bis, zmflxr, zmflxs)
C
cAA--------------------------------------------------------
cAA rem : De la meme facon que l'on effectue le reindicage 
cAA       pour la temperature t et le champ q 
cAA       on reindice les flux necessaires a la convection 
cAA       des traceurs
cAA--------------------------------------------------------
      DO k = 1, klev
      DO i = 1, klon
         d_q(i,klev+1-k) = dtime*d_q_bis(i,k)
         d_t(i,klev+1-k) = dtime*d_t_bis(i,k)
      ENDDO
      ENDDO
c
      DO i = 1, klon
         pmfu(i,1)= 0.
         pmfd(i,1)= 0.
         pen_d(i,1)= 0.
         pde_d(i,1)= 0.
      ENDDO
     
      DO k = 2, klev
      DO i = 1, klon
         pmfu(i,klev+2-k)= zmfu(i,k)
         pmfd(i,klev+2-k)= zmfd(i,k)
      ENDDO
      ENDDO
c
      DO k = 1, klev
      DO i = 1, klon
         pen_u(i,klev+1-k)=  zen_u(i,k)
         pde_u(i,klev+1-k)=  zde_u(i,k)
      ENDDO
      ENDDO
c
      DO k = 1, klev-1
      DO i = 1, klon
         pen_d(i,klev+1-k)= -zen_d(i,k+1)
         pde_d(i,klev+1-k)= -zde_d(i,k+1)
      ENDDO
      ENDDO

      DO k = 1, klev+1
      DO i = 1, klon
         pmflxr(i,klev+2-k)= zmflxr(i,k)
         pmflxs(i,klev+2-k)= zmflxs(i,k)
      ENDDO
      ENDDO

      RETURN
      END
c--------------------------------------------------------------------
      SUBROUTINE flxmain(pdtime, pten, pqen, pqsen, pqhfl, pap, paph,
     .                   pgeo, ldland, ptte, pqte, pvervel,
     .                   prsfc, pssfc, kcbot, kctop, kdtop,
c     *                   ldcum, ktype,
     .                   pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
     .                   dt_con, dq_con, pmflxr, pmflxs)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
C     ------------------------------------------------------------------
C     ----------------------------------------------------------------
      REAL pten(klon,klev), pqen(klon,klev), pqsen(klon,klev)
      REAL ptte(klon,klev)
      REAL pqte(klon,klev)
      REAL pvervel(klon,klev)
      REAL pgeo(klon,klev), pap(klon,klev), paph(klon,klev+1)
      REAL pqhfl(klon)
c
      REAL ptu(klon,klev), pqu(klon,klev), plu(klon,klev)
      REAL plude(klon,klev)
      REAL pmfu(klon,klev)
      REAL prsfc(klon), pssfc(klon)
      INTEGER  kcbot(klon), kctop(klon), ktype(klon)
      LOGICAL  ldland(klon), ldcum(klon)
c
      REAL ztenh(klon,klev), zqenh(klon,klev), zqsenh(klon,klev)
      REAL zgeoh(klon,klev)
      REAL zmfub(klon), zmfub1(klon)
      REAL zmfus(klon,klev), zmfuq(klon,klev), zmful(klon,klev)
      REAL zdmfup(klon,klev), zdpmel(klon,klev)
      REAL zentr(klon), zhcbase(klon)
      REAL zdqpbl(klon), zdqcv(klon), zdhpbl(klon)
      REAL zrfl(klon)
      REAL pmflxr(klon,klev+1)
      REAL pmflxs(klon,klev+1)
      INTEGER  ilab(klon,klev), ictop0(klon)
      LOGICAL  llo1
      REAL dt_con(klon,klev), dq_con(klon,klev)
      REAL zmfmax, zdh
      REAL, intent(in):: pdtime
      real zqumqe, zdqmin, zalvdcp, zhsat, zzz
      REAL zhhat, zpbmpt, zgam, zeps, zfac
      INTEGER i, k, ikb, itopm2, kcum
c
      REAL pen_u(klon,klev), pde_u(klon,klev)
      REAL pen_d(klon,klev), pde_d(klon,klev)
c
      REAL ptd(klon,klev), pqd(klon,klev), pmfd(klon,klev)
      REAL zmfds(klon,klev), zmfdq(klon,klev), zdmfdp(klon,klev)
      INTEGER kdtop(klon)
      LOGICAL lddraf(klon)
C---------------------------------------------------------------------
      LOGICAL firstcal
      SAVE firstcal
      DATA firstcal / .TRUE. /
C---------------------------------------------------------------------
      IF (firstcal) THEN
         CALL flxsetup
         firstcal = .FALSE.
      ENDIF
C---------------------------------------------------------------------
      DO i = 1, klon
         ldcum(i) = .FALSE.
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         dt_con(i,k) = 0.0
         dq_con(i,k) = 0.0
      ENDDO
      ENDDO
c----------------------------------------------------------------------
c initialiser les variables et faire l'interpolation verticale
c----------------------------------------------------------------------
      CALL flxini(pten, pqen, pqsen, pgeo,
     .     paph, zgeoh, ztenh, zqenh, zqsenh,
     .     ptu, pqu, ptd, pqd, pmfd, zmfds, zmfdq, zdmfdp,
     .     pmfu, zmfus, zmfuq, zdmfup,
     .     zdpmel, plu, plude, ilab, pen_u, pde_u, pen_d, pde_d)
c---------------------------------------------------------------------
c determiner les valeurs au niveau de base de la tour convective
c---------------------------------------------------------------------
      CALL flxbase(ztenh, zqenh, zgeoh, paph,
     *            ptu, pqu, plu, ldcum, kcbot, ilab)
c---------------------------------------------------------------------
c calculer la convergence totale de l'humidite et celle en provenance
c de la couche limite, plus precisement, la convergence integree entre
c le sol et la base de la convection. Cette derniere convergence est
c comparee avec l'evaporation obtenue dans la couche limite pour
c determiner le type de la convection
c---------------------------------------------------------------------
      k=1
      DO i = 1, klon
         zdqcv(i) = pqte(i,k)*(paph(i,k+1)-paph(i,k))
         zdhpbl(i) = 0.0
         zdqpbl(i) = 0.0
      ENDDO
c
      DO k=2,klev
      DO i = 1, klon
          zdqcv(i)=zdqcv(i)+pqte(i,k)*(paph(i,k+1)-paph(i,k))
          IF (k.GE.kcbot(i)) THEN
             zdqpbl(i)=zdqpbl(i)+pqte(i,k)*(paph(i,k+1)-paph(i,k))
             zdhpbl(i)=zdhpbl(i)+(RCPD*ptte(i,k)+RLVTT*pqte(i,k))
     .                          *(paph(i,k+1)-paph(i,k))
          ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
         ktype(i) = 2
         if (zdqcv(i).GT.MAX(0.,-1.5*pqhfl(i)*RG)) ktype(i) = 1
ccc         if (zdqcv(i).GT.MAX(0.,-1.1*pqhfl(i)*RG)) ktype(i) = 1
      ENDDO
c
c---------------------------------------------------------------------
c determiner le flux de masse entrant a travers la base.
c on ignore, pour l'instant, l'effet du panache descendant
c---------------------------------------------------------------------
      DO i = 1, klon
         ikb=kcbot(i)
         zqumqe=pqu(i,ikb)+plu(i,ikb)-zqenh(i,ikb)
         zdqmin=MAX(0.01*zqenh(i,ikb),1.E-10)
         IF (zdqpbl(i).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(i)) THEN
            zmfub(i) = zdqpbl(i)/(RG*MAX(zqumqe,zdqmin))
         ELSE
            zmfub(i) = 0.01
            ldcum(i)=.FALSE.
         ENDIF
         IF (ktype(i).EQ.2) THEN
            zdh = RCPD*(ptu(i,ikb)-ztenh(i,ikb)) + RLVTT*zqumqe
            zdh = RG * MAX(zdh,1.0E5*zdqmin)
            IF (zdhpbl(i).GT.0..AND.ldcum(i))zmfub(i)=zdhpbl(i)/zdh
         ENDIF
         zmfmax = (paph(i,ikb)-paph(i,ikb-1)) / (RG*pdtime)
         zmfub(i) = MIN(zmfub(i),zmfmax)
         zentr(i) = ENTRSCV
         IF (ktype(i).EQ.1) zentr(i) = ENTRPEN
      ENDDO
C-----------------------------------------------------------------------
C DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
C-----------------------------------------------------------------------
c (A) calculer d'abord la hauteur "theorique" de la tour convective sans
c     considerer l'entrainement ni le detrainement du panache, sachant
c     ces derniers peuvent abaisser la hauteur theorique.
c
      DO i = 1, klon
         ikb=kcbot(i)
         zhcbase(i)=RCPD*ptu(i,ikb)+zgeoh(i,ikb)+RLVTT*pqu(i,ikb)
         ictop0(i)=kcbot(i)-1
      ENDDO
c
      zalvdcp=RLVTT/RCPD
      DO k=klev-1,3,-1
      DO i = 1, klon
         zhsat=RCPD*ztenh(i,k)+zgeoh(i,k)+RLVTT*zqsenh(i,k)
         zgam=R5LES*zalvdcp*zqsenh(i,k)/
     .        ((1.-RETV  *zqsenh(i,k))*(ztenh(i,k)-R4LES)**2)
         zzz=RCPD*ztenh(i,k)*0.608
         zhhat=zhsat-(zzz+zgam*zzz)/(1.+zgam*zzz/RLVTT)*
     .               MAX(zqsenh(i,k)-zqenh(i,k),0.)
         IF(k.LT.ictop0(i).AND.zhcbase(i).GT.zhhat) ictop0(i)=k
      ENDDO
      ENDDO
c
c (B) calculer le panache ascendant
c
      CALL flxasc(pdtime,ztenh, zqenh, pten, pqen, pqsen,
     .     pgeo, zgeoh, pap, paph, pqte, pvervel,
     .     ldland, ldcum, ktype, ilab,
     .     ptu, pqu, plu, pmfu, zmfub, zentr,
     .     zmfus, zmfuq, zmful, plude, zdmfup,
     .     kcbot, kctop, ictop0, kcum, pen_u, pde_u)
      IF (kcum.EQ.0) GO TO 1000
C
C verifier l'epaisseur de la convection et changer eventuellement
c le taux d'entrainement/detrainement
C
      DO i = 1, klon
         zpbmpt=paph(i,kcbot(i))-paph(i,kctop(i))
         IF(ldcum(i).AND.ktype(i).EQ.1.AND.zpbmpt.LT.2.E4)ktype(i)=2
         IF(ldcum(i)) ictop0(i)=kctop(i)
         IF(ktype(i).EQ.2) zentr(i)=ENTRSCV
      ENDDO
c
      IF (lmfdd) THEN  ! si l'on considere le panache descendant
c
c calculer la precipitation issue du panache ascendant pour 
c determiner l'existence du panache descendant dans la convection
      DO i = 1, klon
         zrfl(i)=zdmfup(i,1)
      ENDDO
      DO k=2,klev
      DO i = 1, klon
         zrfl(i)=zrfl(i)+zdmfup(i,k)
      ENDDO
      ENDDO
c
c determiner le LFS (level of free sinking: niveau de plonge libre)
      CALL flxdlfs(ztenh, zqenh, zgeoh, paph, ptu, pqu,
     *     ldcum,    kcbot,    kctop,    zmfub,    zrfl,
     *     ptd,      pqd,
     *     pmfd,     zmfds,    zmfdq,    zdmfdp,
     *     kdtop,    lddraf)
c
c calculer le panache descendant
      CALL flxddraf(ztenh,    zqenh,
     *     zgeoh,    paph,     zrfl,
     *     ptd,      pqd,
     *     pmfd,     zmfds,    zmfdq,    zdmfdp,
     *     lddraf, pen_d, pde_d)
c
c calculer de nouveau le flux de masse entrant a travers la base
c de la convection, sachant qu'il a ete modifie par le panache
c descendant
      DO i = 1, klon
      IF (lddraf(i)) THEN
         ikb = kcbot(i)
         llo1 = PMFD(i,ikb).LT.0.
         zeps = 0.
         IF ( llo1 ) zeps = CMFDEPS
         zqumqe = pqu(i,ikb)+plu(i,ikb)-
     .            zeps*pqd(i,ikb)-(1.-zeps)*zqenh(i,ikb)
         zdqmin = MAX(0.01*zqenh(i,ikb),1.E-10)
         zmfmax = (paph(i,ikb)-paph(i,ikb-1)) / (RG*pdtime)
         IF (zdqpbl(i).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(i)
     .       .AND.zmfub(i).LT.zmfmax) THEN
            zmfub1(i) = zdqpbl(i) / (RG*MAX(zqumqe,zdqmin))
         ELSE
            zmfub1(i) = zmfub(i)
         ENDIF
         IF (ktype(i).EQ.2) THEN
            zdh = RCPD*(ptu(i,ikb)-zeps*ptd(i,ikb)-
     .            (1.-zeps)*ztenh(i,ikb))+RLVTT*zqumqe
            zdh = RG * MAX(zdh,1.0E5*zdqmin)
            IF (zdhpbl(i).GT.0..AND.ldcum(i))zmfub1(i)=zdhpbl(i)/zdh
         ENDIF
         IF ( .NOT.((ktype(i).EQ.1.OR.ktype(i).EQ.2).AND.
     .              ABS(zmfub1(i)-zmfub(i)).LT.0.2*zmfub(i)) )
     .      zmfub1(i) = zmfub(i)
      ENDIF
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
      IF (lddraf(i)) THEN
         zfac = zmfub1(i)/MAX(zmfub(i),1.E-10)
         pmfd(i,k) = pmfd(i,k)*zfac
         zmfds(i,k) = zmfds(i,k)*zfac
         zmfdq(i,k) = zmfdq(i,k)*zfac
         zdmfdp(i,k) = zdmfdp(i,k)*zfac
         pen_d(i,k) = pen_d(i,k)*zfac
         pde_d(i,k) = pde_d(i,k)*zfac
      ENDIF
      ENDDO
      ENDDO
      DO i = 1, klon
         IF (lddraf(i)) zmfub(i)=zmfub1(i)
      ENDDO
c
      ENDIF   ! fin de test sur lmfdd
c
c-----------------------------------------------------------------------
c calculer de nouveau le panache ascendant
c-----------------------------------------------------------------------
      CALL flxasc(pdtime,ztenh, zqenh, pten, pqen, pqsen,
     .     pgeo, zgeoh, pap, paph, pqte, pvervel,
     .     ldland, ldcum, ktype, ilab,
     .     ptu, pqu, plu, pmfu, zmfub, zentr,
     .     zmfus, zmfuq, zmful, plude, zdmfup,
     .     kcbot, kctop, ictop0, kcum, pen_u, pde_u)
c
c-----------------------------------------------------------------------
c determiner les flux convectifs en forme finale, ainsi que
c la quantite des precipitations
c-----------------------------------------------------------------------
      CALL flxflux(pdtime, pqen, pqsen, ztenh, zqenh, pap, paph, 
     .     ldland, zgeoh, kcbot, kctop, lddraf, kdtop, ktype, ldcum,
     .     pmfu, pmfd, zmfus, zmfds, zmfuq, zmfdq, zmful, plude,
     .     zdmfup, zdmfdp, pten, prsfc, pssfc, zdpmel, itopm2,
     .     pmflxr, pmflxs)
c
c----------------------------------------------------------------------
c calculer les tendances pour T et Q
c----------------------------------------------------------------------
      CALL flxdtdq(itopm2, paph, ldcum, pten,
     e     zmfus, zmfds, zmfuq, zmfdq, zmful, zdmfup, zdmfdp, zdpmel,
     s     dt_con,dq_con)
c
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE flxini(pten, pqen, pqsen, pgeo, paph, pgeoh, ptenh,
     .           pqenh, pqsenh, ptu, pqu, ptd, pqd, pmfd, pmfds, pmfdq,
     .           pdmfdp, pmfu, pmfus, pmfuq, pdmfup, pdpmel, plu, plude,
     .           klab,pen_u, pde_u, pen_d, pde_d)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      IMPLICIT none
C----------------------------------------------------------------------
C THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
C TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
C AND INITIALIZES VALUES FOR UPDRAFTS
C----------------------------------------------------------------------
C
      REAL pten(klon,klev)   ! temperature (environnement)
      REAL pqen(klon,klev)   ! humidite (environnement)
      REAL pqsen(klon,klev)  ! humidite saturante (environnement)
      REAL pgeo(klon,klev)   ! geopotentiel (g * metre)
      REAL pgeoh(klon,klev)  ! geopotentiel aux demi-niveaux
      REAL paph(klon,klev+1) ! pression aux demi-niveaux
      REAL ptenh(klon,klev)  ! temperature aux demi-niveaux
      REAL pqenh(klon,klev)  ! humidite aux demi-niveaux
      REAL pqsenh(klon,klev) ! humidite saturante aux demi-niveaux
C
      REAL ptu(klon,klev)    ! temperature du panache ascendant (p-a)
      REAL pqu(klon,klev)    ! humidite du p-a
      REAL plu(klon,klev)    ! eau liquide du p-a
      REAL pmfu(klon,klev)   ! flux de masse du p-a
      REAL pmfus(klon,klev)  ! flux de l'energie seche dans le p-a
      REAL pmfuq(klon,klev)  ! flux de l'humidite dans le p-a
      REAL pdmfup(klon,klev) ! quantite de l'eau precipitee dans p-a
      REAL plude(klon,klev)  ! quantite de l'eau liquide jetee du
c                              p-a a l'environnement
      REAL pdpmel(klon,klev) ! quantite de neige fondue
c
      REAL ptd(klon,klev)    ! temperature du panache descendant (p-d)
      REAL pqd(klon,klev)    ! humidite du p-d
      REAL pmfd(klon,klev)   ! flux de masse du p-d
      REAL pmfds(klon,klev)  ! flux de l'energie seche dans le p-d
      REAL pmfdq(klon,klev)  ! flux de l'humidite dans le p-d
      REAL pdmfdp(klon,klev) ! quantite de precipitation dans p-d
c
      REAL pen_u(klon,klev) ! quantite de masse entrainee pour p-a
      REAL pde_u(klon,klev) ! quantite de masse detrainee pour p-a
      REAL pen_d(klon,klev) ! quantite de masse entrainee pour p-d
      REAL pde_d(klon,klev) ! quantite de masse detrainee pour p-d
C
      INTEGER  klab(klon,klev)
      LOGICAL  llflag(klon)
      INTEGER k, i, icall
      REAL zzs
C----------------------------------------------------------------------
C SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
C ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
C----------------------------------------------------------------------
      DO 130 k = 2, klev
c
      DO i = 1, klon
         pgeoh(i,k)=pgeo(i,k)+(pgeo(i,k-1)-pgeo(i,k))*0.5
         ptenh(i,k)=(MAX(RCPD*pten(i,k-1)+pgeo(i,k-1),
     .             RCPD*pten(i,k)+pgeo(i,k))-pgeoh(i,k))/RCPD
         pqsenh(i,k)=pqsen(i,k-1)
         llflag(i)=.TRUE.
      ENDDO
c
      icall=0
      CALL flxadjtq(paph(1,k),ptenh(1,k),pqsenh(1,k),llflag,icall)
c
      DO i = 1, klon
         pqenh(i,k)=MIN(pqen(i,k-1),pqsen(i,k-1))
     .               +(pqsenh(i,k)-pqsen(i,k-1))
         pqenh(i,k)=MAX(pqenh(i,k),0.)
      ENDDO
c
  130 CONTINUE
C
      DO 140 i = 1, klon
         ptenh(i,klev)=(RCPD*pten(i,klev)+pgeo(i,klev)-
     1                   pgeoh(i,klev))/RCPD
         pqenh(i,klev)=pqen(i,klev)
         ptenh(i,1)=pten(i,1)
         pqenh(i,1)=pqen(i,1)
         pgeoh(i,1)=pgeo(i,1)
  140 CONTINUE
c
      DO 160 k = klev-1, 2, -1
      DO 150 i = 1, klon
         zzs = MAX(RCPD*ptenh(i,k)+pgeoh(i,k),
     .             RCPD*ptenh(i,k+1)+pgeoh(i,k+1))
         ptenh(i,k) = (zzs-pgeoh(i,k))/RCPD
  150 CONTINUE
  160 CONTINUE
C
C-----------------------------------------------------------------------
C INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
C-----------------------------------------------------------------------
      DO k = 1, klev
      DO i = 1, klon
         ptu(i,k) = ptenh(i,k)
         pqu(i,k) = pqenh(i,k)
         plu(i,k) = 0.
         pmfu(i,k) = 0.
         pmfus(i,k) = 0.
         pmfuq(i,k) = 0.
         pdmfup(i,k) = 0.
         pdpmel(i,k) = 0.
         plude(i,k) = 0.
c
         klab(i,k) = 0
c
         ptd(i,k) = ptenh(i,k)
         pqd(i,k) = pqenh(i,k)
         pmfd(i,k) = 0.0
         pmfds(i,k) = 0.0
         pmfdq(i,k) = 0.0
         pdmfdp(i,k) = 0.0
c
         pen_u(i,k) = 0.0
         pde_u(i,k) = 0.0
         pen_d(i,k) = 0.0
         pde_d(i,k) = 0.0
      ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE flxbase(ptenh, pqenh, pgeoh, paph,
     *     ptu, pqu, plu, ldcum, kcbot, klab)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      IMPLICIT none
C----------------------------------------------------------------------
C THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
C
C INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
C IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
C   klab=1 FOR SUBCLOUD LEVELS
C   klab=2 FOR CONDENSATION LEVEL
C
C LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
C (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
C----------------------------------------------------------------------
C       ----------------------------------------------------------------
      REAL ptenh(klon,klev), pqenh(klon,klev)
      REAL pgeoh(klon,klev), paph(klon,klev+1)
C
      REAL ptu(klon,klev), pqu(klon,klev), plu(klon,klev)
      INTEGER  klab(klon,klev), kcbot(klon)
C
      LOGICAL llflag(klon), ldcum(klon)
      INTEGER i, k, icall, is
      REAL zbuo, zqold(klon)
C----------------------------------------------------------------------
C INITIALIZE VALUES AT LIFTING LEVEL
C----------------------------------------------------------------------
      DO i = 1, klon
         klab(i,klev)=1
         kcbot(i)=klev-1
         ldcum(i)=.FALSE.
      ENDDO
C----------------------------------------------------------------------
C DO ASCENT IN SUBCLOUD LAYER,
C CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
C ADJUST T,Q AND L ACCORDINGLY
C CHECK FOR BUOYANCY AND SET FLAGS
C----------------------------------------------------------------------
      DO 290 k = klev-1, 2, -1
c
      is = 0
      DO i = 1, klon
         IF (klab(i,k+1).EQ.1) is = is + 1
         llflag(i) = .FALSE.
         IF (klab(i,k+1).EQ.1) llflag(i) = .TRUE.
      ENDDO
      IF (is.EQ.0) GOTO 290
c
      DO i = 1, klon
      IF(llflag(i)) THEN
         pqu(i,k) = pqu(i,k+1)
         ptu(i,k) = ptu(i,k+1)+(pgeoh(i,k+1)-pgeoh(i,k))/RCPD
         zbuo = ptu(i,k)*(1.+RETV*pqu(i,k))-
     .          ptenh(i,k)*(1.+RETV*pqenh(i,k))+0.5
         IF (zbuo.GT.0.) klab(i,k) = 1
         zqold(i) = pqu(i,k)
      ENDIF
      ENDDO
c
      icall=1
      CALL flxadjtq(paph(1,k), ptu(1,k), pqu(1,k), llflag, icall)
c
      DO i = 1, klon
      IF (llflag(i).AND.pqu(i,k).NE.zqold(i)) THEN
         klab(i,k) = 2
         plu(i,k) = plu(i,k) + zqold(i)-pqu(i,k)
         zbuo = ptu(i,k)*(1.+RETV*pqu(i,k))-
     .          ptenh(i,k)*(1.+RETV*pqenh(i,k))+0.5
         IF (zbuo.GT.0.) kcbot(i) = k
         IF (zbuo.GT.0.) ldcum(i) = .TRUE.
      ENDIF
      ENDDO
c
  290 CONTINUE
c
      RETURN
      END
      SUBROUTINE flxasc(pdtime, ptenh, pqenh, pten, pqen, pqsen,
     .     pgeo, pgeoh, pap, paph, pqte, pvervel,
     .     ldland, ldcum, ktype, klab, ptu, pqu, plu,
     .     pmfu, pmfub, pentr, pmfus, pmfuq,
     .     pmful, plude, pdmfup, kcbot, kctop, kctop0, kcum,
     .     pen_u, pde_u)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
C----------------------------------------------------------------------
C THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
C FOR CUMULUS PARAMETERIZATION
C----------------------------------------------------------------------
C
      REAL, intent(in):: pdtime
      REAL pten(klon,klev), ptenh(klon,klev)
      REAL pqen(klon,klev), pqenh(klon,klev), pqsen(klon,klev)
      REAL pgeo(klon,klev), pgeoh(klon,klev)
      REAL pap(klon,klev), paph(klon,klev+1)
      REAL pqte(klon,klev)
      REAL pvervel(klon,klev) ! vitesse verticale en Pa/s
C
      REAL pmfub(klon), pentr(klon)
      REAL ptu(klon,klev), pqu(klon,klev), plu(klon,klev)
      REAL plude(klon,klev)
      REAL pmfu(klon,klev), pmfus(klon,klev)
      REAL pmfuq(klon,klev), pmful(klon,klev)
      REAL pdmfup(klon,klev)
      INTEGER ktype(klon), klab(klon,klev), kcbot(klon), kctop(klon)
      INTEGER kctop0(klon)
      LOGICAL ldland(klon), ldcum(klon)
C
      REAL pen_u(klon,klev), pde_u(klon,klev)
      REAL zqold(klon)
      REAL zdland(klon)
      LOGICAL llflag(klon)
      INTEGER k, i, is, icall, kcum
      REAL ztglace, zdphi, zqeen, zseen, zscde, zqude
      REAL zmfusk, zmfuqk, zmfulk, zbuo, zdnoprc, zprcon, zlnew
c
      REAL zpbot(klon), zptop(klon), zrho(klon)
      REAL zdprho, zentr, zpmid, zmftest, zmfmax
      LOGICAL llo1, llo2
c
      REAL zwmax(klon), zzzmb
      INTEGER klwmin(klon) ! level of maximum vertical velocity
C----------------------------------------------------------------------
      ztglace = RTT - 13.
c
c Chercher le niveau ou la vitesse verticale est maximale:
      DO i = 1, klon
         klwmin(i) = klev
         zwmax(i) = 0.0
      ENDDO
      DO k = klev, 3, -1
      DO i = 1, klon
      IF (pvervel(i,k).LT.zwmax(i)) THEN
         zwmax(i) = pvervel(i,k)
         klwmin(i) = k
      ENDIF
      ENDDO
      ENDDO
C----------------------------------------------------------------------
C SET DEFAULT VALUES
C----------------------------------------------------------------------
      DO i = 1, klon
         IF (.NOT.ldcum(i)) ktype(i)=0
      ENDDO
c
      DO k=1,klev
      DO i = 1, klon
         plu(i,k)=0.
         pmfu(i,k)=0.
         pmfus(i,k)=0.
         pmfuq(i,k)=0.
         pmful(i,k)=0.
         plude(i,k)=0.
         pdmfup(i,k)=0.
         IF(.NOT.ldcum(i).OR.ktype(i).EQ.3) klab(i,k)=0
         IF(.NOT.ldcum(i).AND.paph(i,k).LT.4.E4) kctop0(i)=k
      ENDDO
      ENDDO
c
      DO i = 1, klon
      IF (ldland(i)) THEN
         zdland(i)=3.0E4
         zdphi=pgeoh(i,kctop0(i))-pgeoh(i,kcbot(i))
         IF (ptu(i,kctop0(i)).GE.ztglace) zdland(i)=zdphi
         zdland(i)=MAX(3.0E4,zdland(i))
         zdland(i)=MIN(5.0E4,zdland(i))
      ENDIF
      ENDDO
C
C Initialiser les valeurs au niveau d'ascendance
C
      DO i = 1, klon
         kctop(i) = klev-1
         IF (.NOT.ldcum(i)) THEN
            kcbot(i) = klev-1
            pmfub(i) = 0.
            pqu(i,klev) = 0.
         ENDIF
         pmfu(i,klev) = pmfub(i)
         pmfus(i,klev) = pmfub(i)*(RCPD*ptu(i,klev)+pgeoh(i,klev))
         pmfuq(i,klev) = pmfub(i)*pqu(i,klev)
      ENDDO
c
      DO i = 1, klon
         ldcum(i) = .FALSE.
      ENDDO
C----------------------------------------------------------------------
C  DO ASCENT: SUBCLOUD LAYER (klab=1) ,CLOUDS (klab=2)
C  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
C  BY ADJUSTING T,Q AND L ACCORDINGLY IN *flxadjtq*,
C  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
C----------------------------------------------------------------------
      DO 480 k = klev-1,3,-1
c
      IF (LMFMID .AND. k.LT.klev-1 .AND. k.GT.klev/2) THEN
         DO i = 1, klon
         IF (.NOT.ldcum(i) .AND. klab(i,k+1).EQ.0 .AND.
     .       pqen(i,k).GT.0.9*pqsen(i,k)) THEN
            ptu(i,k+1) = pten(i,k) +(pgeo(i,k)-pgeoh(i,k+1))/RCPD
            pqu(i,k+1) = pqen(i,k)
            plu(i,k+1) = 0.0
            zzzmb = MAX(CMFCMIN, -pvervel(i,k)/RG)
            zmfmax = (paph(i,k)-paph(i,k-1))/(RG*pdtime)
            pmfub(i) = MIN(zzzmb,zmfmax)
            pmfu(i,k+1) = pmfub(i)
            pmfus(i,k+1) = pmfub(i)*(RCPD*ptu(i,k+1)+pgeoh(i,k+1))
            pmfuq(i,k+1) = pmfub(i)*pqu(i,k+1)
            pmful(i,k+1) = 0.0
            pdmfup(i,k+1) = 0.0
            kcbot(i) = k
            klab(i,k+1) = 1
            ktype(i) = 3
            pentr(i) = ENTRMID
         ENDIF
         ENDDO
      ENDIF
c
      is = 0
      DO i = 1, klon
         is = is + klab(i,k+1)
         IF (klab(i,k+1) .EQ. 0) klab(i,k) = 0
         llflag(i) = .FALSE.
         IF (klab(i,k+1) .GT. 0) llflag(i) = .TRUE.
      ENDDO
      IF (is .EQ. 0) GOTO 480
c
c calculer le taux d'entrainement et de detrainement
c
      DO i = 1, klon
         pen_u(i,k) = 0.0
         pde_u(i,k) = 0.0
         zrho(i)=paph(i,k+1)/(RD*ptenh(i,k+1))
         zpbot(i)=paph(i,kcbot(i))
         zptop(i)=paph(i,kctop0(i))
      ENDDO
c
      DO 125 i = 1, klon
      IF(ldcum(i)) THEN
         zdprho=(paph(i,k+1)-paph(i,k))/(RG*zrho(i))
         zentr=pentr(i)*pmfu(i,k+1)*zdprho
         llo1=k.LT.kcbot(i)
         IF(llo1) pde_u(i,k)=zentr
         zpmid=0.5*(zpbot(i)+zptop(i))
         llo2=llo1.AND.ktype(i).EQ.2.AND.
     .        (zpbot(i)-paph(i,k).LT.0.2E5.OR.
     .         paph(i,k).GT.zpmid)
         IF(llo2) pen_u(i,k)=zentr
         llo2=llo1.AND.(ktype(i).EQ.1.OR.ktype(i).EQ.3).AND.
     .        (k.GE.MAX(klwmin(i),kctop0(i)+2).OR.pap(i,k).GT.zpmid)
         IF(llo2) pen_u(i,k)=zentr
         llo1=pen_u(i,k).GT.0..AND.(ktype(i).EQ.1.OR.ktype(i).EQ.2)
         IF(llo1) THEN
            zentr=zentr*(1.+3.*(1.-MIN(1.,(zpbot(i)-pap(i,k))/1.5E4)))
            pen_u(i,k)=pen_u(i,k)*(1.+3.*(1.-MIN(1.,
     .                 (zpbot(i)-pap(i,k))/1.5E4)))
            pde_u(i,k)=pde_u(i,k)*(1.+3.*(1.-MIN(1.,
     .                 (zpbot(i)-pap(i,k))/1.5E4)))
         ENDIF
         IF(llo2.AND.pqenh(i,k+1).GT.1.E-5)
     .   pen_u(i,k)=zentr+MAX(pqte(i,k),0.)/pqenh(i,k+1)*
     .              zrho(i)*zdprho
      ENDIF
  125 CONTINUE
c
C----------------------------------------------------------------------
c DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
C----------------------------------------------------------------------
c
      DO 420 i = 1, klon
      IF (llflag(i)) THEN
         IF (k.LT.kcbot(i)) THEN
            zmftest = pmfu(i,k+1)+pen_u(i,k)-pde_u(i,k)
            zmfmax = MIN(zmftest,(paph(i,k)-paph(i,k-1))/(RG*pdtime))
            pen_u(i,k)=MAX(pen_u(i,k)-MAX(0.0,zmftest-zmfmax),0.0)
         ENDIF
         pde_u(i,k)=MIN(pde_u(i,k),0.75*pmfu(i,k+1))
c calculer le flux de masse du niveau k a partir de celui du k+1
         pmfu(i,k)=pmfu(i,k+1)+pen_u(i,k)-pde_u(i,k)
c calculer les valeurs Su, Qu et l du niveau k dans le panache montant
         zqeen=pqenh(i,k+1)*pen_u(i,k)
         zseen=(RCPD*ptenh(i,k+1)+pgeoh(i,k+1))*pen_u(i,k)
         zscde=(RCPD*ptu(i,k+1)+pgeoh(i,k+1))*pde_u(i,k)
         zqude=pqu(i,k+1)*pde_u(i,k)
         plude(i,k)=plu(i,k+1)*pde_u(i,k)
         zmfusk=pmfus(i,k+1)+zseen-zscde
         zmfuqk=pmfuq(i,k+1)+zqeen-zqude
         zmfulk=pmful(i,k+1)    -plude(i,k)
         plu(i,k)=zmfulk*(1./MAX(CMFCMIN,pmfu(i,k)))
         pqu(i,k)=zmfuqk*(1./MAX(CMFCMIN,pmfu(i,k)))
         ptu(i,k)=(zmfusk*(1./MAX(CMFCMIN,pmfu(i,k)))-
     1               pgeoh(i,k))/RCPD
         ptu(i,k)=MAX(100.,ptu(i,k))
         ptu(i,k)=MIN(400.,ptu(i,k))
         zqold(i)=pqu(i,k)
      ELSE
         zqold(i)=0.0
      ENDIF
  420 CONTINUE
c
C----------------------------------------------------------------------
c DO CORRECTIONS FOR MOIST ASCENT BY ADJUSTING T,Q AND L
C----------------------------------------------------------------------
c
      icall = 1
      CALL flxadjtq(paph(1,k), ptu(1,k), pqu(1,k), llflag, icall)
C
      DO 440 i = 1, klon
      IF(llflag(i).AND.pqu(i,k).NE.zqold(i)) THEN
         klab(i,k) = 2
         plu(i,k) = plu(i,k)+zqold(i)-pqu(i,k)
         zbuo = ptu(i,k)*(1.+RETV*pqu(i,k))-
     .          ptenh(i,k)*(1.+RETV*pqenh(i,k))
         IF (klab(i,k+1).EQ.1) zbuo=zbuo+0.5
         IF (zbuo.GT.0..AND.pmfu(i,k).GE.0.1*pmfub(i)) THEN
            kctop(i) = k
            ldcum(i) = .TRUE.
            zdnoprc = 1.5E4
            IF (ldland(i)) zdnoprc = zdland(i)
            zprcon = CPRCON
            IF ((zpbot(i)-paph(i,k)).LT.zdnoprc) zprcon = 0.0
            zlnew=plu(i,k)/(1.+zprcon*(pgeoh(i,k)-pgeoh(i,k+1)))
            pdmfup(i,k)=MAX(0.,(plu(i,k)-zlnew)*pmfu(i,k))
            plu(i,k)=zlnew
         ELSE
            klab(i,k)=0
            pmfu(i,k)=0.
         ENDIF
      ENDIF
  440 CONTINUE
      DO 455 i = 1, klon
      IF (llflag(i)) THEN
         pmful(i,k)=plu(i,k)*pmfu(i,k)
         pmfus(i,k)=(RCPD*ptu(i,k)+pgeoh(i,k))*pmfu(i,k)
         pmfuq(i,k)=pqu(i,k)*pmfu(i,k)
      ENDIF
  455 CONTINUE
C
  480 CONTINUE
C----------------------------------------------------------------------
C DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
C    (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
C           AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
C           FROM PREVIOUS CALCULATIONS ABOVE)
C----------------------------------------------------------------------
      DO i = 1, klon
         IF (kctop(i).EQ.klev-1) ldcum(i) = .FALSE.
         kcbot(i) = MAX(kcbot(i),kctop(i))
      ENDDO
c
      ldcum(1)=ldcum(1)
c
      is = 0
      DO i = 1, klon
         if (ldcum(i)) is = is + 1
      ENDDO
      kcum = is
      IF (is.EQ.0) GOTO 800
c
      DO 530 i = 1, klon
      IF (ldcum(i)) THEN
         k=kctop(i)-1
         pde_u(i,k)=(1.-CMFCTOP)*pmfu(i,k+1)
         plude(i,k)=pde_u(i,k)*plu(i,k+1)
         pmfu(i,k)=pmfu(i,k+1)-pde_u(i,k)
         zlnew=plu(i,k)
         pdmfup(i,k)=MAX(0.,(plu(i,k)-zlnew)*pmfu(i,k))
         plu(i,k)=zlnew
         pmfus(i,k)=(RCPD*ptu(i,k)+pgeoh(i,k))*pmfu(i,k)
         pmfuq(i,k)=pqu(i,k)*pmfu(i,k)
         pmful(i,k)=plu(i,k)*pmfu(i,k)
         plude(i,k-1)=pmful(i,k)
      ENDIF
  530 CONTINUE
C
  800 CONTINUE
      RETURN
      END
      SUBROUTINE flxflux(pdtime, pqen, pqsen, ptenh, pqenh, pap
     .  ,  paph, ldland, pgeoh, kcbot, kctop, lddraf, kdtop
     .  ,  ktype, ldcum, pmfu, pmfd, pmfus, pmfds
     .  ,  pmfuq, pmfdq, pmful, plude, pdmfup, pdmfdp
     .  ,  pten, prfl, psfl, pdpmel, ktopm2
     .  ,  pmflxr, pmflxs)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      use fcttre
            use yoecumf
      IMPLICIT none
C----------------------------------------------------------------------
C THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
C FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
C----------------------------------------------------------------------
C
      REAL cevapcu(klev)
C     -----------------------------------------------------------------
      REAL pqen(klon,klev), pqenh(klon,klev), pqsen(klon,klev)
      REAL pten(klon,klev), ptenh(klon,klev)
      REAL paph(klon,klev+1), pgeoh(klon,klev)
c
      REAL pap(klon,klev)
      REAL ztmsmlt, zdelta, zqsat
C
      REAL pmfu(klon,klev), pmfus(klon,klev)
      REAL pmfd(klon,klev), pmfds(klon,klev)
      REAL pmfuq(klon,klev), pmful(klon,klev)
      REAL pmfdq(klon,klev)
      REAL plude(klon,klev)
      REAL pdmfup(klon,klev), pdpmel(klon,klev)
cjq The variable maxpdmfdp(klon) has been introduced by Olivier Boucher
cjq 14/11/00 to fix the problem with the negative precipitation.      
      REAL pdmfdp(klon,klev), maxpdmfdp(klon,klev) 
      REAL prfl(klon), psfl(klon)
      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1)
      INTEGER  kcbot(klon), kctop(klon), ktype(klon)
      LOGICAL  ldland(klon), ldcum(klon)
      INTEGER k, kp, i
      REAL zcons1, zcons2, zcucov, ztmelp2
      REAL, intent(in):: pdtime
      real zdp, zzp, zfac, zsnmlt, zrfl, zrnew
      REAL zrmin, zrfln, zdrfl
      REAL zpds, zpdr, zdenom
      INTEGER ktopm2, itop, ikb
c
      LOGICAL lddraf(klon)
      INTEGER kdtop(klon)
c
c
      DO 101 k=1,klev
      CEVAPCU(k)=1.93E-6*261.*SQRT(1.E3/(38.3*0.293)
     1 *SQRT(0.5*(paph(1,k)+paph(1,k+1))/paph(1,klev+1)) ) * 0.5/RG
 101  CONTINUE
c
c SPECIFY CONSTANTS
c
      zcons1 = RCPD/(RLMLT*RG*pdtime)
      zcons2 = 1./(RG*pdtime)
      zcucov = 0.05
      ztmelp2 = RTT + 2.
c
c DETERMINE FINAL CONVECTIVE FLUXES
c
      itop=klev
      DO 110 i = 1, klon
         itop=MIN(itop,kctop(i))
         IF (.NOT.ldcum(i) .OR. kdtop(i).LT.kctop(i)) lddraf(i)=.FALSE.
         IF(.NOT.ldcum(i)) ktype(i)=0
  110 CONTINUE
c
      ktopm2=itop-2
      DO 120 k=ktopm2,klev
      DO 115 i = 1, klon
      IF(ldcum(i).AND.k.GE.kctop(i)-1) THEN
         pmfus(i,k)=pmfus(i,k)-pmfu(i,k)*
     .                (RCPD*ptenh(i,k)+pgeoh(i,k))
         pmfuq(i,k)=pmfuq(i,k)-pmfu(i,k)*pqenh(i,k)
         zdp = 1.5E4
         IF ( ldland(i) ) zdp = 3.E4
c
c        l'eau liquide detrainee est precipitee quand certaines
c        conditions sont reunies (sinon, elle est consideree
c        evaporee dans l'environnement)
c
         IF(paph(i,kcbot(i))-paph(i,kctop(i)).GE.zdp.AND.
     .      pqen(i,k-1).GT.0.8*pqsen(i,k-1))
     .      pdmfup(i,k-1)=pdmfup(i,k-1)+plude(i,k-1)
c
         IF(lddraf(i).AND.k.GE.kdtop(i)) THEN
            pmfds(i,k)=pmfds(i,k)-pmfd(i,k)*
     .                   (RCPD*ptenh(i,k)+pgeoh(i,k))
            pmfdq(i,k)=pmfdq(i,k)-pmfd(i,k)*pqenh(i,k)
         ELSE
            pmfd(i,k)=0.
            pmfds(i,k)=0.
            pmfdq(i,k)=0.
            pdmfdp(i,k-1)=0.
         END IF
      ELSE
         pmfu(i,k)=0.
         pmfus(i,k)=0.
         pmfuq(i,k)=0.
         pmful(i,k)=0.
         pdmfup(i,k-1)=0.
         plude(i,k-1)=0.
         pmfd(i,k)=0.
         pmfds(i,k)=0.
         pmfdq(i,k)=0.
         pdmfdp(i,k-1)=0.
      ENDIF
  115 CONTINUE
  120 CONTINUE
c
      DO 130 k=ktopm2,klev
      DO 125 i = 1, klon
      IF(ldcum(i).AND.k.GT.kcbot(i)) THEN
         ikb=kcbot(i)
         zzp=((paph(i,klev+1)-paph(i,k))/
     .        (paph(i,klev+1)-paph(i,ikb)))
         IF (ktype(i).EQ.3) zzp = zzp**2
         pmfu(i,k)=pmfu(i,ikb)*zzp
         pmfus(i,k)=pmfus(i,ikb)*zzp
         pmfuq(i,k)=pmfuq(i,ikb)*zzp
         pmful(i,k)=pmful(i,ikb)*zzp
      ENDIF
  125 CONTINUE
  130 CONTINUE
c
c CALCULATE RAIN/SNOW FALL RATES
c CALCULATE MELTING OF SNOW
c CALCULATE EVAPORATION OF PRECIP
c
      DO k = 1, klev+1
      DO i = 1, klon
         pmflxr(i,k) = 0.0
         pmflxs(i,k) = 0.0
      ENDDO
      ENDDO
      DO k = ktopm2, klev
      DO i = 1, klon
      IF (ldcum(i)) THEN
         IF (pmflxs(i,k).GT.0.0 .AND. pten(i,k).GT.ztmelp2) THEN
            zfac=zcons1*(paph(i,k+1)-paph(i,k))
            zsnmlt=MIN(pmflxs(i,k),zfac*(pten(i,k)-ztmelp2))
            pdpmel(i,k)=zsnmlt
            ztmsmlt=pten(i,k)-zsnmlt/zfac
            zdelta=MAX(0.,SIGN(1.,RTT-ztmsmlt))
            zqsat=R2ES*FOEEW(ztmsmlt, zdelta) / pap(i,k)
            zqsat=MIN(0.5,zqsat)
            zqsat=zqsat/(1.-RETV  *zqsat)
            pqsen(i,k) = zqsat
         ENDIF
         IF (pten(i,k).GT.RTT) THEN
         pmflxr(i,k+1)=pmflxr(i,k)+pdmfup(i,k)+pdmfdp(i,k)+pdpmel(i,k)
         pmflxs(i,k+1)=pmflxs(i,k)-pdpmel(i,k)
         ELSE
           pmflxs(i,k+1)=pmflxs(i,k)+pdmfup(i,k)+pdmfdp(i,k)
           pmflxr(i,k+1)=pmflxr(i,k)
         ENDIF
c        si la precipitation est negative, on ajuste le plux du
c        panache descendant pour eliminer la negativite
         IF ((pmflxr(i,k+1)+pmflxs(i,k+1)).LT.0.0) THEN
            pdmfdp(i,k) = -pmflxr(i,k)-pmflxs(i,k)-pdmfup(i,k)
            pmflxr(i,k+1) = 0.0
            pmflxs(i,k+1) = 0.0
            pdpmel(i,k) = 0.0
         ENDIF
      ENDIF
      ENDDO
      ENDDO
c
cjq The new variable is initialized here.
cjq It contains the humidity which is fed to the downdraft
cjq by evaporation of precipitation in the column below the base
cjq of convection.
cjq 
cjq In the former version, this term has been subtracted from precip
cjq as well as the evaporation.
cjq      
      DO k = 1, klev
      DO i = 1, klon 
         maxpdmfdp(i,k)=0.0
      ENDDO
      ENDDO
      DO k = 1, klev
       DO kp = k, klev
        DO i = 1, klon
         maxpdmfdp(i,k)=maxpdmfdp(i,k)+pdmfdp(i,kp)
        ENDDO
       ENDDO
      ENDDO
cjq End of initialization
c      
      DO k = ktopm2, klev
      DO i = 1, klon
      IF (ldcum(i) .AND. k.GE.kcbot(i)) THEN
         zrfl = pmflxr(i,k) + pmflxs(i,k)
         IF (zrfl.GT.1.0E-20) THEN
            zrnew=(MAX(0.,SQRT(zrfl/zcucov)-
     .            CEVAPCU(k)*(paph(i,k+1)-paph(i,k))*
     .            MAX(0.,pqsen(i,k)-pqen(i,k))))**2*zcucov
            zrmin=zrfl-zcucov*MAX(0.,0.8*pqsen(i,k)-pqen(i,k))
     .            *zcons2*(paph(i,k+1)-paph(i,k))
            zrnew=MAX(zrnew,zrmin)
            zrfln=MAX(zrnew,0.)
            zdrfl=MIN(0.,zrfln-zrfl)
cjq At least the amount of precipiation needed to feed the downdraft
cjq with humidity below the base of convection has to be left and can't
cjq be evaporated (surely the evaporation can't be positive):            
            zdrfl=MAX(zdrfl,
     .            MIN(-pmflxr(i,k)-pmflxs(i,k)-maxpdmfdp(i,k),0.0))
cjq End of insertion
c            
            zdenom=1.0/MAX(1.0E-20,pmflxr(i,k)+pmflxs(i,k))
            IF (pten(i,k).GT.RTT) THEN
               zpdr = pdmfdp(i,k)
               zpds = 0.0
            ELSE
               zpdr = 0.0
               zpds = pdmfdp(i,k)
            ENDIF
            pmflxr(i,k+1) = pmflxr(i,k) + zpdr + pdpmel(i,k)
     .                    + zdrfl*pmflxr(i,k)*zdenom
            pmflxs(i,k+1) = pmflxs(i,k) + zpds - pdpmel(i,k)
     .                    + zdrfl*pmflxs(i,k)*zdenom
            pdmfup(i,k) = pdmfup(i,k) + zdrfl
         ELSE
            pmflxr(i,k+1) = 0.0
            pmflxs(i,k+1) = 0.0
            pdmfdp(i,k) = 0.0
            pdpmel(i,k) = 0.0
         ENDIF         
         if (pmflxr(i,k) + pmflxs(i,k).lt.-1.e-26) 
     .    write(*,*) 'precip. < 1e-16 ',pmflxr(i,k) + pmflxs(i,k)
      ENDIF
      ENDDO
      ENDDO
c
      DO 210 i = 1, klon
         prfl(i) = pmflxr(i,klev+1)
         psfl(i) = pmflxs(i,klev+1)
  210 CONTINUE
c
      RETURN
      END
      SUBROUTINE flxdtdq(ktopm2, paph, ldcum, pten
     .  ,  pmfus, pmfds, pmfuq, pmfdq, pmful, pdmfup, pdmfdp
     .  ,  pdpmel, dt_con, dq_con)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
c----------------------------------------------------------------------
c calculer les tendances T et Q
c----------------------------------------------------------------------
C     -----------------------------------------------------------------
      LOGICAL  llo1
C
      REAL pten(klon,klev), paph(klon,klev+1)
      REAL pmfus(klon,klev), pmfuq(klon,klev), pmful(klon,klev)
      REAL pmfds(klon,klev), pmfdq(klon,klev)
      REAL pdmfup(klon,klev)
      REAL pdmfdp(klon,klev)
      REAL pdpmel(klon,klev)
      LOGICAL ldcum(klon)
      REAL dt_con(klon,klev), dq_con(klon,klev)
c
      INTEGER ktopm2
c
      INTEGER i, k
      REAL zalv, zdtdt, zdqdt
c
      DO 210 k=ktopm2,klev-1
      DO 220 i = 1, klon
      IF (ldcum(i)) THEN
         llo1 = (pten(i,k)-RTT).GT.0.
         zalv = RLSTT
         IF (llo1) zalv = RLVTT
         zdtdt=RG/(paph(i,k+1)-paph(i,k))/RCPD
     .        *(pmfus(i,k+1)-pmfus(i,k)
     .         +pmfds(i,k+1)-pmfds(i,k)
     .          -RLMLT*pdpmel(i,k)
     .          -zalv*(pmful(i,k+1)-pmful(i,k)-pdmfup(i,k)-pdmfdp(i,k))
     .         )
         dt_con(i,k)=zdtdt
         zdqdt=RG/(paph(i,k+1)-paph(i,k))
     .        *(pmfuq(i,k+1)-pmfuq(i,k)
     .         +pmfdq(i,k+1)-pmfdq(i,k)
     .          +pmful(i,k+1)-pmful(i,k)-pdmfup(i,k)-pdmfdp(i,k))
         dq_con(i,k)=zdqdt
      ENDIF
  220 CONTINUE
  210 CONTINUE
C
      k = klev
      DO 230 i = 1, klon
      IF (ldcum(i)) THEN
         llo1 = (pten(i,k)-RTT).GT.0.
         zalv = RLSTT
         IF (llo1) zalv = RLVTT
         zdtdt=-RG/(paph(i,k+1)-paph(i,k))/RCPD
     .         *(pmfus(i,k)+pmfds(i,k)+RLMLT*pdpmel(i,k)
     .           -zalv*(pmful(i,k)+pdmfup(i,k)+pdmfdp(i,k)))
         dt_con(i,k)=zdtdt
         zdqdt=-RG/(paph(i,k+1)-paph(i,k))
     .            *(pmfuq(i,k)+pmfdq(i,k)+pmful(i,k)
     .             +pdmfup(i,k)+pdmfdp(i,k))
         dq_con(i,k)=zdqdt
      ENDIF
  230 CONTINUE
C
      RETURN
      END
      SUBROUTINE flxdlfs(ptenh, pqenh, pgeoh, paph, ptu, pqu,
     .     ldcum, kcbot, kctop, pmfub, prfl, ptd, pqd,
     .     pmfd, pmfds, pmfdq, pdmfdp, kdtop, lddraf)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
C
C----------------------------------------------------------------------
C THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
C CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
C
C TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
C FOR MASSFLUX CUMULUS PARAMETERIZATION
C
C INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
C AND UPDRAFT VALUES T,Q,U AND V AND ALSO
C CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
C IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
C
C CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
C MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
C----------------------------------------------------------------------
C
      REAL ptenh(klon,klev)
      REAL pqenh(klon,klev)
      REAL pgeoh(klon,klev), paph(klon,klev+1)
      REAL ptu(klon,klev), pqu(klon,klev)
      REAL pmfub(klon)
      REAL prfl(klon)
C
      REAL ptd(klon,klev), pqd(klon,klev)
      REAL pmfd(klon,klev), pmfds(klon,klev), pmfdq(klon,klev)
      REAL pdmfdp(klon,klev)
      INTEGER  kcbot(klon), kctop(klon), kdtop(klon)
      LOGICAL  ldcum(klon), lddraf(klon)
C
      REAL ztenwb(klon,klev), zqenwb(klon,klev), zcond(klon)
      REAL zttest, zqtest, zbuo, zmftop
      LOGICAL  llo2(klon)
      INTEGER i, k, is, icall
C----------------------------------------------------------------------
      DO i= 1, klon
         lddraf(i)=.FALSE.
         kdtop(i)=klev+1
      ENDDO
C
C----------------------------------------------------------------------
C DETERMINE LEVEL OF FREE SINKING BY
C DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
C
C FOR EVERY POINT AND PROCEED AS FOLLOWS:
C     (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
C     (2) DO MIXING WITH CUMULUS CLOUD AIR
C     (3) CHECK FOR NEGATIVE BUOYANCY
C
C THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
C OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
C TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
C EVAPORATION OF RAIN AND CLOUD WATER)
C----------------------------------------------------------------------
C
      DO 290 k = 3, klev-3
C
      is=0
      DO 212 i= 1, klon
         ztenwb(i,k)=ptenh(i,k)
         zqenwb(i,k)=pqenh(i,k)
         llo2(i) = ldcum(i).AND.prfl(i).GT.0.
     .             .AND..NOT.lddraf(i)
     .             .AND.(k.LT.kcbot(i).AND.k.GT.kctop(i))
         IF ( llo2(i) ) is = is + 1
  212 CONTINUE
      IF(is.EQ.0) GO TO 290
C
      icall=2
      CALL flxadjtq(paph(1,k), ztenwb(1,k), zqenwb(1,k), llo2, icall)
C
C----------------------------------------------------------------------
C DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
C AND CHECK FOR NEGATIVE BUOYANCY.
C THEN SET VALUES FOR DOWNDRAFT AT LFS.
C----------------------------------------------------------------------
      DO 222 i= 1, klon
      IF (llo2(i)) THEN
         zttest=0.5*(ptu(i,k)+ztenwb(i,k))
         zqtest=0.5*(pqu(i,k)+zqenwb(i,k))
         zbuo=zttest*(1.+RETV*zqtest)-
     .        ptenh(i,k)*(1.+RETV  *pqenh(i,k))
         zcond(i)=pqenh(i,k)-zqenwb(i,k)
         zmftop=-CMFDEPS*pmfub(i)
         IF (zbuo.LT.0..AND.prfl(i).GT.10.*zmftop*zcond(i)) THEN
            kdtop(i)=k
            lddraf(i)=.TRUE.
            ptd(i,k)=zttest
            pqd(i,k)=zqtest
            pmfd(i,k)=zmftop
            pmfds(i,k)=pmfd(i,k)*(RCPD*ptd(i,k)+pgeoh(i,k))
            pmfdq(i,k)=pmfd(i,k)*pqd(i,k)
            pdmfdp(i,k-1)=-0.5*pmfd(i,k)*zcond(i)
            prfl(i)=prfl(i)+pdmfdp(i,k-1)
         ENDIF
      ENDIF
  222 CONTINUE
c
  290 CONTINUE
C
      RETURN
      END
      SUBROUTINE flxddraf(ptenh, pqenh, pgeoh, paph, prfl,
     .           ptd, pqd, pmfd, pmfds, pmfdq, pdmfdp,
     .           lddraf, pen_d, pde_d)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
            use yoecumf
      IMPLICIT none
C
C----------------------------------------------------------------------
C          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
C
C          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
C          (I.E. T,Q,U AND V AND FLUXES)
C
C          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
C          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
C          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
C
C          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
C          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
C          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
C
C----------------------------------------------------------------------
C
      REAL ptenh(klon,klev), pqenh(klon,klev)
      REAL pgeoh(klon,klev), paph(klon,klev+1)
C
      REAL ptd(klon,klev), pqd(klon,klev)
      REAL pmfd(klon,klev), pmfds(klon,klev), pmfdq(klon,klev)
      REAL pdmfdp(klon,klev)
      REAL prfl(klon)
      LOGICAL lddraf(klon)
C
      REAL pen_d(klon,klev), pde_d(klon,klev), zcond(klon)
      LOGICAL llo2(klon), llo1
      INTEGER i, k, is, icall, itopde
      REAL zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zdmfdp
      REAL zbuo
C----------------------------------------------------------------------
C CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
C       (A) CALCULATING ENTRAINMENT RATES, ASSUMING
C           LINEAR DECREASE OF MASSFLUX IN PBL
C       (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
C           AND MOISTENING IS CALCULATED IN *flxadjtq*
C       (C) CHECKING FOR NEGATIVE BUOYANCY AND
C           SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
C
      DO 180 k = 3, klev
c
      is = 0
      DO i = 1, klon
         llo2(i)=lddraf(i).AND.pmfd(i,k-1).LT.0.
         IF (llo2(i)) is = is + 1
      ENDDO
      IF (is.EQ.0) GOTO 180
c
      DO i = 1, klon
      IF (llo2(i)) THEN
         zentr = ENTRDD*pmfd(i,k-1)*RD*ptenh(i,k-1)/
     .           (RG*paph(i,k-1))*(paph(i,k)-paph(i,k-1))
         pen_d(i,k) = zentr
         pde_d(i,k) = zentr
      ENDIF
      ENDDO
c
      itopde = klev-2
      IF (k.GT.itopde) THEN
         DO i = 1, klon
         IF (llo2(i)) THEN
            pen_d(i,k)=0.
            pde_d(i,k)=pmfd(i,itopde)*
     .      (paph(i,k)-paph(i,k-1))/(paph(i,klev+1)-paph(i,itopde))
         ENDIF
         ENDDO
      ENDIF
C
      DO i = 1, klon
      IF (llo2(i)) THEN
         pmfd(i,k) = pmfd(i,k-1)+pen_d(i,k)-pde_d(i,k)
         zseen = (RCPD*ptenh(i,k-1)+pgeoh(i,k-1))*pen_d(i,k)
         zqeen = pqenh(i,k-1)*pen_d(i,k)
         zsdde = (RCPD*ptd(i,k-1)+pgeoh(i,k-1))*pde_d(i,k)
         zqdde = pqd(i,k-1)*pde_d(i,k)
         zmfdsk = pmfds(i,k-1)+zseen-zsdde
         zmfdqk = pmfdq(i,k-1)+zqeen-zqdde
         pqd(i,k) = zmfdqk*(1./MIN(-CMFCMIN,pmfd(i,k)))
         ptd(i,k) = (zmfdsk*(1./MIN(-CMFCMIN,pmfd(i,k)))-
     .               pgeoh(i,k))/RCPD
         ptd(i,k) = MIN(400.,ptd(i,k))
         ptd(i,k) = MAX(100.,ptd(i,k))
         zcond(i) = pqd(i,k)
      ENDIF
      ENDDO
C
      icall = 2
      CALL flxadjtq(paph(1,k), ptd(1,k), pqd(1,k), llo2, icall)
C
      DO i = 1, klon
      IF (llo2(i)) THEN
         zcond(i) = zcond(i)-pqd(i,k)
         zbuo = ptd(i,k)*(1.+RETV  *pqd(i,k))-
     .          ptenh(i,k)*(1.+RETV  *pqenh(i,k))
         llo1 = zbuo.LT.0..AND.(prfl(i)-pmfd(i,k)*zcond(i).GT.0.)
         IF (.not.llo1) pmfd(i,k) = 0.0
         pmfds(i,k) = (RCPD*ptd(i,k)+pgeoh(i,k))*pmfd(i,k)
         pmfdq(i,k) = pqd(i,k)*pmfd(i,k)
         zdmfdp = -pmfd(i,k)*zcond(i)
         pdmfdp(i,k-1) = zdmfdp
         prfl(i) = prfl(i)+zdmfdp
      ENDIF
      ENDDO
c
  180 CONTINUE
      RETURN
      END
      SUBROUTINE flxadjtq(pp, pt, pq, ldflag, kcall)
      use dimens_m
      use dimphy
      use SUPHEC_M
      use yoethf_m
      use fcttre
      IMPLICIT none
c======================================================================
c Objet: ajustement entre T et Q
c======================================================================
C NOTE: INPUT PARAMETER kcall DEFINES CALCULATION AS
C        kcall=0    ENV. T AND QS IN*CUINI*
C        kcall=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
C        kcall=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
C
C
      REAL pt(klon), pq(klon), pp(klon)
      LOGICAL ldflag(klon)
      INTEGER kcall
c
      REAL zcond(klon), zcond1
      REAL Z5alvcp, z5alscp, zalvdcp, zalsdcp
      REAL zdelta, zcvm5, zldcp, zqsat, zcor
      INTEGER is, i
C
      z5alvcp = r5les*RLVTT/RCPD
      z5alscp = r5ies*RLSTT/RCPD
      zalvdcp = rlvtt/RCPD
      zalsdcp = rlstt/RCPD
C

      DO i = 1, klon
         zcond(i) = 0.0
      ENDDO

      DO 210 i =1, klon
      IF (ldflag(i)) THEN
         zdelta = MAX(0.,SIGN(1.,RTT-pt(i)))
         zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
         zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
         zqsat = R2ES*FOEEW(pt(i),zdelta) / pp(i)
         zqsat = MIN(0.5,zqsat)
         zcor = 1./(1.-RETV*zqsat)
         zqsat = zqsat*zcor
         zcond(i) = (pq(i)-zqsat)
     .     / (1. + FOEDE(pt(i), zdelta, zcvm5, zqsat, zcor))
         IF (kcall.EQ.1) zcond(i) = MAX(zcond(i),0.)
         IF (kcall.EQ.2) zcond(i) = MIN(zcond(i),0.)
         pt(i) = pt(i) + zldcp*zcond(i)
         pq(i) = pq(i) - zcond(i)
      ENDIF
  210 CONTINUE
C
      is = 0
      DO i =1, klon
         IF (zcond(i).NE.0.) is = is + 1
      ENDDO
      IF (is.EQ.0) GOTO 230
C
      DO 220 i = 1, klon
      IF(ldflag(i).AND.zcond(i).NE.0.) THEN
         zdelta = MAX(0.,SIGN(1.,RTT-pt(i)))
         zcvm5 = z5alvcp*(1.-zdelta) + zdelta*z5alscp
         zldcp = zalvdcp*(1.-zdelta) + zdelta*zalsdcp
         zqsat = R2ES* FOEEW(pt(i),zdelta) / pp(i)
         zqsat = MIN(0.5,zqsat)
         zcor = 1./(1.-RETV*zqsat)
         zqsat = zqsat*zcor
         zcond1 = (pq(i)-zqsat)
     .     / (1. + FOEDE(pt(i),zdelta,zcvm5,zqsat,zcor))
         pt(i) = pt(i) + zldcp*zcond1
         pq(i) = pq(i) - zcond1
      ENDIF
  220 CONTINUE
C
  230 CONTINUE
      RETURN
      END
      SUBROUTINE flxsetup
            use yoecumf
      IMPLICIT none
C
C     THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME
C
C
      ENTRPEN=1.0E-4  ! ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
      ENTRSCV=3.0E-4  ! ENTRAINMENT RATE FOR SHALLOW CONVECTION
      ENTRMID=1.0E-4  ! ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
      ENTRDD =2.0E-4  ! ENTRAINMENT RATE FOR DOWNDRAFTS
      CMFCTOP=0.33  ! RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUO LEVEL
      CMFCMAX=1.0  ! MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
      CMFCMIN=1.E-10  ! MINIMUM MASSFLUX VALUE (FOR SAFETY)
      CMFDEPS=0.3  ! FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
      CPRCON =2.0E-4  ! CONVERSION FROM CLOUD WATER TO RAIN
      RHCDD=1.  ! RELATIVE SATURATION IN DOWNDRAFRS (NO LONGER USED)
c                 (FORMULATION IMPLIES SATURATION)
      LMFPEN = .TRUE.
      LMFSCV = .TRUE.
      LMFMID = .TRUE.
      LMFDD = .TRUE.
      LMFDUDV = .TRUE.
c
      RETURN
      END
