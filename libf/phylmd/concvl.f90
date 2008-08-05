SUBROUTINE concvl(iflag_con, dtime, paprs, pplay, t, q, u, v, tra,&
     ntra, work1, work2, d_t, d_q, d_u, d_v, d_tra, rain, snow, kbas,&
     ktop, upwd, dnwd, dnwdbis, ma, cape, tvp, iflag, pbase, bbase,&
     dtvpdt1, dtvpdq1, dplcldt, dplcldr, qcondc, wd, pmflxr, pmflxs,&
     da, phi, mp)

  ! From phylmd/concvl.F, v 1.3 2005/04/15 12:36:17
  ! Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: schema de convection de Emanuel (1991) interface

  USE dimens_m
  USE dimphy
  USE yomcst
  USE yoethf
  USE fcttre

  IMPLICIT NONE

  ! Arguments:
  ! dtime--input-R-pas d'integration (s)
  ! s-------input-R-la valeur "s" pour chaque couche
  ! sigs----input-R-la valeur "sigma" de chaque couche
  ! sig-----input-R-la valeur de "sigma" pour chaque niveau
  ! psolpa--input-R-la pression au sol (en Pa)
  ! pskapa--input-R-exponentiel kappa de psolpa
  ! h-------input-R-enthalpie potentielle (Cp*T/P**kappa)
  ! q-------input-R-vapeur d'eau (en kg/kg)

  ! work*: input et output: deux variables de travail,
  !                            on peut les mettre a 0 au debut
  ! ALE-----input-R-energie disponible pour soulevement

  ! d_h-----output-R-increment de l'enthalpie potentielle (h)
  ! d_q-----output-R-increment de la vapeur d'eau
  ! rain----output-R-la pluie (mm/s)
  ! snow----output-R-la neige (mm/s)
  ! upwd----output-R-saturated updraft mass flux (kg/m**2/s)
  ! dnwd----output-R-saturated downdraft mass flux (kg/m**2/s)
  ! dnwd0---output-R-unsaturated downdraft mass flux (kg/m**2/s)
  ! Cape----output-R-CAPE (J/kg)
  ! Tvp-----output-R-Temperature virtuelle d'une parcelle soulevee
  !                  adiabatiquement a partir du niveau 1 (K)
  ! deltapb-output-R-distance entre LCL et base de la colonne (<0 ;
  !  Pa)
  ! Ice_flag-input-L-TRUE->prise en compte de la thermodynamique de
  !  la glace

  INTEGER ntrac
  PARAMETER (ntrac=nqmx-2)

  INTEGER, INTENT (IN) :: iflag_con

  REAL, INTENT (IN) :: dtime
  REAL, INTENT (IN) :: paprs(klon, klev+1)
  REAL, INTENT (IN) :: pplay(klon, klev)
  REAL t(klon, klev), q(klon, klev), u(klon, klev), v(klon, klev)
  REAL, INTENT (IN):: tra(klon, klev, ntrac)
  INTEGER ntra
  REAL work1(klon, klev), work2(klon, klev)
  REAL pmflxr(klon, klev+1), pmflxs(klon, klev+1)

  REAL d_t(klon, klev), d_q(klon, klev), d_u(klon, klev), d_v(klon,&
       klev)
  REAL d_tra(klon, klev, ntrac)
  REAL rain(klon), snow(klon)

  INTEGER kbas(klon), ktop(klon)
  REAL em_ph(klon, klev+1), em_p(klon, klev)
  REAL upwd(klon, klev), dnwd(klon, klev), dnwdbis(klon, klev)
  REAL ma(klon, klev), cape(klon), tvp(klon, klev)
  REAL da(klon, klev), phi(klon, klev, klev), mp(klon, klev)
  INTEGER iflag(klon)
  REAL pbase(klon), bbase(klon)
  REAL dtvpdt1(klon, klev), dtvpdq1(klon, klev)
  REAL dplcldt(klon), dplcldr(klon)
  REAL qcondc(klon, klev)
  REAL wd(klon)

  REAL zx_t, zdelta, zx_qs, zcor

  INTEGER i, k, itra
  REAL qs(klon, klev)
  REAL cbmf(klon)
  SAVE cbmf
  INTEGER ifrst
  SAVE ifrst
  DATA ifrst/0/

  !-----------------------------------------------------------------

  snow(:) = 0

  IF (ifrst==0) THEN
     ifrst = 1
     DO i = 1, klon
        cbmf(i) = 0.
     END DO
  END IF

  DO k = 1, klev + 1
     DO i = 1, klon
        em_ph(i, k) = paprs(i, k)/100.0
        pmflxs(i, k) = 0.
     END DO
  END DO

  DO k = 1, klev
     DO i = 1, klon
        em_p(i, k) = pplay(i, k)/100.0
     END DO
  END DO


  IF (iflag_con==4) THEN
     DO k = 1, klev
        DO i = 1, klon
           zx_t = t(i, k)
           zdelta = max(0., sign(1., rtt-zx_t))
           zx_qs = min(0.5, r2es*foeew(zx_t, zdelta)/em_p(i, k)/100.0)
           zcor = 1./(1.-retv*zx_qs)
           qs(i, k) = zx_qs*zcor
        END DO
     END DO
  ELSE
     ! iflag_con=3 (modif de puristes qui fait la diffce pour la
     ! convergence numerique)
     DO k = 1, klev
        DO i = 1, klon
           zx_t = t(i, k)
           zdelta = max(0., sign(1., rtt-zx_t))
           zx_qs = r2es*foeew(zx_t, zdelta)/em_p(i, k)/100.0
           zx_qs = min(0.5, zx_qs)
           zcor = 1./(1.-retv*zx_qs)
           zx_qs = zx_qs*zcor
           qs(i, k) = zx_qs
        END DO
     END DO
  END IF

  ! Main driver for convection:
  !		iflag_con = 3  -> equivalent to convect3
  !		iflag_con = 4  -> equivalent to convect1/2

  CALL cv_driver(klon, klev, klev+1, ntra, iflag_con, t, q, qs, u, v,&
       tra, em_p, em_ph, iflag, d_t, d_q, d_u, d_v, d_tra, rain,&
       pmflxr, cbmf, work1, work2, kbas, ktop, dtime, ma, upwd, dnwd,&
       dnwdbis, qcondc, wd, cape, da, phi, mp)

  DO i = 1, klon
     rain(i) = rain(i)/86400.
  END DO

  DO k = 1, klev
     DO i = 1, klon
        d_t(i, k) = dtime*d_t(i, k)
        d_q(i, k) = dtime*d_q(i, k)
        d_u(i, k) = dtime*d_u(i, k)
        d_v(i, k) = dtime*d_v(i, k)
     END DO
  END DO
  DO itra = 1, ntra
     DO k = 1, klev
        DO i = 1, klon
           d_tra(i, k, itra) = dtime*d_tra(i, k, itra)
        END DO
     END DO
  END DO
  ! les traceurs ne sont pas mis dans cette version de convect4:
  IF (iflag_con==4) THEN
     DO itra = 1, ntra
        DO k = 1, klev
           DO i = 1, klon
              d_tra(i, k, itra) = 0.
           END DO
        END DO
     END DO
  END IF

END SUBROUTINE concvl
