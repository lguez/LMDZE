module orodrag_m

  IMPLICIT NONE

contains

  SUBROUTINE orodrag(nlon, nlev, ktest, ptsphy, paphm1, papm1, pgeom1, ptm1, &
       pum1, pvm1, pmea, pstd, psig, pgamma, ptheta, ppic, pval, pulow, &
       pvlow, pvom, pvol, pte)

    USE dimens_m
    USE dimphy
    use gwstress_m, only: gwstress
    USE suphec_m
    USE yoegwd
    use gwprofil_m, only: gwprofil
    use orosetup_m, only: orosetup

    !**** *gwdrag* - does the gravity wave parametrization.

    ! purpose.

    ! this routine computes the physical tendencies of the
    ! prognostic variables u, v and t due to vertical transports by
    ! subgridscale orographically excited gravity waves

    !** interface.

    ! called from *callpar*.

    ! the routine takes its input from the long-term storage:
    ! u, v, t and p at t-1.

    ! explicit arguments :

    ! ==== inputs ===
    ! ==== outputs ===

    ! implicit arguments : none

    ! implicit logical (l)

    ! method.

    ! reference.

    ! author.

    ! m.miller + b.ritter e.c.m.w.f. 15/06/86.

    ! f.lott + m. miller e.c.m.w.f. 22/11/94

    !* 0.1 arguments

    INTEGER nlon, nlev
    INTEGER jl, ilevp1, jk, ji
    REAL zdelp, ztemp, zforc, ztend
    REAL rover, zb, zc, zconb, zabsv
    REAL zzd1, ratio, zbet, zust, zvst, zdis
    REAL pte(nlon, nlev), pvol(nlon, nlev), pvom(nlon, nlev), pulow(klon), &
         pvlow(klon)
    REAL pum1(nlon, nlev), pvm1(nlon, nlev), ptm1(nlon, nlev), pmea(nlon)
    REAL, INTENT (IN) :: pstd(nlon)
    REAL, INTENT (IN) :: psig(nlon)
    REAL pgamma(nlon), ptheta(nlon), ppic(nlon), pval(nlon), &
         pgeom1(nlon, nlev), papm1(nlon, nlev), paphm1(nlon, nlev+1)

    INTEGER ktest(nlon)

    !* 0.2 local arrays

    INTEGER icrit(klon), ikcrith(klon), ikenvh(klon), &
         iknu(klon), iknu2(klon), ikcrit(klon)

    REAL ztau(klon, klev+1), zstab(klon, klev+1), &
         zvph(klon, klev+1), zrho(klon, klev+1), zri(klon, klev+1), &
         zpsi(klon, klev+1), zzdep(klon, klev)
    REAL zdudt(klon), zdvdt(klon), zvidis(klon), &
         znu(klon), zd1(klon), zd2(klon), zdmod(klon)
    REAL ztmst
    REAL, INTENT (IN) :: ptsphy

    !------------------------------------------------------------------

    !* 1. initialization

    !* 1.1 computational constants

    ztmst = ptsphy

    !* 1.3 check whether row contains point for printing

    !* 2. precompute basic state variables.

    !* define low level wind, project winds in plane of
    !* low level wind, determine sector in which to take
    !* the variance and set indicator for critical levels.

    CALL orosetup(nlon, ktest, ikcrit, ikcrith, icrit, ikenvh, iknu, iknu2, &
         paphm1, papm1, pum1, pvm1, ptm1, pgeom1, zrho, zri, zstab, ztau, &
         zvph, zpsi, zzdep, pulow, pvlow, ptheta, pgamma, pmea, ppic, pval, &
         znu, zd1, zd2, zdmod)

    !* 3. compute low level stresses using subcritical and
    !* supercritical forms.computes anisotropy coefficient
    !* as measure of orographic twodimensionality.

    CALL gwstress(nlon, nlev, ktest, ikenvh, zrho, zstab, zvph, pstd, &
         psig, pmea, ppic, ztau, pgeom1, zdmod)

    !* 4. compute stress profile.

    CALL gwprofil(nlon, nlev, ktest, ikcrith, icrit, paphm1, zrho, zstab, &
         zvph, zri, ztau, zdmod, psig, pstd)

    !* 5. compute tendencies.

    ! explicit solution at all levels for the gravity wave
    ! implicit solution for the blocked levels

    DO jl = 1, klon
       zvidis(jl) = 0.0
       zdudt(jl) = 0.0
       zdvdt(jl) = 0.0
    end DO

    ilevp1 = klev + 1

    DO jk = 1, klev

       ! Modif vectorisation 02/04/2004
       DO ji = 1, klon
          IF (ktest(ji)==1) THEN

             zdelp = paphm1(ji, jk+1) - paphm1(ji, jk)
             ztemp = -rg*(ztau(ji, jk+1)-ztau(ji, jk))/(zvph(ji, ilevp1)*zdelp)
             zdudt(ji) = (pulow(ji)*zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
             zdvdt(ji) = (pvlow(ji)*zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)

             ! controle des overshoots:

             zforc = sqrt(zdudt(ji)**2+zdvdt(ji)**2) + 1.E-12
             ztend = sqrt(pum1(ji, jk)**2+pvm1(ji, jk)**2)/ztmst + 1.E-12
             rover = 0.25
             IF (zforc>=rover*ztend) THEN
                zdudt(ji) = rover*ztend/zforc*zdudt(ji)
                zdvdt(ji) = rover*ztend/zforc*zdvdt(ji)
             END IF

             ! fin du controle des overshoots

             IF (jk>=ikenvh(ji)) THEN
                zb = 1.0 - 0.18*pgamma(ji) - 0.04*pgamma(ji)**2
                zc = 0.48*pgamma(ji) + 0.3*pgamma(ji)**2
                zconb = 2.*ztmst*gkwake*psig(ji)/(4.*pstd(ji))
                zabsv = sqrt(pum1(ji, jk)**2+pvm1(ji, jk)**2)/2.
                zzd1 = zb*cos(zpsi(ji, jk))**2 + zc*sin(zpsi(ji, jk))**2
                ratio = (cos(zpsi(ji, jk))**2+pgamma(ji)*sin(zpsi(ji, &
                     jk))**2)/(pgamma(ji)*cos(zpsi(ji, jk))**2+sin(zpsi(ji, jk))**2)
                zbet = max(0., 2.-1./ratio)*zconb*zzdep(ji, jk)*zzd1*zabsv

                ! simplement oppose au vent

                zdudt(ji) = -pum1(ji, jk)/ztmst
                zdvdt(ji) = -pvm1(ji, jk)/ztmst

                ! projection dans la direction de l'axe principal de l'orographie
                !mod zdudt(ji)=-(pum1(ji, jk)*cos(ptheta(ji)*rpi/180.)
                !mod * +pvm1(ji, jk)*sin(ptheta(ji)*rpi/180.))
                !mod * *cos(ptheta(ji)*rpi/180.)/ztmst
                !mod zdvdt(ji)=-(pum1(ji, jk)*cos(ptheta(ji)*rpi/180.)
                !mod * +pvm1(ji, jk)*sin(ptheta(ji)*rpi/180.))
                !mod * *sin(ptheta(ji)*rpi/180.)/ztmst
                zdudt(ji) = zdudt(ji)*(zbet/(1.+zbet))
                zdvdt(ji) = zdvdt(ji)*(zbet/(1.+zbet))
             END IF
             pvom(ji, jk) = zdudt(ji)
             pvol(ji, jk) = zdvdt(ji)
             zust = pum1(ji, jk) + ztmst*zdudt(ji)
             zvst = pvm1(ji, jk) + ztmst*zdvdt(ji)
             zdis = 0.5*(pum1(ji, jk)**2+pvm1(ji, jk)**2-zust**2-zvst**2)
             zvidis(ji) = zvidis(ji) + zdis*zdelp

             ! ENCORE UN TRUC POUR EVITER LES EXPLOSIONS

             pte(ji, jk) = 0.0

          END IF
       end DO

    end DO

    RETURN
  END SUBROUTINE orodrag

end module orodrag_m
