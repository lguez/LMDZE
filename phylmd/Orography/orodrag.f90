module orodrag_m

  IMPLICIT NONE

contains

  SUBROUTINE orodrag(ptsphy, paphm1, papm1, pgeom1, ptm1, pum1, pvm1, zmea, &
       zstd, zsig, zgam, zthe, zpic, zval, pulow, pvlow, pvom, pvol, pte)

    USE dimphy, only: klon, klev
    use gwstress_m, only: gwstress
    USE suphec_m, only: rg
    USE yoegwd, only: gkwake
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

    INTEGER jl, ilevp1, jk, ji
    REAL zdelp, ztemp, zforc, ztend
    REAL rover, zb, zc, zconb, zabsv
    REAL zzd1, ratio, zbet, zust, zvst, zdis
    REAL pte(klon, klev), pvol(klon, klev), pvom(klon, klev)
    real, intent(out):: pulow(klon), pvlow(klon)
    REAL, INTENT(IN):: pum1(klon, klev), pvm1(klon, klev), ptm1(klon, klev), &
         zmea(klon)
    REAL, INTENT(IN):: zstd(klon)
    REAL, INTENT(IN):: zsig(klon)
    REAL, intent(inout):: zgam(klon)
    real, INTENT(IN):: zthe(klon), zpic(klon), zval(klon)
    real pgeom1(klon, klev), papm1(klon, klev), paphm1(klon, klev+1)

    !* 0.2 local arrays

    logical ktest(klon) ! points pour lesquels le sch\'ema est actif
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

    ktest = zpic - zmea > 100. .AND. zstd > 10.

    !* 1.1 computational constants

    ztmst = ptsphy

    !* 1.3 check whether row contains point for printing

    !* 2. precompute basic state variables.

    !* define low level wind, project winds in plane of
    !* low level wind, determine sector in which to take
    !* the variance and set indicator for critical levels.

    CALL orosetup(klon, ktest, ikcrit, ikcrith, icrit, ikenvh, iknu, iknu2, &
         paphm1, papm1, pum1, pvm1, ptm1, pgeom1, zrho, zri, zstab, ztau, &
         zvph, zpsi, zzdep, pulow, pvlow, zthe, zgam, zmea, zpic, zval, &
         znu, zd1, zd2, zdmod)

    !* 3. compute low level stresses using subcritical and
    !* supercritical forms.computes anisotropy coefficient
    !* as measure of orographic twodimensionality.

    CALL gwstress(klon, klev, ktest, ikenvh, zrho, zstab, zvph, zstd, &
         zsig, zmea, zpic, ztau, pgeom1, zdmod)

    !* 4. compute stress profile.

    CALL gwprofil(klon, klev, ktest, ikcrith, icrit, paphm1, zrho, zstab, &
         zvph, zri, ztau, zdmod, zsig, zstd)

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
          IF (ktest(ji)) THEN

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
                zb = 1.0 - 0.18*zgam(ji) - 0.04*zgam(ji)**2
                zc = 0.48*zgam(ji) + 0.3*zgam(ji)**2
                zconb = 2.*ztmst*gkwake*zsig(ji)/(4.*zstd(ji))
                zabsv = sqrt(pum1(ji, jk)**2+pvm1(ji, jk)**2)/2.
                zzd1 = zb*cos(zpsi(ji, jk))**2 + zc*sin(zpsi(ji, jk))**2
                ratio = (cos(zpsi(ji, jk))**2+zgam(ji)*sin(zpsi(ji, &
                     jk))**2)/(zgam(ji)*cos(zpsi(ji, jk))**2+sin(zpsi(ji, jk))**2)
                zbet = max(0., 2.-1./ratio)*zconb*zzdep(ji, jk)*zzd1*zabsv

                ! simplement oppose au vent

                zdudt(ji) = -pum1(ji, jk)/ztmst
                zdvdt(ji) = -pvm1(ji, jk)/ztmst

                ! projection dans la direction de l'axe principal de
                ! l'orographie
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
