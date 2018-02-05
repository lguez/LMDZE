module orosetup_m

  IMPLICIT NONE

contains

  SUBROUTINE orosetup(nlon, ktest, kkcrit, kkcrith, kcrit, kkenvh, kknu, &
       kknu2, paphm1, papm1, pum1, pvm1, ptm1, pgeom1, prho, pri, pstab, ptau, &
       pvph, ppsi, pzdep, pulow, pvlow, ptheta, pgamma, pmea, ppic, pval, pnu, &
       pd1, pd2, pdmod)

    ! gwsetup
    ! Interface from orodrag.
    ! See ECMWF research department documentation of the I.F.S.
    ! Modifications F. Lott for the new-gwdrag scheme, november 1993.

    USE dimens_m
    USE dimphy
    use nr_util, only: pi
    USE suphec_m
    USE yoegwd

    ! 0.1   arguments

    INTEGER nlon
    INTEGER jl, jk
    REAL zdelp

    INTEGER kkcrit(nlon), kkcrith(nlon), kcrit(nlon), ktest(nlon), &
         kkenvh(nlon)

    REAL paphm1(nlon, klev+1), papm1(nlon, klev), pum1(nlon, klev), &
         pvm1(nlon, klev), ptm1(nlon, klev), pgeom1(nlon, klev), &
         prho(nlon, klev+1), pri(nlon, klev+1), pstab(nlon, klev+1), &
         ptau(nlon, klev+1), pvph(nlon, klev+1), ppsi(nlon, klev+1), &
         pzdep(nlon, klev)
    REAL pulow(nlon), pvlow(nlon), ptheta(nlon), pgamma(nlon), pnu(nlon), &
         pd1(nlon), pd2(nlon), pdmod(nlon)
    REAL pmea(nlon), ppic(nlon), pval(nlon)

    ! 0.2   local arrays

    INTEGER ilevh
    REAL zcons1, zcons2, zhgeo
    REAL zu, zphi, zvt1, zvt2, zst, zdwind, zwind
    REAL zstabm, zstabp, zrhom, zrhop
    LOGICAL lo
    LOGICAL ll1(klon, klev+1)
    INTEGER kknu(klon), kknu2(klon), kknub(klon), kknul(klon)

    REAL zhcrit(klon, klev), zvpf(klon, klev), zdp(klon, klev)
    REAL znorm(klon), zb(klon), zc(klon), znup(klon), znum(klon)

    !------------------------------------------------------------------

    !!print *, "Call sequence information: orosetup"
    ! 1.    initialization
    ! 1.1   computational constants

    ilevh = klev/3

    zcons1 = 1./rd
    !old  zcons2=g**2/cpd
    zcons2 = rg**2/rcpd

    ! 2.

    ! 2.1     define low level wind, project winds in plane of
    ! low level wind, determine sector in which to take
    ! the variance and set indicator for critical levels.

    DO  jl = 1, klon
       kknu(jl) = klev
       kknu2(jl) = klev
       kknub(jl) = klev
       kknul(jl) = klev
       pgamma(jl) = max(pgamma(jl), gtsec)
       ll1(jl, klev+1) = .FALSE.
    end DO

    ! Ajouter une initialisation (L. Li, le 23fev99):

    DO jk = klev, ilevh, -1
       DO jl = 1, klon
          ll1(jl, jk) = .FALSE.
       END DO
    END DO

    ! define top of low level flow

    DO  jk = klev, ilevh, -1
       DO  jl = 1, klon
          lo = (paphm1(jl, jk)/paphm1(jl, klev+1)) >= gsigcr
          IF (lo) THEN
             kkcrit(jl) = jk
          END IF
          zhcrit(jl, jk) = ppic(jl)
          zhgeo = pgeom1(jl, jk)/rg
          ll1(jl, jk) = (zhgeo>zhcrit(jl, jk))
          IF (ll1(jl, jk) .NEQV. ll1(jl, jk+1)) THEN
             kknu(jl) = jk
          END IF
          IF ( .NOT. ll1(jl, ilevh)) kknu(jl) = ilevh
       end DO
    end DO
    DO  jk = klev, ilevh, -1
       DO  jl = 1, klon
          zhcrit(jl, jk) = ppic(jl) - pval(jl)
          zhgeo = pgeom1(jl, jk)/rg
          ll1(jl, jk) = (zhgeo>zhcrit(jl, jk))
          IF (ll1(jl, jk) .NEQV. ll1(jl, jk+1)) THEN
             kknu2(jl) = jk
          END IF
          IF ( .NOT. ll1(jl, ilevh)) kknu2(jl) = ilevh
       end DO
    end DO
    DO  jk = klev, ilevh, -1
       DO  jl = 1, klon
          zhcrit(jl, jk) = amax1(ppic(jl)-pmea(jl), pmea(jl)-pval(jl))
          zhgeo = pgeom1(jl, jk)/rg
          ll1(jl, jk) = (zhgeo>zhcrit(jl, jk))
          IF (ll1(jl, jk) .NEQV. ll1(jl, jk+1)) THEN
             kknub(jl) = jk
          END IF
          IF ( .NOT. ll1(jl, ilevh)) kknub(jl) = ilevh
       end DO
    end DO

    DO  jl = 1, klon
       kknu(jl) = min(kknu(jl), nktopg)
       kknu2(jl) = min(kknu2(jl), nktopg)
       kknub(jl) = min(kknub(jl), nktopg)
       kknul(jl) = klev
    end DO

    !c*     initialize various arrays

    DO jl = 1, klon
       prho(jl, klev+1) = 0.0
       pstab(jl, klev+1) = 0.0
       pstab(jl, 1) = 0.0
       pri(jl, klev+1) = 9999.0
       ppsi(jl, klev+1) = 0.0
       pri(jl, 1) = 0.0
       pvph(jl, 1) = 0.0
       pulow(jl) = 0.0
       pvlow(jl) = 0.0
       kkcrith(jl) = klev
       kkenvh(jl) = klev
       kcrit(jl) = 1
       ll1(jl, klev+1) = .FALSE.
    end DO

    ! define low-level flow

    DO  jk = klev, 2, -1
       DO  jl = 1, klon
          IF (ktest(jl)==1) THEN
             zdp(jl, jk) = papm1(jl, jk) - papm1(jl, jk-1)
             prho(jl, jk) = 2. * paphm1(jl, jk) * zcons1 &
                  / (ptm1(jl, jk) + ptm1(jl, jk-1))
             pstab(jl, jk) = 2. * zcons2 / (ptm1(jl, jk) + ptm1(jl, jk-1)) &
                  * (1. - rcpd * prho(jl, jk) &
                  * (ptm1(jl, jk) - ptm1(jl, jk - 1)) / zdp(jl, jk))
             pstab(jl, jk) = max(pstab(jl, jk), gssec)
          END IF
       end DO
    end DO

    ! define blocked flow

    DO  jk = klev, ilevh, -1
       DO  jl = 1, klon
          IF (jk>=kknub(jl) .AND. jk<=kknul(jl)) THEN
             pulow(jl) = pulow(jl) + pum1(jl, jk)*(paphm1(jl, jk+1)-paphm1(jl, jk) &
                  )
             pvlow(jl) = pvlow(jl) + pvm1(jl, jk)*(paphm1(jl, jk+1)-paphm1(jl, jk) &
                  )
          END IF
       end DO
    end DO
    DO  jl = 1, klon
       pulow(jl) = pulow(jl)/(paphm1(jl, kknul(jl)+1)-paphm1(jl, kknub(jl)))
       pvlow(jl) = pvlow(jl)/(paphm1(jl, kknul(jl)+1)-paphm1(jl, kknub(jl)))
       znorm(jl) = max(sqrt(pulow(jl)**2+pvlow(jl)**2), gvsec)
       pvph(jl, klev+1) = znorm(jl)
    end DO

    !  setup orography axes and define plane of profiles  

    DO  jl = 1, klon
       lo = (pulow(jl)<gvsec) .AND. (pulow(jl)>=-gvsec)
       IF (lo) THEN
          zu = pulow(jl) + 2.*gvsec
       ELSE
          zu = pulow(jl)
       END IF
       zphi = atan(pvlow(jl)/zu)
       ppsi(jl, klev+1) = ptheta(jl)*pi/180. - zphi
       zb(jl) = 1. - 0.18*pgamma(jl) - 0.04*pgamma(jl)**2
       zc(jl) = 0.48*pgamma(jl) + 0.3*pgamma(jl)**2
       pd1(jl) = zb(jl) - (zb(jl)-zc(jl))*(sin(ppsi(jl, klev+1))**2)
       pd2(jl) = (zb(jl)-zc(jl))*sin(ppsi(jl, klev+1))*cos(ppsi(jl, klev+1))
       pdmod(jl) = sqrt(pd1(jl)**2+pd2(jl)**2)
    end DO

    !   define flow in plane of lowlevel stress 

    DO  jk = 1, klev
       DO  jl = 1, klon
          IF (ktest(jl)==1) THEN
             zvt1 = pulow(jl)*pum1(jl, jk) + pvlow(jl)*pvm1(jl, jk)
             zvt2 = -pvlow(jl)*pum1(jl, jk) + pulow(jl)*pvm1(jl, jk)
             zvpf(jl, jk) = (zvt1*pd1(jl)+zvt2*pd2(jl))/(znorm(jl)*pdmod(jl))
          END IF
          ptau(jl, jk) = 0.0
          pzdep(jl, jk) = 0.0
          ppsi(jl, jk) = 0.0
          ll1(jl, jk) = .FALSE.
       end DO
    end DO
    DO  jk = 2, klev
       DO  jl = 1, klon
          IF (ktest(jl)==1) THEN
             zdp(jl, jk) = papm1(jl, jk) - papm1(jl, jk-1)
             pvph(jl, jk) = ((paphm1(jl, jk)-papm1(jl, jk-1))*zvpf(jl, jk)+(papm1( &
                  jl, jk)-paphm1(jl, jk))*zvpf(jl, jk-1))/zdp(jl, jk)
             IF (pvph(jl, jk)<gvsec) THEN
                pvph(jl, jk) = gvsec
                kcrit(jl) = jk
             END IF
          END IF
       end DO
    end DO

    ! 2.2     brunt-vaisala frequency and density at half levels.

    DO  jk = ilevh, klev
       DO  jl = 1, klon
          IF (ktest(jl)==1) THEN
             IF (jk>=(kknub(jl)+1) .AND. jk<=kknul(jl)) THEN
                zst = zcons2/ptm1(jl, jk)*(1.-rcpd*prho(jl, jk)*(ptm1(jl, &
                     jk)-ptm1(jl, jk-1))/zdp(jl, jk))
                pstab(jl, klev+1) = pstab(jl, klev+1) + zst*zdp(jl, jk)
                pstab(jl, klev+1) = max(pstab(jl, klev+1), gssec)
                prho(jl, klev+1) = prho(jl, klev+1) + paphm1(jl, jk)*2.*zdp(jl, jk)* &
                     zcons1/(ptm1(jl, jk)+ptm1(jl, jk-1))
             END IF
          END IF
       end DO
    end DO

    DO  jl = 1, klon
       pstab(jl, klev + 1) = pstab(jl, klev + 1) &
            / (papm1(jl, kknul(jl)) - papm1(jl, kknub(jl)))
       prho(jl, klev + 1) = prho(jl, klev + 1) &
            / (papm1(jl, kknul(jl)) - papm1(jl, kknub(jl)))
    end DO

    ! 2.3     mean flow richardson number.
    ! and critical height for froude layer

    DO  jk = 2, klev
       DO  jl = 1, klon
          IF (ktest(jl)==1) THEN
             zdwind = max(abs(zvpf(jl, jk)-zvpf(jl, jk-1)), gvsec)
             pri(jl, jk) = pstab(jl, jk)*(zdp(jl, jk)/(rg*prho(jl, jk)*zdwind))**2
             pri(jl, jk) = max(pri(jl, jk), grcrit)
          END IF
       end DO
    end do

    ! define top of 'envelope' layer

    DO  jl = 1, klon
       pnu(jl) = 0.0
       znum(jl) = 0.0
    end DO

    DO  jk = 2, klev - 1
       DO jl = 1, klon

          IF (ktest(jl)==1) THEN

             IF (jk>=kknub(jl)) THEN

                znum(jl) = pnu(jl)
                zwind = (pulow(jl)*pum1(jl, jk)+pvlow(jl)*pvm1(jl, jk))/ &
                     max(sqrt(pulow(jl)**2+pvlow(jl)**2), gvsec)
                zwind = max(sqrt(zwind**2), gvsec)
                zdelp = paphm1(jl, jk+1) - paphm1(jl, jk)
                zstabm = sqrt(max(pstab(jl, jk), gssec))
                zstabp = sqrt(max(pstab(jl, jk+1), gssec))
                zrhom = prho(jl, jk)
                zrhop = prho(jl, jk+1)
                pnu(jl) = pnu(jl) + (zdelp/rg)*((zstabp/zrhop+zstabm/zrhom)/2.)/ &
                     zwind
                IF ((znum(jl)<=gfrcrit) .AND. (pnu(jl)>gfrcrit) .AND. (kkenvh( &
                     jl)==klev)) kkenvh(jl) = jk

             END IF

          END IF

       end DO
    end do

    !  calculation of a dynamical mixing height for the breaking
    !  of gravity waves:

    DO  jl = 1, klon
       znup(jl) = 0.0
       znum(jl) = 0.0
    end DO

    DO  jk = klev - 1, 2, -1
       DO  jl = 1, klon

          IF (ktest(jl)==1) THEN

             znum(jl) = znup(jl)
             zwind = (pulow(jl)*pum1(jl, jk)+pvlow(jl)*pvm1(jl, jk))/ &
                  max(sqrt(pulow(jl)**2+pvlow(jl)**2), gvsec)
             zwind = max(sqrt(zwind**2), gvsec)
             zdelp = paphm1(jl, jk+1) - paphm1(jl, jk)
             zstabm = sqrt(max(pstab(jl, jk), gssec))
             zstabp = sqrt(max(pstab(jl, jk+1), gssec))
             zrhom = prho(jl, jk)
             zrhop = prho(jl, jk+1)
             znup(jl) = znup(jl) + (zdelp/rg)*((zstabp/zrhop+zstabm/zrhom)/2.)/ &
                  zwind
             IF ((znum(jl)<=pi/2.) .AND. (znup(jl)>pi/2.) .AND. (kkcrith( &
                  jl)==klev)) kkcrith(jl) = jk

          END IF

       end DO
    end DO

    DO  jl = 1, klon
       kkcrith(jl) = min0(kkcrith(jl), kknu2(jl))
       kkcrith(jl) = max0(kkcrith(jl), ilevh*2)
    end DO

    !     directional info for flow blocking 

    DO  jk = ilevh, klev
       DO  jl = 1, klon
          IF (jk>=kkenvh(jl)) THEN
             lo = (pum1(jl, jk)<gvsec) .AND. (pum1(jl, jk)>=-gvsec)
             IF (lo) THEN
                zu = pum1(jl, jk) + 2.*gvsec
             ELSE
                zu = pum1(jl, jk)
             END IF
             zphi = atan(pvm1(jl, jk)/zu)
             ppsi(jl, jk) = ptheta(jl)*pi/180. - zphi
          END IF
       end DO
    end DO
    !      forms the vertical 'leakiness' 

    DO  jk = ilevh, klev
       DO  jl = 1, klon
          IF (jk>=kkenvh(jl)) THEN
             pzdep(jl, jk) = (pgeom1(jl, kkenvh(jl)-1)-pgeom1(jl, jk))/ &
                  (pgeom1(jl, kkenvh(jl)-1)-pgeom1(jl, klev))
          END IF
       end DO
    end DO

  END SUBROUTINE orosetup

end module orosetup_m
