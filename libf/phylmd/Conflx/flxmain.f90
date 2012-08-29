module flxmain_m

  IMPLICIT none

contains

  SUBROUTINE flxmain(pdtime, pten, pqen, pqsen, pqhfl, pap, paph, pgeo, &
       ldland, ptte, pqte, pvervel, prsfc, pssfc, kcbot, kctop, kdtop, pmfu, &
       pmfd, pen_u, pde_u, pen_d, pde_d, dt_con, dq_con, pmflxr, pmflxs)

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rcpd, retv, rg, rlvtt
    USE yoethf_m, ONLY: r4les, r5les
    USE yoecumf, ONLY: flxsetup, cmfdeps, entrpen, entrscv, lmfdd

    REAL, intent(in):: pdtime
    REAL pten(klon,klev), pqen(klon,klev), pqsen(klon,klev)
    REAL ptte(klon,klev)
    REAL pqte(klon,klev)
    REAL pvervel(klon,klev)
    REAL pgeo(klon,klev), pap(klon,klev), paph(klon,klev+1)
    REAL pqhfl(klon)

    REAL ptu(klon,klev), pqu(klon,klev), plu(klon,klev)
    REAL plude(klon,klev)
    REAL pmfu(klon,klev)
    REAL prsfc(klon), pssfc(klon)
    INTEGER kcbot(klon), kctop(klon), ktype(klon)
    LOGICAL ldland(klon), ldcum(klon)

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
    INTEGER ilab(klon,klev), ictop0(klon)
    LOGICAL llo1
    REAL dt_con(klon,klev), dq_con(klon,klev)
    REAL zmfmax, zdh
    real zqumqe, zdqmin, zalvdcp, zhsat, zzz
    REAL zhhat, zpbmpt, zgam, zeps, zfac
    INTEGER i, k, ikb, itopm2, kcum

    REAL pen_u(klon,klev), pde_u(klon,klev)
    REAL pen_d(klon,klev), pde_d(klon,klev)

    REAL ptd(klon,klev), pqd(klon,klev), pmfd(klon,klev)
    REAL zmfds(klon,klev), zmfdq(klon,klev), zdmfdp(klon,klev)
    INTEGER kdtop(klon)
    LOGICAL lddraf(klon)

    LOGICAL:: firstcal = .TRUE.

    !---------------------------------------------------------------------

    IF (firstcal) THEN
       CALL flxsetup
       firstcal = .FALSE.
    ENDIF

    DO i = 1, klon
       ldcum(i) = .FALSE.
    ENDDO
    DO k = 1, klev
       DO i = 1, klon
          dt_con(i,k) = 0.0
          dq_con(i,k) = 0.0
       ENDDO
    ENDDO

    ! initialiser les variables et faire l'interpolation verticale

    CALL flxini(pten, pqen, pqsen, pgeo, &
         paph, zgeoh, ztenh, zqenh, zqsenh, &
         ptu, pqu, ptd, pqd, pmfd, zmfds, zmfdq, zdmfdp, &
         pmfu, zmfus, zmfuq, zdmfup, &
         zdpmel, plu, plude, ilab, pen_u, pde_u, pen_d, pde_d)

    ! determiner les valeurs au niveau de base de la tour convective

    CALL flxbase(ztenh, zqenh, zgeoh, paph, &
         ptu, pqu, plu, ldcum, kcbot, ilab)

    ! calculer la convergence totale de l'humidite et celle en provenance
    ! de la couche limite, plus precisement, la convergence integree entre
    ! le sol et la base de la convection. Cette derniere convergence est
    ! comparee avec l'evaporation obtenue dans la couche limite pour
    ! determiner le type de la convection

    k=1
    DO i = 1, klon
       zdqcv(i) = pqte(i,k)*(paph(i,k+1)-paph(i,k))
       zdhpbl(i) = 0.0
       zdqpbl(i) = 0.0
    ENDDO

    DO k=2,klev
       DO i = 1, klon
          zdqcv(i)=zdqcv(i)+pqte(i,k)*(paph(i,k+1)-paph(i,k))
          IF (k.GE.kcbot(i)) THEN
             zdqpbl(i)=zdqpbl(i)+pqte(i,k)*(paph(i,k+1)-paph(i,k))
             zdhpbl(i)=zdhpbl(i)+(RCPD*ptte(i,k)+RLVTT*pqte(i,k)) &
                  *(paph(i,k+1)-paph(i,k))
          ENDIF
       ENDDO
    ENDDO

    DO i = 1, klon
       ktype(i) = 2
       if (zdqcv(i).GT.MAX(0.,-1.5*pqhfl(i)*RG)) ktype(i) = 1
       !cc if (zdqcv(i).GT.MAX(0.,-1.1*pqhfl(i)*RG)) ktype(i) = 1
    ENDDO

    ! determiner le flux de masse entrant a travers la base.
    ! on ignore, pour l'instant, l'effet du panache descendant

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

    ! DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME

    ! (A) calculer d'abord la hauteur "theorique" de la tour convective sans
    ! considerer l'entrainement ni le detrainement du panache, sachant
    ! ces derniers peuvent abaisser la hauteur theorique.

    DO i = 1, klon
       ikb=kcbot(i)
       zhcbase(i)=RCPD*ptu(i,ikb)+zgeoh(i,ikb)+RLVTT*pqu(i,ikb)
       ictop0(i)=kcbot(i)-1
    ENDDO

    zalvdcp=RLVTT/RCPD
    DO k=klev-1,3,-1
       DO i = 1, klon
          zhsat=RCPD*ztenh(i,k)+zgeoh(i,k)+RLVTT*zqsenh(i,k)
          zgam=R5LES*zalvdcp*zqsenh(i,k)/ &
               ((1.-RETV *zqsenh(i,k))*(ztenh(i,k)-R4LES)**2)
          zzz=RCPD*ztenh(i,k)*0.608
          zhhat=zhsat-(zzz+zgam*zzz)/(1.+zgam*zzz/RLVTT)* &
               MAX(zqsenh(i,k)-zqenh(i,k),0.)
          IF(k.LT.ictop0(i).AND.zhcbase(i).GT.zhhat) ictop0(i)=k
       ENDDO
    ENDDO

    ! (B) calculer le panache ascendant

    CALL flxasc(pdtime,ztenh, zqenh, pten, pqen, pqsen, &
         pgeo, zgeoh, pap, paph, pqte, pvervel, &
         ldland, ldcum, ktype, ilab, &
         ptu, pqu, plu, pmfu, zmfub, zentr, &
         zmfus, zmfuq, zmful, plude, zdmfup, &
         kcbot, kctop, ictop0, kcum, pen_u, pde_u)

    IF (kcum /= 0) then
       ! verifier l'epaisseur de la convection et changer eventuellement
       ! le taux d'entrainement/detrainement

       DO i = 1, klon
          zpbmpt=paph(i,kcbot(i))-paph(i,kctop(i))
          IF(ldcum(i).AND.ktype(i).EQ.1.AND.zpbmpt.LT.2.E4)ktype(i)=2
          IF(ldcum(i)) ictop0(i)=kctop(i)
          IF(ktype(i).EQ.2) zentr(i)=ENTRSCV
       ENDDO

       IF (lmfdd) THEN ! si l'on considere le panache descendant
          ! calculer la precipitation issue du panache ascendant pour
          ! determiner l'existence du panache descendant dans la convection
          DO i = 1, klon
             zrfl(i)=zdmfup(i,1)
          ENDDO
          DO k=2,klev
             DO i = 1, klon
                zrfl(i)=zrfl(i)+zdmfup(i,k)
             ENDDO
          ENDDO

          ! determiner le LFS (level of free sinking: niveau de plonge libre)
          CALL flxdlfs(ztenh, zqenh, zgeoh, paph, ptu, pqu, &
               ldcum, kcbot, kctop, zmfub, zrfl, &
               ptd, pqd, &
               pmfd, zmfds, zmfdq, zdmfdp, &
               kdtop, lddraf)

          ! calculer le panache descendant
          CALL flxddraf(ztenh, zqenh, &
               zgeoh, paph, zrfl, &
               ptd, pqd, &
               pmfd, zmfds, zmfdq, zdmfdp, &
               lddraf, pen_d, pde_d)

          ! calculer de nouveau le flux de masse entrant a travers la base
          ! de la convection, sachant qu'il a ete modifie par le panache
          ! descendant
          DO i = 1, klon
             IF (lddraf(i)) THEN
                ikb = kcbot(i)
                llo1 = PMFD(i,ikb).LT.0.
                zeps = 0.
                IF ( llo1 ) zeps = CMFDEPS
                zqumqe = pqu(i,ikb)+plu(i,ikb)- &
                     zeps*pqd(i,ikb)-(1.-zeps)*zqenh(i,ikb)
                zdqmin = MAX(0.01*zqenh(i,ikb),1.E-10)
                zmfmax = (paph(i,ikb)-paph(i,ikb-1)) / (RG*pdtime)
                IF (zdqpbl(i).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(i) &
                     .AND.zmfub(i).LT.zmfmax) THEN
                   zmfub1(i) = zdqpbl(i) / (RG*MAX(zqumqe,zdqmin))
                ELSE
                   zmfub1(i) = zmfub(i)
                ENDIF
                IF (ktype(i).EQ.2) THEN
                   zdh = RCPD*(ptu(i,ikb)-zeps*ptd(i,ikb)- &
                        (1.-zeps)*ztenh(i,ikb))+RLVTT*zqumqe
                   zdh = RG * MAX(zdh,1.0E5*zdqmin)
                   IF (zdhpbl(i).GT.0..AND.ldcum(i))zmfub1(i)=zdhpbl(i)/zdh
                ENDIF
                IF ( .NOT.((ktype(i).EQ.1.OR.ktype(i).EQ.2).AND. &
                     ABS(zmfub1(i)-zmfub(i)).LT.0.2*zmfub(i)) ) &
                     zmfub1(i) = zmfub(i)
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
       ENDIF ! fin de test sur lmfdd

       ! calculer de nouveau le panache ascendant

       CALL flxasc(pdtime,ztenh, zqenh, pten, pqen, pqsen, &
            pgeo, zgeoh, pap, paph, pqte, pvervel, &
            ldland, ldcum, ktype, ilab, &
            ptu, pqu, plu, pmfu, zmfub, zentr, &
            zmfus, zmfuq, zmful, plude, zdmfup, &
            kcbot, kctop, ictop0, kcum, pen_u, pde_u)

       ! determiner les flux convectifs en forme finale, ainsi que
       ! la quantite des precipitations

       CALL flxflux(pdtime, pqen, pqsen, ztenh, zqenh, pap, paph, &
            ldland, zgeoh, kcbot, kctop, lddraf, kdtop, ktype, ldcum, &
            pmfu, pmfd, zmfus, zmfds, zmfuq, zmfdq, zmful, plude, &
            zdmfup, zdmfdp, pten, prsfc, pssfc, zdpmel, itopm2, &
            pmflxr, pmflxs)

       ! calculer les tendances pour T et Q

       CALL flxdtdq(itopm2, paph, ldcum, pten, &
            zmfus, zmfds, zmfuq, zmfdq, zmful, zdmfup, zdmfdp, zdpmel, &
            dt_con,dq_con)
    end IF

  END SUBROUTINE flxmain

end module flxmain_m
