module flxmain_m

  IMPLICIT none

contains

  SUBROUTINE flxmain(dtime, ten, qen, qsen, pqhfl, pap, paph, pgeo, ldland, &
       ptte, pqte, pvervel, prsfc, pssfc, kcbot, kctop, kdtop, mfu, mfd, &
       pen_u, pde_u, pen_d, pde_d, dt_con, dq_con, pmflxr, pmflxs)

    USE dimphy, ONLY: klev, klon
    use flxasc_m, only: flxasc
    use flxdtdq_m, only: flxdtdq
    use flxflux_m, only: flxflux
    use flxini_m, only: flxini
    USE suphec_m, ONLY: rcpd, retv, rg, rlvtt
    USE yoecumf, ONLY: flxsetup, cmfdeps, entrpen, entrscv, lmfdd
    USE yoethf_m, ONLY: r4les, r5les

    REAL, intent(in):: dtime
    REAL, intent(in):: ten(klon, klev)
    real, intent(in):: qen(klon, klev)
    real, intent(inout):: qsen(klon, klev)
    REAL, intent(in):: pqhfl(klon)
    real pap(klon, klev), paph(klon, klev+1)
    REAL, intent(in):: pgeo(klon, klev)
    LOGICAL ldland(klon)
    REAL ptte(klon, klev)
    REAL pqte(klon, klev)
    REAL pvervel(klon, klev)
    REAL prsfc(klon), pssfc(klon)
    INTEGER kcbot(klon), kctop(klon)
    INTEGER kdtop(klon)
    REAL, intent(out):: mfu(klon, klev)
    real, intent(out):: mfd(klon, klev)
    REAL pen_u(klon, klev), pde_u(klon, klev)
    REAL pen_d(klon, klev), pde_d(klon, klev)
    REAL dt_con(klon, klev), dq_con(klon, klev)
    REAL pmflxr(klon, klev+1)
    REAL pmflxs(klon, klev+1)

    ! Local:
    REAL ptu(klon, klev), pqu(klon, klev), plu(klon, klev)
    REAL plude(klon, klev)
    INTEGER ktype(klon)
    LOGICAL ldcum(klon)

    REAL ztenh(klon, klev), zqenh(klon, klev), zqsenh(klon, klev)
    REAL zgeoh(klon, klev)
    REAL mfub(klon), mfub1(klon)
    REAL mfus(klon, klev), mfuq(klon, klev), mful(klon, klev)
    REAL zdmfup(klon, klev), zdpmel(klon, klev)
    REAL zentr(klon), zhcbase(klon)
    REAL zdqpbl(klon), zdqcv(klon), zdhpbl(klon)
    REAL zrfl(klon)
    INTEGER ilab(klon, klev), ictop0(klon)
    LOGICAL llo1
    REAL zmfmax, zdh
    real zqumqe, zdqmin, zalvdcp, zhsat, zzz
    REAL zhhat, zpbmpt, zgam, zeps, zfac
    INTEGER i, k, ikb, itopm2, kcum


    REAL ptd(klon, klev), pqd(klon, klev)
    REAL zmfds(klon, klev), zmfdq(klon, klev), zdmfdp(klon, klev)
    LOGICAL lddraf(klon)

    LOGICAL:: firstcal = .TRUE.

    !---------------------------------------------------------------------

    IF (firstcal) THEN
       CALL flxsetup
       firstcal = .FALSE.
    ENDIF

    ldcum = .FALSE.
    dt_con = 0.
    dq_con = 0.

    ! Initialiser les variables et faire l'interpolation verticale :
    CALL flxini(ten, qen, qsen, pgeo, paph, zgeoh, ztenh, zqenh, zqsenh, &
         ptu, pqu, ptd, pqd, mfd, zmfds, zmfdq, zdmfdp, mfu, mfus, mfuq, &
         zdmfup, zdpmel, plu, plude, ilab, pen_u, pde_u, pen_d, pde_d)

    ! Déterminer les valeurs au niveau de base de la tour convective :
    CALL flxbase(ztenh, zqenh, zgeoh, paph, ptu, pqu, plu, ldcum, kcbot, ilab)

    ! Calculer la convergence totale de l'humidité et celle en
    ! provenance de la couche limite, plus précisément, la convergence
    ! intégrée entre le sol et la base de la convection. Cette
    ! dernière convergence est comparée avec l'&vaporation obtenue
    ! dans la couche limite pour déterminer le type de la convection.

    zdqcv = pqte(:, 1) * (paph(:, 2) - paph(:, 1))
    zdhpbl = 0.
    zdqpbl = 0.

    DO k=2, klev
       DO i = 1, klon
          zdqcv(i)=zdqcv(i)+pqte(i, k)*(paph(i, k+1)-paph(i, k))
          IF (k.GE.kcbot(i)) THEN
             zdqpbl(i)=zdqpbl(i)+pqte(i, k)*(paph(i, k+1)-paph(i, k))
             zdhpbl(i)=zdhpbl(i)+(RCPD*ptte(i, k)+RLVTT*pqte(i, k)) &
                  *(paph(i, k+1)-paph(i, k))
          ENDIF
       ENDDO
    ENDDO

    DO i = 1, klon
       if (zdqcv(i) > MAX(0., - 1.5 * pqhfl(i) * RG)) then
          ktype(i) = 1
       else
          ktype(i) = 2
       end if
    ENDDO

    ! Déterminer le flux de masse entrant à travers la base. On
    ! ignore, pour l'instant, l'effet du panache descendant

    DO i = 1, klon
       ikb=kcbot(i)
       zqumqe=pqu(i, ikb)+plu(i, ikb)-zqenh(i, ikb)
       zdqmin=MAX(0.01*zqenh(i, ikb), 1.E-10)
       IF (zdqpbl(i) > 0..AND.zqumqe > zdqmin.AND.ldcum(i)) THEN
          mfub(i) = zdqpbl(i)/(RG*MAX(zqumqe, zdqmin))
       ELSE
          mfub(i) = 0.01
          ldcum(i)=.FALSE.
       ENDIF
       IF (ktype(i) == 2) THEN
          zdh = RCPD*(ptu(i, ikb)-ztenh(i, ikb)) + RLVTT*zqumqe
          zdh = RG * MAX(zdh, 1.0E5*zdqmin)
          IF (zdhpbl(i) > 0..AND.ldcum(i))mfub(i)=zdhpbl(i)/zdh
       ENDIF
       zmfmax = (paph(i, ikb)-paph(i, ikb-1)) / (RG*dtime)
       mfub(i) = MIN(mfub(i), zmfmax)
       zentr(i) = ENTRSCV
       IF (ktype(i) == 1) zentr(i) = ENTRPEN
    ENDDO

    ! DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME

    ! (A) calculer d'abord la hauteur "theorique" de la tour convective sans
    ! considerer l'entrainement ni le detrainement du panache, sachant
    ! ces derniers peuvent abaisser la hauteur theorique.

    DO i = 1, klon
       ikb=kcbot(i)
       zhcbase(i)=RCPD*ptu(i, ikb)+zgeoh(i, ikb)+RLVTT*pqu(i, ikb)
       ictop0(i)=kcbot(i)-1
    ENDDO

    zalvdcp=RLVTT/RCPD
    DO k=klev-1, 3, -1
       DO i = 1, klon
          zhsat=RCPD*ztenh(i, k)+zgeoh(i, k)+RLVTT*zqsenh(i, k)
          zgam=R5LES*zalvdcp*zqsenh(i, k)/ &
               ((1.-RETV *zqsenh(i, k))*(ztenh(i, k)-R4LES)**2)
          zzz=RCPD*ztenh(i, k)*0.608
          zhhat=zhsat-(zzz+zgam*zzz)/(1.+zgam*zzz/RLVTT)* &
               MAX(zqsenh(i, k)-zqenh(i, k), 0.)
          IF(k < ictop0(i).AND.zhcbase(i) > zhhat) ictop0(i)=k
       ENDDO
    ENDDO

    ! (B) calculer le panache ascendant

    CALL flxasc(dtime, ztenh, zqenh, ten, qen, qsen, pgeo, zgeoh, pap, &
         paph, pqte, pvervel, ldland, ldcum, ktype, ilab, ptu, pqu, plu, &
         mfu, mfub, zentr, mfus, mfuq, mful, plude, zdmfup, kcbot, &
         kctop, ictop0, kcum, pen_u, pde_u)

    kcum_not_zero: IF (kcum /= 0) then
       ! verifier l'epaisseur de la convection et changer eventuellement
       ! le taux d'entrainement/detrainement

       DO i = 1, klon
          zpbmpt=paph(i, kcbot(i))-paph(i, kctop(i))
          IF(ldcum(i) .AND. ktype(i) == 1 .AND. zpbmpt < 2E4) ktype(i) = 2
          IF(ldcum(i)) ictop0(i)=kctop(i)
          IF(ktype(i) == 2) zentr(i)=ENTRSCV
       ENDDO

       downdraft: IF (lmfdd) THEN
          ! si l'on considere le panache descendant
          ! calculer la precipitation issue du panache ascendant pour
          ! determiner l'existence du panache descendant dans la convection
          DO i = 1, klon
             zrfl(i)=zdmfup(i, 1)
          ENDDO
          DO k=2, klev
             DO i = 1, klon
                zrfl(i)=zrfl(i)+zdmfup(i, k)
             ENDDO
          ENDDO

          ! determiner le LFS (level of free sinking: niveau de plonge libre)
          CALL flxdlfs(ztenh, zqenh, zgeoh, paph, ptu, pqu, &
               ldcum, kcbot, kctop, mfub, zrfl, &
               ptd, pqd, &
               mfd, zmfds, zmfdq, zdmfdp, &
               kdtop, lddraf)

          ! calculer le panache descendant
          CALL flxddraf(ztenh, zqenh, &
               zgeoh, paph, zrfl, &
               ptd, pqd, &
               mfd, zmfds, zmfdq, zdmfdp, &
               lddraf, pen_d, pde_d)

          ! calculer de nouveau le flux de masse entrant a travers la base
          ! de la convection, sachant qu'il a ete modifie par le panache
          ! descendant
          DO i = 1, klon
             IF (lddraf(i)) THEN
                ikb = kcbot(i)
                llo1 = MFD(i, ikb) < 0.
                zeps = 0.
                IF (llo1) zeps = CMFDEPS
                zqumqe = pqu(i, ikb)+plu(i, ikb)- &
                     zeps*pqd(i, ikb)-(1.-zeps)*zqenh(i, ikb)
                zdqmin = MAX(0.01*zqenh(i, ikb), 1.E-10)
                zmfmax = (paph(i, ikb)-paph(i, ikb-1)) / (RG*dtime)
                IF (zdqpbl(i) > 0..AND.zqumqe > zdqmin.AND.ldcum(i) &
                     .AND.mfub(i) < zmfmax) THEN
                   mfub1(i) = zdqpbl(i) / (RG*MAX(zqumqe, zdqmin))
                ELSE
                   mfub1(i) = mfub(i)
                ENDIF
                IF (ktype(i) == 2) THEN
                   zdh = RCPD*(ptu(i, ikb)-zeps*ptd(i, ikb)- &
                        (1.-zeps)*ztenh(i, ikb))+RLVTT*zqumqe
                   zdh = RG * MAX(zdh, 1.0E5*zdqmin)
                   IF (zdhpbl(i) > 0..AND.ldcum(i))mfub1(i)=zdhpbl(i)/zdh
                ENDIF
                IF (.NOT. ((ktype(i) == 1 .OR. ktype(i) == 2) .AND. &
                     ABS(mfub1(i)-mfub(i)) < 0.2*mfub(i))) &
                     mfub1(i) = mfub(i)
             ENDIF
          ENDDO
          DO k = 1, klev
             DO i = 1, klon
                IF (lddraf(i)) THEN
                   zfac = mfub1(i)/MAX(mfub(i), 1.E-10)
                   mfd(i, k) = mfd(i, k)*zfac
                   zmfds(i, k) = zmfds(i, k)*zfac
                   zmfdq(i, k) = zmfdq(i, k)*zfac
                   zdmfdp(i, k) = zdmfdp(i, k)*zfac
                   pen_d(i, k) = pen_d(i, k)*zfac
                   pde_d(i, k) = pde_d(i, k)*zfac
                ENDIF
             ENDDO
          ENDDO
          DO i = 1, klon
             IF (lddraf(i)) mfub(i)=mfub1(i)
          ENDDO
       ENDIF downdraft

       ! calculer de nouveau le panache ascendant

       CALL flxasc(dtime, ztenh, zqenh, ten, qen, qsen, pgeo, zgeoh, pap, &
            paph, pqte, pvervel, ldland, ldcum, ktype, ilab, ptu, pqu, plu, &
            mfu, mfub, zentr, mfus, mfuq, mful, plude, zdmfup, kcbot, &
            kctop, ictop0, kcum, pen_u, pde_u)

       ! Déterminer les flux convectifs en forme finale, ainsi que la
       ! quantité des précipitations

       CALL flxflux(dtime, qen, qsen, ztenh, zqenh, pap, paph, &
            ldland, zgeoh, kcbot, kctop, lddraf, kdtop, ktype, ldcum, &
            mfu, mfd, mfus, zmfds, mfuq, zmfdq, mful, plude, &
            zdmfup, zdmfdp, ten, prsfc, pssfc, zdpmel, itopm2, &
            pmflxr, pmflxs)

       ! calculer les tendances pour T et Q

       CALL flxdtdq(itopm2, paph, ldcum, ten, mfus, zmfds, mfuq, zmfdq, &
            mful, zdmfup, zdmfdp, zdpmel, dt_con, dq_con)
    end IF kcum_not_zero

  END SUBROUTINE flxmain

end module flxmain_m
