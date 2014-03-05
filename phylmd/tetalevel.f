! ================================================================
! ================================================================
SUBROUTINE tetalevel(ilon, ilev, lnew, pgcm, pres, qgcm, qpres)
  ! ================================================================
  ! ================================================================

  USE dimens_m
  USE paramet_m
  USE dimphy
  IMPLICIT NONE


  ! ================================================================

  ! Interpoler des champs 3-D u, v et g du modele a un niveau de
  ! pression donnee (pres)

  ! INPUT:  ilon ----- nombre de points
  ! ilev ----- nombre de couches
  ! lnew ----- true si on doit reinitialiser les poids
  ! pgcm ----- pressions modeles
  ! pres ----- pression vers laquelle on interpolle
  ! Qgcm ----- champ GCM
  ! Qpres ---- champ interpolle au niveau pres

  ! ================================================================

  ! arguments :
  ! -----------

  INTEGER ilon, ilev
  LOGICAL lnew

  REAL pgcm(ilon, ilev)
  REAL, INTENT (IN) :: qgcm(ilon, ilev)
  REAL pres
  REAL qpres(ilon)

  ! local :
  ! -------


  INTEGER lt(ip1jmp1), lb(ip1jmp1)
  REAL ptop, pbot, aist(ip1jmp1), aisb(ip1jmp1)
  SAVE lt, lb, ptop, pbot, aist, aisb

  INTEGER i, k

  ! PRINT*,'tetalevel pres=',pres
  ! =====================================================================
  IF (lnew) THEN
    ! on réinitialise les réindicages et les poids
    ! =====================================================================


    ! Chercher les 2 couches les plus proches du niveau a obtenir

    ! Eventuellement, faire l'extrapolation a partir des deux couches
    ! les plus basses ou les deux couches les plus hautes:
    DO i = 1, ilon
      ! IM      IF ( ABS(pres-pgcm(i,ilev) ) .LT.
      IF (abs(pres-pgcm(i,ilev))>abs(pres-pgcm(i,1))) THEN
        lt(i) = ilev ! 2
        lb(i) = ilev - 1 ! 1
      ELSE
        lt(i) = 2
        lb(i) = 1
      END IF
      ! IM   PRINT*,'i, ABS(pres-pgcm),ABS(pres-pgcm)',
      ! IM  .i, ABS(pres-pgcm(i,ilev)),ABS(pres-pgcm(i,1))
    END DO
    DO k = 1, ilev - 1
      DO i = 1, ilon
        pbot = pgcm(i, k)
        ptop = pgcm(i, k+1)
        ! IM         IF (ptop.LE.pres .AND. pbot.GE.pres) THEN
        IF (ptop>=pres .AND. pbot<=pres) THEN
          lt(i) = k + 1
          lb(i) = k
        END IF
      END DO
    END DO

    ! Interpolation lineaire:

    DO i = 1, ilon
      ! interpolation en logarithme de pression:

      ! ...   Modif . P. Le Van    ( 20/01/98) ....
      ! Modif Frédéric Hourdin (3/01/02)

      ! IF(pgcm(i,lb(i)).NE.0.OR.
      ! $     pgcm(i,lt(i)).NE.0.) THEN

      ! PRINT*,'i,lb,lt,2pgcm,pres',i,lb(i),
      ! .  lt(i),pgcm(i,lb(i)),pgcm(i,lt(i)),pres

      aist(i) = log(pgcm(i,lb(i))/pres)/log(pgcm(i,lb(i))/pgcm(i,lt(i)))
      aisb(i) = log(pres/pgcm(i,lt(i)))/log(pgcm(i,lb(i))/pgcm(i,lt(i)))
    END DO


  END IF ! lnew

  ! ======================================================================
  ! inteprollation
  ! ======================================================================

  DO i = 1, ilon
    qpres(i) = qgcm(i, lb(i))*aisb(i) + qgcm(i, lt(i))*aist(i)
  END DO

  ! Je mets les vents a zero quand je rencontre une montagne
  DO i = 1, ilon
    ! IM      if (pgcm(i,1).LT.pres) THEN
    IF (pgcm(i,1)>pres) THEN
      ! Qpres(i)=1e33
      qpres(i) = 1E+20
      ! IM         PRINT*,'i,pgcm(i,1),pres =',i,pgcm(i,1),pres
    END IF
  END DO


  RETURN
END SUBROUTINE tetalevel
