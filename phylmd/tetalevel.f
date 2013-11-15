c================================================================
c================================================================
      SUBROUTINE tetalevel(ilon,ilev,lnew,pgcm,pres,Qgcm,Qpres)
c================================================================
c================================================================

      use dimens_m
      use paramet_m
      use dimphy
      IMPLICIT none


c================================================================
c
c Interpoler des champs 3-D u, v et g du modele a un niveau de
c pression donnee (pres)
c
c INPUT:  ilon ----- nombre de points
c         ilev ----- nombre de couches
c         lnew ----- true si on doit reinitialiser les poids
c         pgcm ----- pressions modeles
c         pres ----- pression vers laquelle on interpolle
c         Qgcm ----- champ GCM
c         Qpres ---- champ interpolle au niveau pres
c
c================================================================
c
c   arguments :
c   -----------

      INTEGER ilon, ilev
      logical lnew

      REAL pgcm(ilon,ilev)
      REAL, intent(in):: Qgcm(ilon,ilev)
      real pres
      REAL Qpres(ilon)

c   local :
c   -------
c
c
      INTEGER lt(ip1jmp1), lb(ip1jmp1)
      REAL ptop, pbot, aist(ip1jmp1), aisb(ip1jmp1)
      save lt,lb,ptop,pbot,aist,aisb

      INTEGER i, k
c
c     PRINT*,'tetalevel pres=',pres
c=====================================================================
      if (lnew) then
c   on réinitialise les réindicages et les poids
c=====================================================================


c Chercher les 2 couches les plus proches du niveau a obtenir
c
c Eventuellement, faire l'extrapolation a partir des deux couches
c les plus basses ou les deux couches les plus hautes:
      DO 130 i = 1, ilon
cIM      IF ( ABS(pres-pgcm(i,ilev) ) .LT.
         IF ( ABS(pres-pgcm(i,ilev) ) .GT.
     .        ABS(pres-pgcm(i,1)) ) THEN
            lt(i) = ilev     ! 2
            lb(i) = ilev-1   ! 1
         ELSE
            lt(i) = 2
            lb(i) = 1
         ENDIF
cIM   PRINT*,'i, ABS(pres-pgcm),ABS(pres-pgcm)',
cIM  .i, ABS(pres-pgcm(i,ilev)),ABS(pres-pgcm(i,1))
  130 CONTINUE
      DO 150 k = 1, ilev-1
         DO 140 i = 1, ilon
            pbot = pgcm(i,k)
            ptop = pgcm(i,k+1)
cIM         IF (ptop.LE.pres .AND. pbot.GE.pres) THEN
            IF (ptop.GE.pres .AND. pbot.LE.pres) THEN
               lt(i) = k+1
               lb(i) = k
            ENDIF
  140    CONTINUE
  150 CONTINUE
c
c Interpolation lineaire:
c
      DO i = 1, ilon
c interpolation en logarithme de pression:
c
c ...   Modif . P. Le Van    ( 20/01/98) ....
c       Modif Frédéric Hourdin (3/01/02)

c       IF(pgcm(i,lb(i)).NE.0.OR.
c    $     pgcm(i,lt(i)).NE.0.) THEN
c
c       PRINT*,'i,lb,lt,2pgcm,pres',i,lb(i),
c    .  lt(i),pgcm(i,lb(i)),pgcm(i,lt(i)),pres
c
        aist(i) = LOG( pgcm(i,lb(i))/ pres )
     .       / LOG( pgcm(i,lb(i))/ pgcm(i,lt(i)) )
        aisb(i) = LOG( pres / pgcm(i,lt(i)) )
     .       / LOG( pgcm(i,lb(i))/ pgcm(i,lt(i)))
      enddo


      endif ! lnew

c======================================================================
c    inteprollation
c======================================================================

      do i=1,ilon
         Qpres(i)= Qgcm(i,lb(i))*aisb(i)+Qgcm(i,lt(i))*aist(i)
      enddo
c
c Je mets les vents a zero quand je rencontre une montagne
      do i = 1, ilon
cIM      if (pgcm(i,1).LT.pres) THEN
         if (pgcm(i,1).GT.pres) THEN
c           Qpres(i)=1e33
            Qpres(i)=1e+20
cIM         PRINT*,'i,pgcm(i,1),pres =',i,pgcm(i,1),pres
         endif
      enddo

c
      RETURN
      END
