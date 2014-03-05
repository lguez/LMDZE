
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/fyhyp.F,v 1.2 2005/06/03 09:11:32
! fairhead Exp $



SUBROUTINE fyhyp(yzoomdeg, grossism, dzooma, tau, rrlatu, yyprimu, rrlatv, &
    yyprimv, rlatu2, yprimu2, rlatu1, yprimu1, champmin, champmax)

  ! c    ...  Version du 01/04/2001 ....

  USE dimens_m
  USE paramet_m
  IMPLICIT NONE

  ! ...   Auteur :  P. Le Van  ...

  ! .......    d'apres  formulations  de R. Sadourny  .......

  ! Calcule les latitudes et derivees dans la grille du GCM pour une
  ! fonction f(y) a tangente  hyperbolique  .

  ! grossism etant le grossissement ( = 2 si 2 fois, = 3 si 3 fois , etc)
  ! dzoom  etant  la distance totale de la zone du zoom ( en radians )
  ! tau  la raideur de la transition de l'interieur a l'exterieur du zoom


  ! N.B : Il vaut mieux avoir : grossism * dzoom  <  pi/2  (radians) ,en
  ! lati.
  ! ********************************************************************



  INTEGER nmax, nmax2
  PARAMETER (nmax=30000, nmax2=2*nmax)


  ! .......  arguments  d'entree    .......

  REAL yzoomdeg, grossism, dzooma, tau
  ! ( rentres  par  run.def )

  ! .......  arguments  de sortie   .......

  REAL rrlatu(jjp1), yyprimu(jjp1), rrlatv(jjm), yyprimv(jjm), rlatu1(jjm), &
    yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)


  ! .....     champs  locaux    .....


  REAL dzoom
  DOUBLE PRECISION ylat(jjp1), yprim(jjp1)
  DOUBLE PRECISION yuv
  DOUBLE PRECISION yt(0:nmax2)
  DOUBLE PRECISION fhyp(0:nmax2), beta, ytprim(0:nmax2), fxm(0:nmax2)
  SAVE ytprim, yt, yf
  DOUBLE PRECISION yf(0:nmax2), yypr(0:nmax2)
  DOUBLE PRECISION yvrai(jjp1), yprimm(jjp1), ylatt(jjp1)
  DOUBLE PRECISION pi, depi, pis2, epsilon, y0, pisjm
  DOUBLE PRECISION yo1, yi, ylon2, ymoy, yprimin, champmin, champmax
  DOUBLE PRECISION yfi, yf1, ffdy
  DOUBLE PRECISION ypn, deply, y00
  SAVE y00, deply

  INTEGER i, j, it, ik, iter, jlat
  INTEGER jpn, jjpn
  SAVE jpn
  DOUBLE PRECISION a0, a1, a2, a3, yi2, heavyy0, heavyy0m
  DOUBLE PRECISION fa(0:nmax2), fb(0:nmax2)
  REAL y0min, y0max

  DOUBLE PRECISION heavyside

  pi = 2.*asin(1.)
  depi = 2.*pi
  pis2 = pi/2.
  pisjm = pi/float(jjm)
  epsilon = 1.E-3
  y0 = yzoomdeg*pi/180.

  IF (dzooma<1.) THEN
    dzoom = dzooma*pi
  ELSE IF (dzooma<12.) THEN
    WRITE (6, *) ' Le param. dzoomy pour fyhyp est trop petit &
      &! L aug                                                &
      &            menter et relancer'
    STOP 1
  ELSE
    dzoom = dzooma*pi/180.
  END IF

  WRITE (6, 18)
  WRITE (6, *) ' yzoom( rad.),grossism,tau,dzoom (radians)'
  WRITE (6, 24) y0, grossism, tau, dzoom

  DO i = 0, nmax2
    yt(i) = -pis2 + float(i)*pi/nmax2
  END DO

  heavyy0m = heavyside(-y0)
  heavyy0 = heavyside(y0)
  y0min = 2.*y0*heavyy0m - pis2
  y0max = 2.*y0*heavyy0 + pis2

  fa = 999.999
  fb = 999.999

  DO i = 0, nmax2
    IF (yt(i)<y0) THEN
      fa(i) = tau*(yt(i)-y0+dzoom/2.)
      fb(i) = (yt(i)-2.*y0*heavyy0m+pis2)*(y0-yt(i))
    ELSE IF (yt(i)>y0) THEN
      fa(i) = tau*(y0-yt(i)+dzoom/2.)
      fb(i) = (2.*y0*heavyy0-yt(i)+pis2)*(yt(i)-y0)
    END IF

    IF (200.*fb(i)<-fa(i)) THEN
      fhyp(i) = -1.
    ELSE IF (200.*fb(i)<fa(i)) THEN
      fhyp(i) = 1.
    ELSE
      fhyp(i) = tanh(fa(i)/fb(i))
    END IF

    IF (yt(i)==y0) fhyp(i) = 1.
    IF (yt(i)==y0min .OR. yt(i)==y0max) fhyp(i) = -1.

  END DO

  ! c  ....  Calcul  de  beta  ....

  ffdy = 0.

  DO i = 1, nmax2
    ymoy = 0.5*(yt(i-1)+yt(i))
    IF (ymoy<y0) THEN
      fa(i) = tau*(ymoy-y0+dzoom/2.)
      fb(i) = (ymoy-2.*y0*heavyy0m+pis2)*(y0-ymoy)
    ELSE IF (ymoy>y0) THEN
      fa(i) = tau*(y0-ymoy+dzoom/2.)
      fb(i) = (2.*y0*heavyy0-ymoy+pis2)*(ymoy-y0)
    END IF

    IF (200.*fb(i)<-fa(i)) THEN
      fxm(i) = -1.
    ELSE IF (200.*fb(i)<fa(i)) THEN
      fxm(i) = 1.
    ELSE
      fxm(i) = tanh(fa(i)/fb(i))
    END IF
    IF (ymoy==y0) fxm(i) = 1.
    IF (ymoy==y0min .OR. yt(i)==y0max) fxm(i) = -1.
    ffdy = ffdy + fxm(i)*(yt(i)-yt(i-1))

  END DO

  beta = (grossism*ffdy-pi)/(ffdy-pi)

  IF (2.*beta-grossism<=0.) THEN

    WRITE (6, *) ' **  Attention ! La valeur beta calculee dans &
      &la rou                                                 &
      &           tine fyhyp est mauvaise'
    WRITE (6, *) 'Modifier les valeurs de  grossismy ,tauy ou dzoomy', &
      ' et relancer ! ***  '
    STOP 1

  END IF

  ! .....  calcul  de  Ytprim   .....


  DO i = 0, nmax2
    ytprim(i) = beta + (grossism-beta)*fhyp(i)
  END DO

  ! .....  Calcul  de  Yf     ........

  yf(0) = -pis2
  DO i = 1, nmax2
    yypr(i) = beta + (grossism-beta)*fxm(i)
  END DO

  DO i = 1, nmax2
    yf(i) = yf(i-1) + yypr(i)*(yt(i)-yt(i-1))
  END DO

  ! ****************************************************************

  ! .....   yuv  = 0.   si calcul des latitudes  aux pts.  U  .....
  ! .....   yuv  = 0.5  si calcul des latitudes  aux pts.  V  .....

  WRITE (6, 18)

  DO ik = 1, 4

    IF (ik==1) THEN
      yuv = 0.
      jlat = jjm + 1
    ELSE IF (ik==2) THEN
      yuv = 0.5
      jlat = jjm
    ELSE IF (ik==3) THEN
      yuv = 0.25
      jlat = jjm
    ELSE IF (ik==4) THEN
      yuv = 0.75
      jlat = jjm
    END IF

    yo1 = 0.
    DO j = 1, jlat
      yo1 = 0.
      ylon2 = -pis2 + pisjm*(float(j)+yuv-1.)
      yfi = ylon2

      DO it = nmax2, 0, -1
        IF (yfi>=yf(it)) GO TO 350
      END DO
      it = 0
350   CONTINUE

      yi = yt(it)
      IF (it==nmax2) THEN
        it = nmax2 - 1
        yf(it+1) = pis2
      END IF
      ! .................................................................
      ! ....  Interpolation entre  yi(it) et yi(it+1)   pour avoir Y(yi)
      ! .....           et   Y'(yi)                             .....
      ! .................................................................

      CALL coefpoly(yf(it), yf(it+1), ytprim(it), ytprim(it+1), yt(it), &
        yt(it+1), a0, a1, a2, a3)

      yf1 = yf(it)
      yprimin = a1 + 2.*a2*yi + 3.*a3*yi*yi

      DO iter = 1, 300
        yi = yi - (yf1-yfi)/yprimin

        IF (abs(yi-yo1)<=epsilon) GO TO 550
        yo1 = yi
        yi2 = yi*yi
        yf1 = a0 + a1*yi + a2*yi2 + a3*yi2*yi
        yprimin = a1 + 2.*a2*yi + 3.*a3*yi2
      END DO
      WRITE (6, *) ' Pas de solution ***** ', j, ylon2, iter
      STOP 2
550   CONTINUE

      yprimin = a1 + 2.*a2*yi + 3.*a3*yi*yi
      yprim(j) = pi/(jjm*yprimin)
      yvrai(j) = yi

    END DO

    DO j = 1, jlat - 1
      IF (yvrai(j+1)<yvrai(j)) THEN
        WRITE (6, *) ' PBS. avec  rlat(', j + 1, ') plus petit que rlat(', j, &
          ')'
        STOP 3
      END IF
    END DO

    WRITE (6, *) 'Reorganisation des latitudes pour avoir entre - pi/2', &
      ' et  pi/2 '

    IF (ik==1) THEN
      ypn = pis2
      DO j = jlat, 1, -1
        IF (yvrai(j)<=ypn) GO TO 1502
      END DO
1502  CONTINUE

      jpn = j
      y00 = yvrai(jpn)
      deply = pis2 - y00
    END IF

    DO j = 1, jjm + 1 - jpn
      ylatt(j) = -pis2 - y00 + yvrai(jpn+j-1)
      yprimm(j) = yprim(jpn+j-1)
    END DO

    jjpn = jpn
    IF (jlat==jjm) jjpn = jpn - 1

    DO j = 1, jjpn
      ylatt(j+jjm+1-jpn) = yvrai(j) + deply
      yprimm(j+jjm+1-jpn) = yprim(j)
    END DO

    ! ***********   Fin de la reorganisation     *************


    DO j = 1, jlat
      ylat(j) = ylatt(jlat+1-j)
      yprim(j) = yprimm(jlat+1-j)
    END DO

    DO j = 1, jlat
      yvrai(j) = ylat(j)*180./pi
    END DO

    IF (ik==1) THEN
      ! WRITE(6,18)
      ! WRITE(6,*)  ' YLAT  en U   apres ( en  deg. ) '
      ! WRITE(6,68) (yvrai(j),j=1,jlat)
      ! c         WRITE(6,*) ' YPRIM '
      ! c         WRITE(6,445) ( yprim(j),j=1,jlat)

      DO j = 1, jlat
        rrlatu(j) = ylat(j)
        yyprimu(j) = yprim(j)
      END DO

    ELSE IF (ik==2) THEN
      ! WRITE(6,18)
      ! WRITE(6,*) ' YLAT   en V  apres ( en  deg. ) '
      ! WRITE(6,68) (yvrai(j),j=1,jlat)
      ! c         WRITE(6,*)' YPRIM '
      ! c         WRITE(6,445) ( yprim(j),j=1,jlat)

      DO j = 1, jlat
        rrlatv(j) = ylat(j)
        yyprimv(j) = yprim(j)
      END DO

    ELSE IF (ik==3) THEN
      ! WRITE(6,18)
      ! WRITE(6,*)  ' YLAT  en U + 0.75  apres ( en  deg. ) '
      ! WRITE(6,68) (yvrai(j),j=1,jlat)
      ! c         WRITE(6,*) ' YPRIM '
      ! c         WRITE(6,445) ( yprim(j),j=1,jlat)

      DO j = 1, jlat
        rlatu2(j) = ylat(j)
        yprimu2(j) = yprim(j)
      END DO

    ELSE IF (ik==4) THEN
      ! WRITE(6,18)
      ! WRITE(6,*)  ' YLAT en U + 0.25  apres ( en  deg. ) '
      ! WRITE(6,68)(yvrai(j),j=1,jlat)
      ! c         WRITE(6,*) ' YPRIM '
      ! c         WRITE(6,68) ( yprim(j),j=1,jlat)

      DO j = 1, jlat
        rlatu1(j) = ylat(j)
        yprimu1(j) = yprim(j)
      END DO

    END IF

  END DO

  WRITE (6, 18)

  ! .....     fin de la boucle  do 5000 .....

  DO j = 1, jjm
    ylat(j) = rrlatu(j) - rrlatu(j+1)
  END DO
  champmin = 1.E12
  champmax = -1.E12
  DO j = 1, jjm
    champmin = min(champmin, ylat(j))
    champmax = max(champmax, ylat(j))
  END DO
  champmin = champmin*180./pi
  champmax = champmax*180./pi

24 FORMAT (2X, 'Parametres yzoom,gross,tau ,dzoom pour fyhyp ', 4F8.3)
18 FORMAT (/)
68 FORMAT (1X, 7F9.2)

  RETURN
END SUBROUTINE fyhyp
