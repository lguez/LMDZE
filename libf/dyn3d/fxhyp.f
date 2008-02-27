!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/fxhyp.F,v 1.2 2005/06/03 09:11:32 fairhead Exp $
!
c
c
       SUBROUTINE fxhyp ( xzoomdeg,grossism,dzooma,tau ,
     , rlonm025,xprimm025,rlonv,xprimv,rlonu,xprimu,rlonp025,xprimp025,
     , champmin,champmax                                               )

c      Auteur :  P. Le Van 

       use dimens_m
      use paramet_m
       IMPLICIT NONE

c    Calcule les longitudes et derivees dans la grille du GCM pour une
c     fonction f(x) a tangente  hyperbolique  .
c
c     grossism etant le grossissement ( = 2 si 2 fois, = 3 si 3 fois,etc.)
c     dzoom  etant  la distance totale de la zone du zoom
c     tau  la raideur de la transition de l'interieur a l'exterieur du zoom
c
c    On doit avoir grossism x dzoom <  pi ( radians )   , en longitude.
c   ********************************************************************


       INTEGER nmax, nmax2
       PARAMETER (  nmax = 30000, nmax2 = 2*nmax )
c
       LOGICAL scal180
       PARAMETER ( scal180 = .TRUE. )

c      scal180 = .TRUE.  si on veut avoir le premier point scalaire pour   
c      une grille reguliere ( grossism = 1.,tau=0.,clon=0. ) a -180. degres.
c      sinon scal180 = .FALSE.

       
c     ......  arguments  d'entree   .......
c
       REAL xzoomdeg,dzooma,tau,grossism

c    ......   arguments  de  sortie  ......

       REAL rlonm025(iip1),xprimm025(iip1),rlonv(iip1),xprimv(iip1),
     ,  rlonu(iip1),xprimu(iip1),rlonp025(iip1),xprimp025(iip1)

c     .... variables locales  ....
c
       REAL   dzoom
       REAL*8 xlon(iip1),xprimm(iip1),xuv
       REAL*8 xtild(0:nmax2)
       REAL*8 fhyp(0:nmax2),ffdx,beta,Xprimt(0:nmax2)
       REAL*8 Xf(0:nmax2),xxpr(0:nmax2)
       REAL*8 xvrai(iip1),xxprim(iip1) 
       REAL*8 pi,depi,epsilon,xzoom,fa,fb
       REAL*8 Xf1, Xfi , a0,a1,a2,a3,xi2
       INTEGER i,it,ik,iter,ii,idif,ii1,ii2
       REAL*8 xi,xo1,xmoy,xlon2,fxm,Xprimin
       REAL*8 champmin,champmax,decalx
       INTEGER is2
       SAVE is2

       REAL*8 heavyside

       pi       = 2. * ASIN(1.)
       depi     = 2. * pi
       epsilon  = 1.e-3
       xzoom    = xzoomdeg * pi/180. 
c
           decalx   = .75
       IF( grossism.EQ.1..AND.scal180 )  THEN
           decalx   = 1.
       ENDIF

       WRITE(6,*) 'FXHYP scal180,decalx', scal180,decalx
c
       IF( dzooma.LT.1.)  THEN
         dzoom = dzooma * depi
       ELSEIF( dzooma.LT. 25. ) THEN
         WRITE(6,*) ' Le param. dzoomx pour fxhyp est trop petit ! L aug
     ,menter et relancer ! '
         STOP 1
       ELSE
         dzoom = dzooma * pi/180.
       ENDIF

       WRITE(6,*) ' xzoom( rad.),grossism,tau,dzoom (radians)'
       WRITE(6,24) xzoom,grossism,tau,dzoom

       DO i = 0, nmax2 
        xtild(i) = - pi + FLOAT(i) * depi /nmax2
       ENDDO

       DO i = nmax, nmax2

       fa  = tau*  ( dzoom/2.  - xtild(i) )
       fb  = xtild(i) *  ( pi - xtild(i) )

         IF( 200.* fb .LT. - fa )   THEN
           fhyp ( i) = - 1.
         ELSEIF( 200. * fb .LT. fa ) THEN
           fhyp ( i) =   1.
         ELSE
            IF( ABS(fa).LT.1.e-13.AND.ABS(fb).LT.1.e-13)  THEN
                IF(   200.*fb + fa.LT.1.e-10 )  THEN
                    fhyp ( i ) = - 1.
                ELSEIF( 200.*fb - fa.LT.1.e-10 )  THEN
                    fhyp ( i )  =   1.
                ENDIF
            ELSE
                    fhyp ( i )  =  TANH ( fa/fb )
            ENDIF
         ENDIF
        IF ( xtild(i).EQ. 0. )  fhyp(i) =  1.
        IF ( xtild(i).EQ. pi )  fhyp(i) = -1.

       ENDDO

cc  ....  Calcul  de  beta  ....

       ffdx = 0.

       DO i = nmax +1,nmax2

       xmoy    = 0.5 * ( xtild(i-1) + xtild( i ) )
       fa  = tau*  ( dzoom/2.  - xmoy )
       fb  = xmoy *  ( pi - xmoy )

       IF( 200.* fb .LT. - fa )   THEN
         fxm = - 1.
       ELSEIF( 200. * fb .LT. fa ) THEN
         fxm =   1.
       ELSE
            IF( ABS(fa).LT.1.e-13.AND.ABS(fb).LT.1.e-13)  THEN
                IF(   200.*fb + fa.LT.1.e-10 )  THEN
                    fxm   = - 1.
                ELSEIF( 200.*fb - fa.LT.1.e-10 )  THEN
                    fxm   =   1.
                ENDIF
            ELSE
                    fxm   =  TANH ( fa/fb )
            ENDIF
       ENDIF

       IF ( xmoy.EQ. 0. )  fxm  =  1.
       IF ( xmoy.EQ. pi )  fxm  = -1.

       ffdx = ffdx + fxm * ( xtild(i) - xtild(i-1) )

       ENDDO

        beta  = ( grossism * ffdx - pi ) / ( ffdx - pi )

       IF( 2.*beta - grossism.LE. 0.)  THEN
        WRITE(6,*) ' **  Attention ! La valeur beta calculee dans la rou
     ,tine fxhyp est mauvaise ! '
        WRITE(6,*)'Modifier les valeurs de  grossismx ,tau ou dzoomx ',
     , ' et relancer ! ***  '
        STOP 1
       ENDIF
c
c   .....  calcul  de  Xprimt   .....
c
       
       DO i = nmax, nmax2
        Xprimt(i) = beta  + ( grossism - beta ) * fhyp(i)
       ENDDO
c   
       DO i =  nmax+1, nmax2
        Xprimt( nmax2 - i ) = Xprimt( i )
       ENDDO
c

c   .....  Calcul  de  Xf     ........

       Xf(0) = - pi

       DO i =  nmax +1, nmax2

       xmoy    = 0.5 * ( xtild(i-1) + xtild( i ) )
       fa  = tau*  ( dzoom/2.  - xmoy )
       fb  = xmoy *  ( pi - xmoy )

       IF( 200.* fb .LT. - fa )   THEN
         fxm = - 1.
       ELSEIF( 200. * fb .LT. fa ) THEN
         fxm =   1.
       ELSE
         fxm =  TANH ( fa/fb )
       ENDIF

       IF ( xmoy.EQ. 0. )  fxm =  1.
       IF ( xmoy.EQ. pi )  fxm = -1.
       xxpr(i)    = beta + ( grossism - beta ) * fxm

       ENDDO

       DO i = nmax+1, nmax2
        xxpr(nmax2-i+1) = xxpr(i)
       ENDDO

        DO i=1,nmax2
         Xf(i)   = Xf(i-1) + xxpr(i) * ( xtild(i) - xtild(i-1) )
        ENDDO


c    *****************************************************************
c

c     .....  xuv = 0.   si  calcul  aux pts   scalaires   ........
c     .....  xuv = 0.5  si  calcul  aux pts      U        ........
c
      WRITE(6,18)
c
      DO 5000  ik = 1, 4

       IF( ik.EQ.1 )        THEN
         xuv =  -0.25
       ELSE IF ( ik.EQ.2 )  THEN
         xuv =   0.
       ELSE IF ( ik.EQ.3 )  THEN
         xuv =   0.50
       ELSE IF ( ik.EQ.4 )  THEN
         xuv =   0.25
       ENDIF

      xo1   = 0.

      ii1=1
      ii2=iim
      IF(ik.EQ.1.and.grossism.EQ.1.) THEN
        ii1 = 2 
        ii2 = iim+1
      ENDIF
      DO 1500 i = ii1, ii2

      xlon2 = - pi + (FLOAT(i) + xuv - decalx) * depi / FLOAT(iim) 

      Xfi    = xlon2
c
      DO 250 it =  nmax2,0,-1
      IF( Xfi.GE.Xf(it))  GO TO 350
250   CONTINUE

      it = 0

350   CONTINUE

c    ......  Calcul de   Xf(xi)    ...... 
c
      xi  = xtild(it)

      IF(it.EQ.nmax2)  THEN
       it       = nmax2 -1
       Xf(it+1) = pi
      ENDIF
c  .....................................................................
c
c   Appel de la routine qui calcule les coefficients a0,a1,a2,a3 d'un
c   polynome de degre 3  qui passe  par les points (Xf(it),xtild(it) )
c          et (Xf(it+1),xtild(it+1) )

       CALL coefpoly ( Xf(it),Xf(it+1),Xprimt(it),Xprimt(it+1),
     ,                xtild(it),xtild(it+1),  a0, a1, a2, a3  )

       Xf1     = Xf(it)
       Xprimin = a1 + 2.* a2 * xi + 3.*a3 * xi *xi

       DO 500 iter = 1,300
        xi = xi - ( Xf1 - Xfi )/ Xprimin

        IF( ABS(xi-xo1).LE.epsilon)  GO TO 550
         xo1      = xi
         xi2      = xi * xi
         Xf1      = a0 +  a1 * xi +     a2 * xi2  +     a3 * xi2 * xi
         Xprimin  =       a1      + 2.* a2 *  xi  + 3.* a3 * xi2
500   CONTINUE
        WRITE(6,*) ' Pas de solution ***** ',i,xlon2,iter
          STOP 6
550   CONTINUE

       xxprim(i) = depi/ ( FLOAT(iim) * Xprimin )
       xvrai(i)  =  xi + xzoom

1500   CONTINUE


       IF(ik.EQ.1.and.grossism.EQ.1.)  THEN
         xvrai(1)    = xvrai(iip1)-depi
         xxprim(1)   = xxprim(iip1)
       ENDIF
       DO i = 1 , iim
        xlon(i)     = xvrai(i)
        xprimm(i)   = xxprim(i)
       ENDDO
       DO i = 1, iim -1
        IF( xvrai(i+1). LT. xvrai(i) )  THEN
         WRITE(6,*) ' PBS. avec rlonu(',i+1,') plus petit que rlonu(',i,
     ,  ')'
        STOP 7
        ENDIF
       ENDDO
c
c   ... Reorganisation  des  longitudes  pour les avoir  entre - pi et pi ..
c   ........................................................................

       champmin =  1.e12
       champmax = -1.e12
       DO i = 1, iim
        champmin = MIN( champmin,xvrai(i) )
        champmax = MAX( champmax,xvrai(i) )
       ENDDO

      IF(champmin .GE.-pi-0.10.and.champmax.LE.pi+0.10 )  THEN
                GO TO 1600
      ELSE
       WRITE(6,*) 'Reorganisation des longitudes pour avoir entre - pi',
     ,  ' et pi '
c
        IF( xzoom.LE.0.)  THEN
          IF( ik.EQ. 1 )  THEN
          DO i = 1, iim
           IF( xvrai(i).GE. - pi )  GO TO 80
          ENDDO
            WRITE(6,*)  ' PBS. 1 !  Xvrai plus petit que  - pi ! '
            STOP 8
 80       CONTINUE
          is2 = i
          ENDIF

          IF( is2.NE. 1 )  THEN
            DO ii = is2 , iim
             xlon  (ii-is2+1) = xvrai(ii)
             xprimm(ii-is2+1) = xxprim(ii)
            ENDDO
            DO ii = 1 , is2 -1
             xlon  (ii+iim-is2+1) = xvrai(ii) + depi
             xprimm(ii+iim-is2+1) = xxprim(ii) 
            ENDDO
          ENDIF
        ELSE 
          IF( ik.EQ.1 )  THEN
           DO i = iim,1,-1
             IF( xvrai(i).LE. pi ) GO TO 90
           ENDDO
             WRITE(6,*) ' PBS.  2 ! Xvrai plus grand  que   pi ! '
              STOP 9
 90        CONTINUE
            is2 = i
          ENDIF
           idif = iim -is2
           DO ii = 1, is2
            xlon  (ii+idif) = xvrai(ii)
            xprimm(ii+idif) = xxprim(ii)
           ENDDO
           DO ii = 1, idif
            xlon (ii)  = xvrai (ii+is2) - depi
            xprimm(ii) = xxprim(ii+is2) 
           ENDDO
         ENDIF
      ENDIF
c
c     .........   Fin  de la reorganisation   ............................

 1600    CONTINUE


         xlon  ( iip1)  = xlon(1) + depi
         xprimm( iip1 ) = xprimm (1 )
       
         DO i = 1, iim+1
         xvrai(i) = xlon(i)*180./pi
         ENDDO

         IF( ik.EQ.1 )  THEN
c          WRITE(6,*)  ' XLON aux pts. V-0.25   apres ( en  deg. ) '
c          WRITE(6,18) 
c          WRITE(6,68) xvrai
c          WRITE(6,*) ' XPRIM k ',ik
c          WRITE(6,566)  xprimm

           DO i = 1,iim +1
             rlonm025(i) = xlon( i )
            xprimm025(i) = xprimm(i)
           ENDDO
         ELSE IF( ik.EQ.2 )  THEN
c          WRITE(6,18) 
c          WRITE(6,*)  ' XLON aux pts. V   apres ( en  deg. ) '
c          WRITE(6,68) xvrai
c          WRITE(6,*) ' XPRIM k ',ik
c          WRITE(6,566)  xprimm

           DO i = 1,iim + 1
             rlonv(i) = xlon( i )
            xprimv(i) = xprimm(i)
           ENDDO

         ELSE IF( ik.EQ.3)   THEN
c          WRITE(6,18) 
c          WRITE(6,*)  ' XLON aux pts. U   apres ( en  deg. ) '
c          WRITE(6,68) xvrai
c          WRITE(6,*) ' XPRIM ik ',ik
c          WRITE(6,566)  xprimm

           DO i = 1,iim + 1
             rlonu(i) = xlon( i )
            xprimu(i) = xprimm(i)
           ENDDO

         ELSE IF( ik.EQ.4 )  THEN
c          WRITE(6,18) 
c          WRITE(6,*)  ' XLON aux pts. V+0.25   apres ( en  deg. ) '
c          WRITE(6,68) xvrai
c          WRITE(6,*) ' XPRIM ik ',ik
c          WRITE(6,566)  xprimm

           DO i = 1,iim + 1
             rlonp025(i) = xlon( i )
            xprimp025(i) = xprimm(i)
           ENDDO

         ENDIF

5000    CONTINUE
c
       WRITE(6,18)
c
c    ...........  fin  de la boucle  do 5000      ............

        DO i = 1, iim
         xlon(i) = rlonv(i+1) - rlonv(i)
        ENDDO
        champmin =  1.e12
        champmax = -1.e12
        DO i = 1, iim
         champmin = MIN( champmin, xlon(i) )
         champmax = MAX( champmax, xlon(i) )
        ENDDO
         champmin = champmin * 180./pi
         champmax = champmax * 180./pi

18     FORMAT(/)
24     FORMAT(2x,'Parametres xzoom,gross,tau ,dzoom pour fxhyp ',4f8.3)
68     FORMAT(1x,7f9.2)
566    FORMAT(1x,7f9.4)

       RETURN
       END
