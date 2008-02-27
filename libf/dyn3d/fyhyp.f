!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/fyhyp.F,v 1.2 2005/06/03 09:11:32 fairhead Exp $
!
c
c
       SUBROUTINE fyhyp ( yzoomdeg, grossism, dzooma,tau  ,  
     ,  rrlatu,yyprimu,rrlatv,yyprimv,rlatu2,yprimu2,rlatu1,yprimu1 ,
     ,  champmin,champmax                                            ) 

cc    ...  Version du 01/04/2001 ....

       use dimens_m
      use paramet_m
       IMPLICIT NONE
c
c    ...   Auteur :  P. Le Van  ... 
c
c    .......    d'apres  formulations  de R. Sadourny  .......
c
c     Calcule les latitudes et derivees dans la grille du GCM pour une
c     fonction f(y) a tangente  hyperbolique  .
c
c     grossism etant le grossissement ( = 2 si 2 fois, = 3 si 3 fois , etc)
c     dzoom  etant  la distance totale de la zone du zoom ( en radians )
c     tau  la raideur de la transition de l'interieur a l'exterieur du zoom   
c
c
c N.B : Il vaut mieux avoir : grossism * dzoom  <  pi/2  (radians) ,en lati.
c      ********************************************************************
c
c

       INTEGER      nmax , nmax2
       PARAMETER (  nmax = 30000, nmax2 = 2*nmax )
c
c
c     .......  arguments  d'entree    .......
c
       REAL yzoomdeg, grossism,dzooma,tau 
c         ( rentres  par  run.def )

c     .......  arguments  de sortie   .......
c
       REAL rrlatu(jjp1), yyprimu(jjp1),rrlatv(jjm), yyprimv(jjm),
     , rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)

c
c     .....     champs  locaux    .....
c
     
       REAL   dzoom
       REAL*8 ylat(jjp1), yprim(jjp1)
       REAL*8 yuv
       REAL*8 yt(0:nmax2)
       REAL*8 fhyp(0:nmax2),beta,Ytprim(0:nmax2),fxm(0:nmax2)
       SAVE Ytprim, yt,Yf
       REAL*8 Yf(0:nmax2),yypr(0:nmax2)
       REAL*8 yvrai(jjp1), yprimm(jjp1),ylatt(jjp1)
       REAL*8 pi,depi,pis2,epsilon,y0,pisjm
       REAL*8 yo1,yi,ylon2,ymoy,Yprimin,champmin,champmax
       REAL*8 yfi,Yf1,ffdy
       REAL*8 ypn,deply,y00
       SAVE y00, deply

       INTEGER i,j,it,ik,iter,jlat
       INTEGER jpn,jjpn
       SAVE jpn
       REAL*8 a0,a1,a2,a3,yi2,heavyy0,heavyy0m
       REAL*8 fa(0:nmax2),fb(0:nmax2)
       REAL y0min,y0max

       REAL*8     heavyside

       pi       = 2. * ASIN(1.)
       depi     = 2. * pi
       pis2     = pi/2.
       pisjm    = pi/ FLOAT(jjm)
       epsilon  = 1.e-3
       y0       =  yzoomdeg * pi/180. 

       IF( dzooma.LT.1.)  THEN
         dzoom = dzooma * pi
       ELSEIF( dzooma.LT. 12. ) THEN
         WRITE(6,*) ' Le param. dzoomy pour fyhyp est trop petit ! L aug
     ,menter et relancer ! '
         STOP 1
       ELSE
         dzoom = dzooma * pi/180.
       ENDIF

       WRITE(6,18)
       WRITE(6,*) ' yzoom( rad.),grossism,tau,dzoom (radians)'
       WRITE(6,24) y0,grossism,tau,dzoom

       DO i = 0, nmax2 
        yt(i) = - pis2  + FLOAT(i)* pi /nmax2
       ENDDO

       heavyy0m = heavyside( -y0 )
       heavyy0  = heavyside(  y0 )
       y0min    = 2.*y0*heavyy0m - pis2
       y0max    = 2.*y0*heavyy0  + pis2

       fa = 999.999
       fb = 999.999
       
       DO i = 0, nmax2 
        IF( yt(i).LT.y0 )  THEN
         fa (i) = tau*  (yt(i)-y0+dzoom/2. )
         fb(i) =   (yt(i)-2.*y0*heavyy0m +pis2) * ( y0 - yt(i) )
        ELSEIF ( yt(i).GT.y0 )  THEN
         fa(i) =   tau *(y0-yt(i)+dzoom/2. )
         fb(i) = (2.*y0*heavyy0 -yt(i)+pis2) * ( yt(i) - y0 ) 
       ENDIF
        
       IF( 200.* fb(i) .LT. - fa(i) )   THEN
         fhyp ( i) = - 1.
       ELSEIF( 200. * fb(i) .LT. fa(i) ) THEN
         fhyp ( i) =   1.
       ELSE  
         fhyp(i) =  TANH ( fa(i)/fb(i) )
       ENDIF

       IF( yt(i).EQ.y0 )  fhyp(i) = 1.
       IF(yt(i).EQ. y0min. OR.yt(i).EQ. y0max ) fhyp(i) = -1.

       ENDDO

cc  ....  Calcul  de  beta  ....
c
       ffdy   = 0.

       DO i = 1, nmax2
        ymoy    = 0.5 * ( yt(i-1) + yt( i ) )
        IF( ymoy.LT.y0 )  THEN
         fa(i)= tau * ( ymoy-y0+dzoom/2.) 
         fb(i) = (ymoy-2.*y0*heavyy0m +pis2) * ( y0 - ymoy )
        ELSEIF ( ymoy.GT.y0 )  THEN
         fa(i)= tau * ( y0-ymoy+dzoom/2. ) 
         fb(i) = (2.*y0*heavyy0 -ymoy+pis2) * ( ymoy - y0 )
        ENDIF

        IF( 200.* fb(i) .LT. - fa(i) )    THEN
         fxm ( i) = - 1.
        ELSEIF( 200. * fb(i) .LT. fa(i) ) THEN
         fxm ( i) =   1.
        ELSE
         fxm(i) =  TANH ( fa(i)/fb(i) )
        ENDIF
         IF( ymoy.EQ.y0 )  fxm(i) = 1.
         IF (ymoy.EQ. y0min. OR.yt(i).EQ. y0max ) fxm(i) = -1.
         ffdy = ffdy + fxm(i) * ( yt(i) - yt(i-1) )

        ENDDO

        beta  = ( grossism * ffdy - pi ) / ( ffdy - pi )

       IF( 2.*beta - grossism.LE. 0.)  THEN

        WRITE(6,*) ' **  Attention ! La valeur beta calculee dans la rou
     ,tine fyhyp est mauvaise ! '
        WRITE(6,*)'Modifier les valeurs de  grossismy ,tauy ou dzoomy',
     , ' et relancer ! ***  '
        STOP 1

       ENDIF
c
c   .....  calcul  de  Ytprim   .....
c
       
       DO i = 0, nmax2
        Ytprim(i) = beta  + ( grossism - beta ) * fhyp(i)
       ENDDO

c   .....  Calcul  de  Yf     ........

       Yf(0) = - pis2
       DO i = 1, nmax2
        yypr(i)    = beta + ( grossism - beta ) * fxm(i)
       ENDDO

       DO i=1,nmax2
        Yf(i)   = Yf(i-1) + yypr(i) * ( yt(i) - yt(i-1) )
       ENDDO

c    ****************************************************************
c
c   .....   yuv  = 0.   si calcul des latitudes  aux pts.  U  .....
c   .....   yuv  = 0.5  si calcul des latitudes  aux pts.  V  .....
c
      WRITE(6,18)
c
      DO 5000  ik = 1,4

       IF( ik.EQ.1 )  THEN
         yuv  = 0.
         jlat = jjm + 1
       ELSE IF ( ik.EQ.2 )  THEN
         yuv  = 0.5
         jlat = jjm 
       ELSE IF ( ik.EQ.3 )  THEN
         yuv  = 0.25
         jlat = jjm 
       ELSE IF ( ik.EQ.4 )  THEN
         yuv  = 0.75
         jlat = jjm 
       ENDIF
c
       yo1   = 0.
       DO 1500 j =  1,jlat
        yo1   = 0.
        ylon2 =  - pis2 + pisjm * ( FLOAT(j)  + yuv  -1.)  
        yfi    = ylon2
c
       DO 250 it =  nmax2,0,-1
        IF( yfi.GE.Yf(it))  GO TO 350
250    CONTINUE
       it = 0
350    CONTINUE

       yi = yt(it)
       IF(it.EQ.nmax2)  THEN
        it       = nmax2 -1
        Yf(it+1) = pis2
       ENDIF
c  .................................................................
c  ....  Interpolation entre  yi(it) et yi(it+1)   pour avoir Y(yi)  
c      .....           et   Y'(yi)                             .....
c  .................................................................

       CALL coefpoly ( Yf(it),Yf(it+1),Ytprim(it), Ytprim(it+1),   
     ,                  yt(it),yt(it+1) ,   a0,a1,a2,a3   )      

       Yf1     = Yf(it)
       Yprimin = a1 + 2.* a2 * yi + 3.*a3 * yi *yi

       DO 500 iter = 1,300
         yi = yi - ( Yf1 - yfi )/ Yprimin

        IF( ABS(yi-yo1).LE.epsilon)  GO TO 550
         yo1      = yi
         yi2      = yi * yi
         Yf1      = a0 +  a1 * yi +     a2 * yi2  +     a3 * yi2 * yi
         Yprimin  =       a1      + 2.* a2 *  yi  + 3.* a3 * yi2
500   CONTINUE
        WRITE(6,*) ' Pas de solution ***** ',j,ylon2,iter
         STOP 2
550   CONTINUE
c
       Yprimin   = a1  + 2.* a2 *  yi   + 3.* a3 * yi* yi
       yprim(j)  = pi / ( jjm * Yprimin )
       yvrai(j)  = yi 

1500    CONTINUE

       DO j = 1, jlat -1
        IF( yvrai(j+1). LT. yvrai(j) )  THEN
         WRITE(6,*) ' PBS. avec  rlat(',j+1,') plus petit que rlat(',j,
     ,  ')'
         STOP 3
        ENDIF
       ENDDO

       WRITE(6,*) 'Reorganisation des latitudes pour avoir entre - pi/2'
     , ,' et  pi/2 '
c
        IF( ik.EQ.1 )   THEN
           ypn = pis2 
          DO j = jlat,1,-1
           IF( yvrai(j).LE. ypn ) GO TO 1502
          ENDDO
1502     CONTINUE

         jpn   = j
         y00   = yvrai(jpn)
         deply = pis2 -  y00
        ENDIF

         DO  j = 1, jjm +1 - jpn
           ylatt (j)  = -pis2 - y00  + yvrai(jpn+j-1)
           yprimm(j)  = yprim(jpn+j-1)
         ENDDO

         jjpn  = jpn
         IF( jlat.EQ. jjm ) jjpn = jpn -1

         DO j = 1,jjpn 
          ylatt (j + jjm+1 -jpn) = yvrai(j) + deply
          yprimm(j + jjm+1 -jpn) = yprim(j)
         ENDDO

c      ***********   Fin de la reorganisation     *************
c
 1600   CONTINUE

       DO j = 1, jlat
          ylat(j) =  ylatt( jlat +1 -j )
         yprim(j) = yprimm( jlat +1 -j )
       ENDDO
  
        DO j = 1, jlat
         yvrai(j) = ylat(j)*180./pi
        ENDDO

        IF( ik.EQ.1 )  THEN
c         WRITE(6,18) 
c         WRITE(6,*)  ' YLAT  en U   apres ( en  deg. ) '
c         WRITE(6,68) (yvrai(j),j=1,jlat)
cc         WRITE(6,*) ' YPRIM '
cc         WRITE(6,445) ( yprim(j),j=1,jlat)

          DO j = 1, jlat
            rrlatu(j) =  ylat( j )
           yyprimu(j) = yprim( j )
          ENDDO

        ELSE IF ( ik.EQ. 2 )  THEN
c         WRITE(6,18) 
c         WRITE(6,*) ' YLAT   en V  apres ( en  deg. ) '
c         WRITE(6,68) (yvrai(j),j=1,jlat)
cc         WRITE(6,*)' YPRIM '
cc         WRITE(6,445) ( yprim(j),j=1,jlat)

          DO j = 1, jlat
            rrlatv(j) =  ylat( j )
           yyprimv(j) = yprim( j )
          ENDDO

        ELSE IF ( ik.EQ. 3 )  THEN
c         WRITE(6,18) 
c         WRITE(6,*)  ' YLAT  en U + 0.75  apres ( en  deg. ) '
c         WRITE(6,68) (yvrai(j),j=1,jlat)
cc         WRITE(6,*) ' YPRIM '
cc         WRITE(6,445) ( yprim(j),j=1,jlat)

          DO j = 1, jlat
            rlatu2(j) =  ylat( j )
           yprimu2(j) = yprim( j )
          ENDDO

        ELSE IF ( ik.EQ. 4 )  THEN
c         WRITE(6,18) 
c         WRITE(6,*)  ' YLAT en U + 0.25  apres ( en  deg. ) '
c         WRITE(6,68)(yvrai(j),j=1,jlat)
cc         WRITE(6,*) ' YPRIM '
cc         WRITE(6,68) ( yprim(j),j=1,jlat)

          DO j = 1, jlat
            rlatu1(j) =  ylat( j )
           yprimu1(j) = yprim( j )
          ENDDO

        ENDIF

5000   CONTINUE
c
        WRITE(6,18)
c
c  .....     fin de la boucle  do 5000 .....

        DO j = 1, jjm
         ylat(j) = rrlatu(j) - rrlatu(j+1)
        ENDDO
        champmin =  1.e12
        champmax = -1.e12
        DO j = 1, jjm
         champmin = MIN( champmin, ylat(j) )
         champmax = MAX( champmax, ylat(j) )
        ENDDO
         champmin = champmin * 180./pi
         champmax = champmax * 180./pi

24     FORMAT(2x,'Parametres yzoom,gross,tau ,dzoom pour fyhyp ',4f8.3)
18      FORMAT(/)
68      FORMAT(1x,7f9.2)

        RETURN
        END
