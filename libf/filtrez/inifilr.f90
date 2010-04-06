SUBROUTINE inifilr

  ! From filtrez/inifilr.F,v 1.1.1.1 2004/05/19 12:53:09
  ! H. Upadhyaya, O.Sharma

  !   This routine computes the eigenfunctions of the laplacien           
  !   on the stretched grid, and the filtering coefficients               

  !  We designate:                                                        
  !   eignfn   eigenfunctions of the discrete laplacien                   
  !   eigenvl  eigenvalues                                                
  !   jfiltn   indexof the last scalar line filtered in NH                
  !   jfilts   index of the first line filtered in SH                     
  !   modfrst  index of the mode from where modes are filtered            
  !   modemax  maximum number of modes ( im )                             
  !   coefil   filtering coefficients ( lamda_max*cos(rlat)/lamda )       
  !   sdd      SQRT( dx )                                                 

  !     the modes are filtered from modfrst to modemax                    

  USE dimens_m
  USE paramet_m
  USE logic
  USE comgeom
  USE serre
  USE parafilt
  USE coefils

  IMPLICIT NONE

  REAL dlonu(iim), dlatu(jjm)
  REAL rlamda(iim), eignvl(iim)

  REAL lamdamax, pi, cof
  INTEGER i, j, modemax, imx, k, kf, ii
  REAL dymin, dxmin, colat0
  REAL eignft(iim,iim), coff
  EXTERNAL inifgn

  !-----------------------------------------------------------            


  pi = 2.*asin(1.)

  DO i = 1, iim
     dlonu(i) = xprimu(i)
  END DO

  CALL inifgn(eignvl)

  PRINT *, ' EIGNVL '
  PRINT 250, eignvl
250 FORMAT (1X,5E13.6)

  ! compute eigenvalues and eigenfunctions                                


  !.................................................................      

  !  compute the filtering coefficients for scalar lines and              
  !  meridional wind v-lines                                              

  !  we filter all those latitude lines where coefil < 1                  
  !  NO FILTERING AT POLES                                                

  !  colat0 is to be used  when alpha (stretching coefficient)            
  !  is set equal to zero for the regular grid case                       

  !    .......   Calcul  de  colat0   .........                           
  !     .....  colat0 = minimum de ( 0.5, min dy/ min dx )   ...          


  DO  j = 1, jjm
     dlatu(j) = rlatu(j) - rlatu(j+1)
  END DO

  dxmin = dlonu(1)
  DO i = 2, iim
     dxmin = min(dxmin,dlonu(i))
  END DO
  dymin = dlatu(1)
  DO j = 2, jjm
     dymin = min(dymin,dlatu(j))
  END DO


  colat0 = min(0.5,dymin/dxmin)

  IF ( .NOT. fxyhypb .AND. ysinus) THEN
     colat0 = 0.6
     !         ...... a revoir  pour  ysinus !   .......                     
     alphax = 0.
  END IF

  PRINT 50, colat0, alphax
50 FORMAT (/15X,' Inifilr colat0 alphax ',2E16.7)

  IF (alphax==1.) THEN
     PRINT *, ' Inifilr  alphax doit etre  <  a 1.  Corriger '
     STOP 1
  END IF

  lamdamax = iim/(pi*colat0*(1.-alphax))

  DO  i = 2, iim
     rlamda(i) = lamdamax/sqrt(abs(eignvl(i)))
  END DO


  DO  j = 1, jjm
     DO  i = 1, iim
        coefilu(i,j) = 0.0
        coefilv(i,j) = 0.0
        coefilu2(i,j) = 0.0
        coefilv2(i,j) = 0.0
     end DO
  END DO


  !    ... Determination de jfiltnu,jfiltnv,jfiltsu,jfiltsv ....          
  !    .........................................................          

  modemax = iim

  !ccc    imx = modemax - 4 * (modemax/iim)                               

  imx = iim

  PRINT *, ' TRUNCATION AT ', imx

  DO  j = 2, jjm/2 + 1
     cof = cos(rlatu(j))/colat0
     IF (cof<1.) THEN
        IF (rlamda(imx)*cos(rlatu(j))<1.) jfiltnu = j
     END IF

     cof = cos(rlatu(jjp1-j+1))/colat0
     IF (cof<1.) THEN
        IF (rlamda(imx)*cos(rlatu(jjp1-j+1))<1.) jfiltsu = jjp1 - j + 1
     END IF
  END DO

  DO  j = 1, jjm/2
     cof = cos(rlatv(j))/colat0
     IF (cof<1.) THEN
        IF (rlamda(imx)*cos(rlatv(j))<1.) jfiltnv = j
     END IF

     cof = cos(rlatv(jjm-j+1))/colat0
     IF (cof<1.) THEN
        IF (rlamda(imx)*cos(rlatv(jjm-j+1))<1.) jfiltsv = jjm - j + 1
     END IF
  END DO


  IF (jfiltnu<=0) jfiltnu = 1
  IF (jfiltnu>jjm/2+1) THEN
     PRINT *, ' jfiltnu en dehors des valeurs acceptables ', jfiltnu
     STOP 1
  END IF

  IF (jfiltsu<=0) jfiltsu = 1
  IF (jfiltsu>jjm+1) THEN
     PRINT *, ' jfiltsu en dehors des valeurs acceptables ', jfiltsu
     STOP 1
  END IF

  IF (jfiltnv<=0) jfiltnv = 1
  IF (jfiltnv>jjm/2) THEN
     PRINT *, ' jfiltnv en dehors des valeurs acceptables ', jfiltnv
     STOP 1
  END IF

  IF (jfiltsv<=0) jfiltsv = 1
  IF (jfiltsv>jjm) THEN
     PRINT *, ' jfiltsv en dehors des valeurs acceptables ', jfiltsv
     STOP 1
  END IF

  PRINT *, ' jfiltnv jfiltsv jfiltnu jfiltsu ', jfiltnv, jfiltsv, jfiltnu, &
       jfiltsu


  !   ... Determination de coefilu,coefilv,n=modfrstu,modfrstv ....       
  !................................................................       


  DO  j = 1, jjm
     modfrstu(j) = iim
     modfrstv(j) = iim
  END DO

  DO  j = 2, jfiltnu
     DO  k = 2, modemax
        cof = rlamda(k)*cos(rlatu(j))
        IF (cof<1.) GO TO 82
     end DO
     cycle
82   modfrstu(j) = k

     kf = modfrstu(j)
     DO  k = kf, modemax
        cof = rlamda(k)*cos(rlatu(j))
        coefilu(k,j) = cof - 1.
        coefilu2(k,j) = cof*cof - 1.
     end DO
  END DO


  DO j = 1, jfiltnv
     DO  k = 2, modemax
        cof = rlamda(k)*cos(rlatv(j))
        IF (cof<1.) GO TO 87
     end DO
     cycle
87   modfrstv(j) = k

     kf = modfrstv(j)
     DO  k = kf, modemax
        cof = rlamda(k)*cos(rlatv(j))
        coefilv(k,j) = cof - 1.
        coefilv2(k,j) = cof*cof - 1.
     end DO
  end DO

  DO  j = jfiltsu, jjm
     DO  k = 2, modemax
        cof = rlamda(k)*cos(rlatu(j))
        IF (cof<1.) GO TO 92
     end DO
     cycle
92   modfrstu(j) = k

     kf = modfrstu(j)
     DO  k = kf, modemax
        cof = rlamda(k)*cos(rlatu(j))
        coefilu(k,j) = cof - 1.
        coefilu2(k,j) = cof*cof - 1.
     end DO
  end DO

  DO  j = jfiltsv, jjm
     DO  k = 2, modemax
        cof = rlamda(k)*cos(rlatv(j))
        IF (cof<1.) GO TO 97
     end DO
     cycle
97   modfrstv(j) = k

     kf = modfrstv(j)
     DO  k = kf, modemax
        cof = rlamda(k)*cos(rlatv(j))
        coefilv(k,j) = cof - 1.
        coefilv2(k,j) = cof*cof - 1.
     end DO
  END DO


  IF (jfiltnv>=jjm/2 .OR. jfiltnu>=jjm/2) THEN

     IF (jfiltnv==jfiltsv) jfiltsv = 1 + jfiltnv
     IF (jfiltnu==jfiltsu) jfiltsu = 1 + jfiltnu

     PRINT *, 'jfiltnv jfiltsv jfiltnu jfiltsu', jfiltnv, jfiltsv, jfiltnu, &
          jfiltsu
  END IF

  PRINT *, '   Modes premiers  v  '
  PRINT 334, modfrstv
  PRINT *, '   Modes premiers  u  '
  PRINT 334, modfrstu


  IF (nfilun<jfiltnu) THEN
     PRINT *, ' le parametre nfilun utilise pour la matrice ', &
          ' matriceun  est trop petit ! '
     PRINT *, 'Le changer dans parafilt.h et le mettre a  ', jfiltnu
     PRINT *, 'Pour information, nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
     STOP 1
  END IF
  IF (nfilun>jfiltnu+2) THEN
     PRINT *, ' le parametre nfilun utilise pour la matrice ', &
          ' matriceun est trop grand ! Gachis de memoire ! '
     PRINT *, 'Le changer dans parafilt.h et le mettre a  ', jfiltnu
     PRINT *, 'Pour information, nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
  END IF
  IF (nfilus<jjm-jfiltsu+1) THEN
     PRINT *, ' le parametre nfilus utilise pour la matrice ', &
          ' matriceus  est trop petit !  '
     PRINT *, ' Le changer dans parafilt.h et le mettre a  ', &
          jjm - jfiltsu + 1
     PRINT *, ' Pour information , nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
     STOP 1
  END IF
  IF (nfilus>jjm-jfiltsu+3) THEN
     PRINT *, ' le parametre nfilus utilise pour la matrice ', &
          ' matriceus  est trop grand ! '
     PRINT *, ' Le changer dans parafilt.h et le mettre a  ', &
          jjm - jfiltsu + 1
     PRINT *, ' Pour information , nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
  END IF
  IF (nfilvn<jfiltnv) THEN
     PRINT *, ' le parametre nfilvn utilise pour la matrice ', &
          ' matricevn  est trop petit ! '
     PRINT *, 'Le changer dans parafilt.h et le mettre a  ', jfiltnv
     PRINT *, ' Pour information , nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
     STOP 1
  END IF
  IF (nfilvn>jfiltnv+2) THEN
     PRINT *, ' le parametre nfilvn utilise pour la matrice ', &
          ' matricevn est trop grand !  Gachis de memoire ! '
     PRINT *, 'Le changer dans parafilt.h et le mettre a  ', jfiltnv
     PRINT *, ' Pour information , nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
  END IF
  IF (nfilvs<jjm-jfiltsv+1) THEN
     PRINT *, ' le parametre nfilvs utilise pour la matrice ', &
          ' matricevs  est trop petit !  Le changer dans parafilt.h '
     PRINT *, ' Le changer dans parafilt.h et le mettre a  ', &
          jjm - jfiltsv + 1
     PRINT *, ' Pour information , nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
     STOP 1
  END IF
  IF (nfilvs>jjm-jfiltsv+3) THEN
     PRINT *, ' le parametre nfilvs utilise pour la matrice ', &
          ' matricevs  est trop grand ! Gachis de memoire ! '
     PRINT *, ' Le changer dans parafilt.h et le mettre a  ', &
          jjm - jfiltsv + 1
     PRINT *, ' Pour information , nfilun,nfilus,nfilvn,nfilvs ', &
          'doivent etre egaux successivement a ', jfiltnu, jjm - jfiltsu + 1, &
          jfiltnv, jjm - jfiltsv + 1
  END IF

  !   ... Calcul de la matrice filtre 'matriceu'  pour les champs situes  
  !                       sur la grille scalaire                 ........ 

  DO j = 2, jfiltnu

     DO i = 1, iim
        coff = coefilu(i,j)
        IF (i<modfrstu(j)) coff = 0.
        DO k = 1, iim
           eignft(i,k) = eignfnv(k,i)*coff
        END DO
     END DO
     DO k = 1, iim
        DO i = 1, iim
           matriceun(i,k,j) = 0.0
           DO ii = 1, iim
              matriceun(i,k,j) = matriceun(i,k,j) + eignfnv(i,ii)*eignft(ii,k)
           END DO
        END DO
     END DO

  END DO

  DO j = jfiltsu, jjm

     DO i = 1, iim
        coff = coefilu(i,j)
        IF (i<modfrstu(j)) coff = 0.
        DO k = 1, iim
           eignft(i,k) = eignfnv(k,i)*coff
        END DO
     END DO
     DO k = 1, iim
        DO i = 1, iim
           matriceus(i,k,j-jfiltsu+1) = 0.0
           DO ii = 1, iim
              matriceus(i,k,j-jfiltsu+1) = matriceus(i,k,j-jfiltsu+1) + &
                   eignfnv(i,ii)*eignft(ii,k)
           END DO
        END DO
     END DO

  END DO

  !   ................................................................... 

  !   ... Calcul de la matrice filtre 'matricev'  pour les champs situes  
  !                       sur la grille   de V ou de Z           ........ 
  !   ................................................................... 

  DO j = 1, jfiltnv

     DO i = 1, iim
        coff = coefilv(i,j)
        IF (i<modfrstv(j)) coff = 0.
        DO k = 1, iim
           eignft(i,k) = eignfnu(k,i)*coff
        END DO
     END DO
     DO k = 1, iim
        DO i = 1, iim
           matricevn(i,k,j) = 0.0
           DO ii = 1, iim
              matricevn(i,k,j) = matricevn(i,k,j) + eignfnu(i,ii)*eignft(ii,k)
           END DO
        END DO
     END DO

  END DO

  DO j = jfiltsv, jjm

     DO i = 1, iim
        coff = coefilv(i,j)
        IF (i<modfrstv(j)) coff = 0.
        DO k = 1, iim
           eignft(i,k) = eignfnu(k,i)*coff
        END DO
     END DO
     DO k = 1, iim
        DO i = 1, iim
           matricevs(i,k,j-jfiltsv+1) = 0.0
           DO ii = 1, iim
              matricevs(i,k,j-jfiltsv+1) = matricevs(i,k,j-jfiltsv+1) + &
                   eignfnu(i,ii)*eignft(ii,k)
           END DO
        END DO
     END DO

  END DO

  !   ................................................................... 

  !   ... Calcul de la matrice filtre 'matrinv'  pour les champs situes   
  !              sur la grille scalaire , pour le filtre inverse ........ 
  !   ................................................................... 

  DO j = 2, jfiltnu

     DO i = 1, iim
        coff = coefilu(i,j)/(1.+coefilu(i,j))
        IF (i<modfrstu(j)) coff = 0.
        DO k = 1, iim
           eignft(i,k) = eignfnv(k,i)*coff
        END DO
     END DO
     DO k = 1, iim
        DO i = 1, iim
           matrinvn(i,k,j) = 0.0
           DO ii = 1, iim
              matrinvn(i,k,j) = matrinvn(i,k,j) + eignfnv(i,ii)*eignft(ii,k)
           END DO
        END DO
     END DO

  END DO

  DO j = jfiltsu, jjm

     DO i = 1, iim
        coff = coefilu(i,j)/(1.+coefilu(i,j))
        IF (i<modfrstu(j)) coff = 0.
        DO k = 1, iim
           eignft(i,k) = eignfnv(k,i)*coff
        END DO
     END DO
     DO k = 1, iim
        DO i = 1, iim
           matrinvs(i,k,j-jfiltsu+1) = 0.0
           DO ii = 1, iim
              matrinvs(i,k,j-jfiltsu+1) = matrinvs(i,k,j-jfiltsu+1) + &
                   eignfnv(i,ii)*eignft(ii,k)
           END DO
        END DO
     END DO

  END DO

334 FORMAT (1X,24I3)
755 FORMAT (1X,6F10.3,I3)

END SUBROUTINE inifilr
