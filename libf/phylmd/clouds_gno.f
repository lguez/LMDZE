!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/clouds_gno.F,v 1.2 2004/11/09 16:55:40 lmdzadmin Exp $
!
C
C================================================================================
C
      SUBROUTINE CLOUDS_GNO(klon,ND,R,RS,QSUB,PTCONV,RATQSC,CLDF)
      use numer_rec_95, only: nr_erf
      IMPLICIT NONE
C     
C--------------------------------------------------------------------------------
C
C Inputs:
C
C  ND----------: Number of vertical levels
C  R--------ND-: Domain-averaged mixing ratio of total water 
C  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox
C  QSUB-----ND-: Mixing ratio of condensed water within clouds associated
C                with SUBGRID-SCALE condensation processes (here, it is
C                predicted by the convection scheme)
C Outputs:
C
C  PTCONV-----ND-: Point convectif = TRUE
C  RATQSC-----ND-: Largeur normalisee de la distribution
C  CLDF-----ND-: Fraction nuageuse
C
C--------------------------------------------------------------------------------


      INTEGER klon,ND
      REAL R(klon,ND),  RS(klon,ND), QSUB(klon,ND)
      LOGICAL PTCONV(klon,ND)
      REAL RATQSC(klon,ND)
      REAL CLDF(klon,ND)

c -- parameters controlling the iteration:
c --    nmax    : maximum nb of iterations (hopefully never reached)
c --    epsilon : accuracy of the numerical resolution 
c --    vmax    : v-value above which we use an asymptotic expression for ERF(v)

      INTEGER nmax
      PARAMETER ( nmax = 10) 
      REAL epsilon, vmax0, vmax(klon)
      PARAMETER ( epsilon = 0.02, vmax0 = 2.0 ) 

      REAL min_mu, min_Q
      PARAMETER ( min_mu =  1.e-12, min_Q=1.e-12 )
     
      INTEGER i,K, n, m
      REAL mu(klon), qsat(klon), delta(klon), beta(klon) 
      real zu2(klon),zv2(klon)
      REAL xx(klon), aux(klon), coeff(klon), block(klon)
      REAL  dist(klon), fprime(klon), det(klon)
      REAL pi, u(klon), v(klon), erfcu(klon), erfcv(klon)
      REAL  xx1(klon), xx2(klon)
      real kkk
      real sqrtpi,sqrt2,zx1,zx2,exdel
c lconv = true si le calcul a converge (entre autre si qsub < min_q)
       LOGICAL lconv(klon)

cym
      cldf(:,:)=0.0
      
      pi = ACOS(-1.)
      sqrtpi=sqrt(pi)
      sqrt2=sqrt(2.)

      ptconv=.false.
      ratqsc=0.


      DO 500 K = 1, ND

                                    do i=1,klon ! vector
      mu(i) = R(i,K)
      mu(i) = MAX(mu(i),min_mu)
      qsat(i) = RS(i,K) 
      qsat(i) = MAX(qsat(i),min_mu)
      delta(i) = log(mu(i)/qsat(i))
                                    enddo ! vector

C
C ***          There is no subgrid-scale condensation;        ***
C ***   the scheme becomes equivalent to an "all-or-nothing"  *** 
C ***             large-scale condensation scheme.            ***
C

C
C ***     Some condensation is produced at the subgrid-scale       ***
C ***                                                              ***
C ***       PDF = generalized log-normal distribution (GNO)        ***
C ***   (k<0 because a lower bound is considered for the PDF)      ***
C ***                                                              ***
C ***  -> Determine x (the parameter k of the GNO PDF) such        ***
C ***  that the contribution of subgrid-scale processes to         ***
C ***  the in-cloud water content is equal to QSUB(K)              ***
C ***  (equations (13), (14), (15) + Appendix B of the paper)      ***
C ***                                                              ***
C ***    Here, an iterative method is used for this purpose        ***
C ***    (other numerical methods might be more efficient)         ***
C ***                                                              ***
C ***          NB: the "error function" is called ERF              ***
C ***                 (ERF in double precision)                   ***
C

c  On commence par eliminer les cas pour lesquels on n'a pas
c  suffisamment d'eau nuageuse.

                                    do i=1,klon ! vector

      IF ( QSUB(i,K) .lt. min_Q ) THEN
        ptconv(i,k)=.false.
        ratqsc(i,k)=0.
        lconv(i)  = .true.

c   Rien on a deja initialise

      ELSE 

        lconv(i)  = .FALSE. 
        vmax(i) = vmax0

        beta(i) = QSUB(i,K)/mu(i) + EXP( -MIN(0.0,delta(i)) )

c --  roots of equation v > vmax:

        det(i) = delta(i) + vmax(i)**2.
        if (det(i).LE.0.0) vmax(i) = vmax0 + 1.0
        det(i) = delta(i) + vmax(i)**2.

        if (det(i).LE.0.) then
          xx(i) = -0.0001
        else 
         zx1=-sqrt2*vmax(i)
         zx2=SQRT(1.0+delta(i)/(vmax(i)**2.))
         xx1(i)=zx1*(1.0-zx2)
         xx2(i)=zx1*(1.0+zx2)
         xx(i) = 1.01 * xx1(i)
         if ( xx1(i) .GE. 0.0 ) xx(i) = 0.5*xx2(i)
        endif
        if (delta(i).LT.0.) xx(i) = -0.5*SQRT(log(2.)) 

      ENDIF

                                    enddo       ! vector

c----------------------------------------------------------------------
c   Debut des nmax iterations pour trouver la solution.
c----------------------------------------------------------------------

      DO n = 1, nmax 

                                    do i=1,klon ! vector
        if (.not.lconv(i)) then

          u(i) = delta(i)/(xx(i)*sqrt2) + xx(i)/(2.*sqrt2)
          v(i) = delta(i)/(xx(i)*sqrt2) - xx(i)/(2.*sqrt2)

          IF ( v(i) .GT. vmax(i) ) THEN 

            IF (     ABS(u(i))  .GT. vmax(i) 
     :          .AND.  delta(i) .LT. 0. ) THEN

c -- use asymptotic expression of erf for u and v large:
c ( -> analytic solution for xx )
             exdel=beta(i)*EXP(delta(i))
             aux(i) = 2.0*delta(i)*(1.-exdel)
     :                       /(1.+exdel)
             if (aux(i).lt.0.) then
c                print*,'AUX(',i,',',k,')<0',aux(i),delta(i),beta(i)
                aux(i)=0.
             endif
             xx(i) = -SQRT(aux(i))
             block(i) = EXP(-v(i)*v(i)) / v(i) / sqrtpi
             dist(i) = 0.0
             fprime(i) = 1.0

            ELSE

c -- erfv -> 1.0, use an asymptotic expression of erfv for v large:

             erfcu(i) = 1.0-NR_ERF(u(i))
c  !!! ATTENTION : rajout d'un seuil pour l'exponentiel
             aux(i) = sqrtpi*erfcu(i)*EXP(min(v(i)*v(i),100.))
             coeff(i) = 1.0 - 1./2./(v(i)**2.) + 3./4./(v(i)**4.)
             block(i) = coeff(i) * EXP(-v(i)*v(i)) / v(i) / sqrtpi
             dist(i) = v(i) * aux(i) / coeff(i) - beta(i)
             fprime(i) = 2.0 / xx(i) * (v(i)**2.)
     :           * ( coeff(i)*EXP(-delta(i)) - u(i) * aux(i) )
     :           / coeff(i) / coeff(i)
            
            ENDIF ! ABS(u)

          ELSE

c -- general case:

           erfcu(i) = 1.0-NR_ERF(u(i))
           erfcv(i) = 1.0-NR_ERF(v(i))
           block(i) = erfcv(i)
           dist(i) = erfcu(i) / erfcv(i) - beta(i)
           zu2(i)=u(i)*u(i)
           zv2(i)=v(i)*v(i)
           if(zu2(i).gt.20..or. zv2(i).gt.20.) then
c              print*,'ATTENTION !!! xx(',i,') =', xx(i)
c           print*,'ATTENTION !!! klon,ND,R,RS,QSUB,PTCONV,RATQSC,CLDF',
c     .klon,ND,R(i,k),RS(i,k),QSUB(i,k),PTCONV(i,k),RATQSC(i,k),
c     .CLDF(i,k)
c              print*,'ATTENTION !!! zu2 zv2 =',zu2(i),zv2(i)
              zu2(i)=20.
              zv2(i)=20.
             fprime(i) = 0.
           else
             fprime(i) = 2. /sqrtpi /xx(i) /erfcv(i)**2.
     :           * (   erfcv(i)*v(i)*EXP(-zu2(i))
     :               - erfcu(i)*u(i)*EXP(-zv2(i)) )
           endif
          ENDIF ! x

c -- test numerical convergence:

c         print*,'avant test ',i,k,lconv(i),u(i),v(i)
          if ( ABS(dist(i)/beta(i)) .LT. epsilon ) then 
c           print*,'v-u **2',(v(i)-u(i))**2
c           print*,'exp v-u **2',exp((v(i)-u(i))**2)
            ptconv(i,K) = .TRUE. 
            lconv(i)=.true.
c  borne pour l'exponentielle
            ratqsc(i,k)=min(2.*(v(i)-u(i))**2,20.)
            ratqsc(i,k)=sqrt(exp(ratqsc(i,k))-1.)
            CLDF(i,K) = 0.5 * block(i)
          else
            xx(i) = xx(i) - dist(i)/fprime(i)
          endif
c         print*,'apres test ',i,k,lconv(i)

        endif ! lconv
                                    enddo       ! vector

c----------------------------------------------------------------------
c   Fin des nmax iterations pour trouver la solution.
        ENDDO ! n
c----------------------------------------------------------------------

500   CONTINUE  ! K

       RETURN
       END



