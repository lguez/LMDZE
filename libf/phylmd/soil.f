!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/soil.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE soil(ptimestep, indice, knon, snow, ptsrf, ptsoil,
     s          pcapcal, pfluxgrd)
      use dimens_m
      use indicesol
      use dimphy
      use dimsoil
      use YOMCST
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  Frederic Hourdin     30/01/92
c   -------
c
c   objet:  computation of : the soil temperature evolution
c   ------                   the surfacic heat capacity "Capcal"
c                            the surface conduction flux pcapcal
c
c
c   Method: implicit time integration
c   -------
c   Consecutive ground temperatures are related by:
c           T(k+1) = C(k) + D(k)*T(k)  (1)
c   the coefficients C and D are computed at the t-dt time-step.
c   Routine structure:
c   1)new temperatures are computed  using (1)
c   2)C and D coefficients are computed from the new temperature
c     profile for the t+dt time-step
c   3)the coefficients A and B are computed where the diffusive
c     fluxes at the t+dt time-step is given by
c            Fdiff = A + B Ts(t+dt)
c     or     Fdiff = F0 + Capcal (Ts(t+dt)-Ts(t))/dt
c            with F0 = A + B (Ts(t))
c                 Capcal = B*dt
c           
c   Interface:
c   ----------
c
c   Arguments:
c   ----------
c   ptimestep            physical timestep (s)
c   indice               sub-surface index
c   snow(klon,nbsrf)     snow
c   ptsrf(klon)          surface temperature at time-step t (K)
c   ptsoil(klon,nsoilmx) temperature inside the ground (K)
c   pcapcal(klon)        surfacic specific heat (W*m-2*s*K-1)
c   pfluxgrd(klon)       surface diffusive flux from ground (Wm-2)
c   
c=======================================================================
c   declarations:
c   -------------


c-----------------------------------------------------------------------
c  arguments
c  ---------

      REAL ptimestep
      INTEGER indice, knon
      REAL ptsrf(klon),ptsoil(klon,nsoilmx),snow(klon)
      REAL pcapcal(klon),pfluxgrd(klon)

c-----------------------------------------------------------------------
c  local arrays
c  ------------

      INTEGER ig,jk
c$$$      REAL zdz2(nsoilmx),z1(klon)
      REAL zdz2(nsoilmx),z1(klon,nbsrf)
      REAL min_period,dalph_soil
      REAL ztherm_i(klon)

c   local saved variables:
c   ----------------------
      REAL dz1(nsoilmx),dz2(nsoilmx)
c$$$          REAL zc(klon,nsoilmx),zd(klon,nsoilmx)
      REAL zc(klon,nsoilmx,nbsrf),zd(klon,nsoilmx,nbsrf)
      REAL lambda
      SAVE dz1,dz2,zc,zd,lambda
      LOGICAL firstcall, firstsurf(nbsrf)
      SAVE firstcall, firstsurf
      REAL isol,isno,iice
      SAVE isol,isno,iice

      DATA firstcall/.true./
      DATA firstsurf/.TRUE.,.TRUE.,.TRUE.,.TRUE./

      DATA isol,isno,iice/2000.,2000.,2000./

c-----------------------------------------------------------------------
c   Depthts:
c   --------

      REAL fz,rk,fz1,rk1,rk2
      fz(rk)=fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)
      pfluxgrd(:) = 0.
c   calcul de l'inertie thermique a partir de la variable rnat.
c   on initialise a iice meme au-dessus d'un point de mer au cas 
c   ou le point de mer devienne point de glace au pas suivant
c   on corrige si on a un point de terre avec ou sans glace
c
      IF (indice.EQ.is_sic) THEN
         DO ig = 1, knon
            ztherm_i(ig)   = iice
            IF (snow(ig).GT.0.0) ztherm_i(ig)   = isno
         ENDDO
      ELSE IF (indice.EQ.is_lic) THEN
         DO ig = 1, knon
            ztherm_i(ig)   = iice
            IF (snow(ig).GT.0.0) ztherm_i(ig)   = isno
         ENDDO
      ELSE IF (indice.EQ.is_ter) THEN
         DO ig = 1, knon
            ztherm_i(ig)   = isol
            IF (snow(ig).GT.0.0) ztherm_i(ig)   = isno
         ENDDO
      ELSE IF (indice.EQ.is_oce) THEN
         DO ig = 1, knon
            ztherm_i(ig)   = iice
         ENDDO
      ELSE
         PRINT*, "valeur d indice non prevue", indice
         stop 1
      ENDIF


c$$$      IF (firstcall) THEN
      IF (firstsurf(indice)) THEN 

c-----------------------------------------------------------------------
c   ground levels 
c   grnd=z/l where l is the skin depth of the diurnal cycle:
c   --------------------------------------------------------

         min_period=1800. ! en secondes
         dalph_soil=2.    ! rapport entre les epaisseurs de 2 couches succ.

         OPEN(99,file='soil.def',status='old',form='formatted',err=9999)
         READ(99,*) min_period
         READ(99,*) dalph_soil
         PRINT*,'Discretization for the soil model'
         PRINT*,'First level e-folding depth',min_period,
     s   '   dalph',dalph_soil
         CLOSE(99)
9999     CONTINUE

c   la premiere couche represente un dixieme de cycle diurne
         fz1=sqrt(min_period/3.14)

         DO jk=1,nsoilmx
            rk1=jk
            rk2=jk-1
            dz2(jk)=fz(rk1)-fz(rk2)
         ENDDO
         DO jk=1,nsoilmx-1
            rk1=jk+.5
            rk2=jk-.5
            dz1(jk)=1./(fz(rk1)-fz(rk2))
         ENDDO
         lambda=fz(.5)*dz1(1)
         PRINT*,'full layers, intermediate layers (seconds)'
         DO jk=1,nsoilmx
            rk=jk
            rk1=jk+.5
            rk2=jk-.5
            PRINT *,'fz=',
     .               fz(rk1)*fz(rk2)*3.14,fz(rk)*fz(rk)*3.14
         ENDDO
C PB
         firstsurf(indice) = .FALSE. 
c$$$         firstcall =.false.

c   Initialisations:
c   ----------------

      ELSE   !--not firstcall
c-----------------------------------------------------------------------
c   Computation of the soil temperatures using the Cgrd and Dgrd
c  coefficient computed at the previous time-step:
c  -----------------------------------------------

c    surface temperature
         DO ig=1,knon
            ptsoil(ig,1)=(lambda*zc(ig,1,indice)+ptsrf(ig))/
     s      (lambda*(1.-zd(ig,1,indice))+1.)
         ENDDO

c   other temperatures
         DO jk=1,nsoilmx-1
            DO ig=1,knon
               ptsoil(ig,jk+1)=zc(ig,jk,indice)+zd(ig,jk,indice)
     $            *ptsoil(ig,jk)
            ENDDO
         ENDDO

      ENDIF !--not firstcall
c-----------------------------------------------------------------------
c   Computation of the Cgrd and Dgrd coefficient for the next step:
c   ---------------------------------------------------------------

c$$$  PB ajout pour cas glace de mer
      IF (indice .EQ. is_sic) THEN
          DO ig = 1 , knon
            ptsoil(ig,nsoilmx) = RTT - 1.8
          END DO 
      ENDIF 

      DO jk=1,nsoilmx
         zdz2(jk)=dz2(jk)/ptimestep
      ENDDO

      DO ig=1,knon
         z1(ig,indice)=zdz2(nsoilmx)+dz1(nsoilmx-1)
         zc(ig,nsoilmx-1,indice)=
     $       zdz2(nsoilmx)*ptsoil(ig,nsoilmx)/z1(ig,indice)
         zd(ig,nsoilmx-1,indice)=dz1(nsoilmx-1)/z1(ig,indice)
      ENDDO

      DO jk=nsoilmx-1,2,-1
         DO ig=1,knon
            z1(ig,indice)=1./(zdz2(jk)+dz1(jk-1)+dz1(jk)
     $         *(1.-zd(ig,jk,indice)))
            zc(ig,jk-1,indice)=
     s      (ptsoil(ig,jk)*zdz2(jk)+dz1(jk)*zc(ig,jk,indice))
     $          *z1(ig,indice)
            zd(ig,jk-1,indice)=dz1(jk-1)*z1(ig,indice)
         ENDDO
      ENDDO

c-----------------------------------------------------------------------
c   computation of the surface diffusive flux from ground and
c   calorific capacity of the ground:
c   ---------------------------------

      DO ig=1,knon
         pfluxgrd(ig)=ztherm_i(ig)*dz1(1)*
     s   (zc(ig,1,indice)+(zd(ig,1,indice)-1.)*ptsoil(ig,1))
         pcapcal(ig)=ztherm_i(ig)*
     s   (dz2(1)+ptimestep*(1.-zd(ig,1,indice))*dz1(1))
         z1(ig,indice)=lambda*(1.-zd(ig,1,indice))+1.
         pcapcal(ig)=pcapcal(ig)/z1(ig,indice)
         pfluxgrd(ig) = pfluxgrd(ig)
     s   + pcapcal(ig) * (ptsoil(ig,1) * z1(ig,indice)
     $       - lambda * zc(ig,1,indice)
     $       - ptsrf(ig))
     s   /ptimestep
      ENDDO

      RETURN
      END
