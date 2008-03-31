      subroutine aaam_bud (iam,nlon,nlev,rsec,
     i                   rea,rg,ome,      
     i                   plat,plon,phis,
     i                   dragu,liftu,phyu,
     i                   dragv,liftv,phyv,
     i                   p, u, v,
     o                   aam, torsfc)
c
      use dimens_m
      use dimphy
      implicit none
c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 20031020
c Object: Compute different terms of the axial AAAM Budget.
C No outputs, every AAM quantities are written on the IAM
C File. 
c
c Modif : I.Musat (LMD/CNRS) date : 20041020
c Outputs : axial components of wind AAM "aam" and total surface torque "torsfc",
c but no write in the iam file.
c
C WARNING: Only valid for regular rectangular grids.
C REMARK: CALL DANS PHYSIQ AFTER lift_noro:
C        CALL aaam_bud (27,klon,klev,rjourvrai,gmtime,
C    C               ra,rg,romega,
C    C               rlat,rlon,pphis,
C    C               zustrdr,zustrli,zustrph,
C    C               zvstrdr,zvstrli,zvstrph,
C    C               paprs,u,v)
C
C======================================================================
c Explicit Arguments:
c ==================
c iam-----input-I-File number where AAMs and torques are written
c                 It is a formatted file that has been opened
c                 in physiq.F
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c rsec----input-R-Seconde de la journee
c rea-----input-R-Earth radius
c rg------input-R-gravity constant
c ome-----input-R-Earth rotation rate
c plat ---input-R-Latitude en degres
c plon ---input-R-Longitude en degres
c phis ---input-R-Geopotential at the ground
c dragu---input-R-orodrag stress (zonal)
c liftu---input-R-orolift stress (zonal)
c phyu----input-R-Stress total de la physique (zonal)
c dragv---input-R-orodrag stress (Meridional)
c liftv---input-R-orolift stress (Meridional)
c phyv----input-R-Stress total de la physique (Meridional)
c p-------input-R-Pressure (Pa) at model half levels
c u-------input-R-Horizontal wind (m/s)
c v-------input-R-Meridional wind (m/s)
c aam-----output-R-Axial Wind AAM (=raam(3))
c torsfc--output-R-Total surface torque (=tmou(3)+tsso(3)+tbls(3))
c
c Implicit Arguments:
c ===================
c
c iim--common-I: Number of longitude intervals
c jjm--common-I: Number of latitude intervals
c klon-common-I: Number of points seen by the physics
c                iim*(jjm-1)+2 for instance
c klev-common-I: Number of vertical layers
c======================================================================
c Local Variables:
c ================
c dlat-----R: Latitude increment (Radians)
c dlon-----R: Longitude increment (Radians)
c raam  ---R: Wind AAM (3 Components, 1 & 2 Equatoriales; 3 Axiale)
c oaam  ---R: Mass AAM (3 Components, 1 & 2 Equatoriales; 3 Axiale)
c tmou-----R: Resolved Mountain torque (3 components)
c tsso-----R: Parameterised Moutain drag torque (3 components)
c tbls-----R: Parameterised Boundary layer torque (3 components)
c
c LOCAL ARRAY:
c ===========
c zs    ---R: Topographic height
c ps    ---R: Surface Pressure  
c ub    ---R: Barotropic wind zonal
c vb    ---R: Barotropic wind meridional
c zlat  ---R: Latitude in radians
c zlon  ---R: Longitude in radians
c======================================================================

c
c ARGUMENTS
c
      INTEGER iam,nlon,nlev
      real, intent(in):: rsec
      real rea
      real, intent(in):: rg
      real ome
      REAL, intent(in):: plat(nlon),plon(nlon)
      real phis(nlon)
      REAL dragu(nlon),liftu(nlon),phyu(nlon)             
      REAL dragv(nlon),liftv(nlon),phyv(nlon)             
      REAL, intent(in):: p(nlon,nlev+1)
      real u(nlon,nlev), v(nlon,nlev)
c
c Variables locales:
c
      INTEGER i,j,k,l
      REAL xpi,hadley,hadday
      REAL dlat,dlon
      REAL raam(3),oaam(3),tmou(3),tsso(3),tbls(3)
      integer iax
cIM ajout aam, torsfc
c aam = composante axiale du Wind AAM raam
c torsfc = composante axiale de (tmou+tsso+tbls)
      REAL aam, torsfc

      REAL ZS(801,401),PS(801,401)
      REAL UB(801,401),VB(801,401)
      REAL SSOU(801,401),SSOV(801,401)
      REAL BLSU(801,401),BLSV(801,401)
      REAL ZLON(801),ZLAT(401)
C
C  PUT AAM QUANTITIES AT ZERO:
C
      if(iim+1.gt.801.or.jjm+1.gt.401)then
      print *,' Pb de dimension dans aaam_bud'
      stop
      endif

      xpi=acos(-1.)
      hadley=1.e18
      hadday=1.e18*24.*3600.
      dlat=xpi/float(jjm)
      dlon=2.*xpi/float(iim) 
      
      do iax=1,3
      oaam(iax)=0.
      raam(iax)=0.
      tmou(iax)=0.
      tsso(iax)=0.
      tbls(iax)=0.
      enddo

C MOUNTAIN HEIGHT, PRESSURE AND BAROTROPIC WIND:

C North pole values (j=1):
 
      l=1

        ub(1,1)=0.
        vb(1,1)=0.
        do k=1,nlev
          ub(1,1)=ub(1,1)+u(l,k)*(p(l,k)-p(l,k+1))/rg
          vb(1,1)=vb(1,1)+v(l,k)*(p(l,k)-p(l,k+1))/rg
        enddo

          zlat(1)=plat(l)*xpi/180.

        do i=1,iim+1

          zs(i,1)=phis(l)/rg
          ps(i,1)=p(l,1)
          ub(i,1)=ub(1,1)                             
          vb(i,1)=vb(1,1)                             
          ssou(i,1)=dragu(l)+liftu(l)
          ssov(i,1)=dragv(l)+liftv(l)
          blsu(i,1)=phyu(l)-dragu(l)-liftu(l)
          blsv(i,1)=phyv(l)-dragv(l)-liftv(l)

        enddo


      do j = 2,jjm

C Values at Greenwich (Periodicity)

      zs(iim+1,j)=phis(l+1)/rg
      ps(iim+1,j)=p(l+1,1)
          ssou(iim+1,j)=dragu(l+1)+liftu(l+1)
          ssov(iim+1,j)=dragv(l+1)+liftv(l+1)
          blsu(iim+1,j)=phyu(l+1)-dragu(l+1)-liftu(l+1)
          blsv(iim+1,j)=phyv(l+1)-dragv(l+1)-liftv(l+1)
      zlon(iim+1)=-plon(l+1)*xpi/180.
      zlat(j)=plat(l+1)*xpi/180.

      ub(iim+1,j)=0.
      vb(iim+1,j)=0.
         do k=1,nlev
         ub(iim+1,j)=ub(iim+1,j)+u(l+1,k)*(p(l+1,k)-p(l+1,k+1))/rg
         vb(iim+1,j)=vb(iim+1,j)+v(l+1,k)*(p(l+1,k)-p(l+1,k+1))/rg
         enddo
      

      do i=1,iim

      l=l+1
      zs(i,j)=phis(l)/rg
      ps(i,j)=p(l,1)
          ssou(i,j)=dragu(l)+liftu(l)
          ssov(i,j)=dragv(l)+liftv(l)
          blsu(i,j)=phyu(l)-dragu(l)-liftu(l)
          blsv(i,j)=phyv(l)-dragv(l)-liftv(l)
      zlon(i)=plon(l)*xpi/180.

      ub(i,j)=0.
      vb(i,j)=0.
         do k=1,nlev
         ub(i,j)=ub(i,j)+u(l,k)*(p(l,k)-p(l,k+1))/rg
         vb(i,j)=vb(i,j)+v(l,k)*(p(l,k)-p(l,k+1))/rg
         enddo

      enddo

      enddo


C South Pole

      l=l+1
      ub(1,jjm+1)=0.
      vb(1,jjm+1)=0.
      do k=1,nlev
         ub(1,jjm+1)=ub(1,jjm+1)+u(l,k)*(p(l,k)-p(l,k+1))/rg
         vb(1,jjm+1)=vb(1,jjm+1)+v(l,k)*(p(l,k)-p(l,k+1))/rg
      enddo
      zlat(jjm+1)=plat(l)*xpi/180.

      do i=1,iim+1
      zs(i,jjm+1)=phis(l)/rg
      ps(i,jjm+1)=p(l,1)
          ssou(i,jjm+1)=dragu(l)+liftu(l)
          ssov(i,jjm+1)=dragv(l)+liftv(l)
          blsu(i,jjm+1)=phyu(l)-dragu(l)-liftu(l)
          blsv(i,jjm+1)=phyv(l)-dragv(l)-liftv(l)
      ub(i,jjm+1)=ub(1,jjm+1)                               
      vb(i,jjm+1)=vb(1,jjm+1)                                
      enddo

C
C  MOMENT ANGULAIRE 
C
        DO j=1,jjm    
        DO i=1,iim

           raam(1)=raam(1)-rea**3*dlon*dlat*0.5*
     c    (cos(zlon(i  ))*sin(zlat(j  ))*cos(zlat(j  ))*ub(i  ,j  )
     c    +cos(zlon(i  ))*sin(zlat(j+1))*cos(zlat(j+1))*ub(i  ,j+1))
     c                    +rea**3*dlon*dlat*0.5*
     c    (sin(zlon(i  ))*cos(zlat(j  ))*vb(i  ,j  )
     c    +sin(zlon(i  ))*cos(zlat(j+1))*vb(i  ,j+1))

           oaam(1)=oaam(1)-ome*rea**4*dlon*dlat/rg*0.5*
     c   (cos(zlon(i  ))*cos(zlat(j  ))**2*sin(zlat(j  ))*ps(i  ,j  )
     c   +cos(zlon(i  ))*cos(zlat(j+1))**2*sin(zlat(j+1))*ps(i  ,j+1))

           raam(2)=raam(2)-rea**3*dlon*dlat*0.5*
     c    (sin(zlon(i  ))*sin(zlat(j  ))*cos(zlat(j  ))*ub(i  ,j  )
     c    +sin(zlon(i  ))*sin(zlat(j+1))*cos(zlat(j+1))*ub(i  ,j+1))
     c                    -rea**3*dlon*dlat*0.5*
     c    (cos(zlon(i  ))*cos(zlat(j  ))*vb(i  ,j  )
     c    +cos(zlon(i  ))*cos(zlat(j+1))*vb(i  ,j+1))

           oaam(2)=oaam(2)-ome*rea**4*dlon*dlat/rg*0.5*
     c   (sin(zlon(i  ))*cos(zlat(j  ))**2*sin(zlat(j  ))*ps(i  ,j  )
     c   +sin(zlon(i  ))*cos(zlat(j+1))**2*sin(zlat(j+1))*ps(i  ,j+1))

           raam(3)=raam(3)+rea**3*dlon*dlat*0.5*
     c           (cos(zlat(j))**2*ub(i,j)+cos(zlat(j+1))**2*ub(i,j+1))

           oaam(3)=oaam(3)+ome*rea**4*dlon*dlat/rg*0.5*
     c        (cos(zlat(j))**3*ps(i,j)+cos(zlat(j+1))**3*ps(i,j+1))

        ENDDO
        ENDDO

C
C COUPLE DES MONTAGNES:
C

        DO j=1,jjm
        DO i=1,iim
           tmou(1)=tmou(1)-rea**2*dlon*0.5*sin(zlon(i))
     c  *(zs(i,j)-zs(i,j+1))
     c  *(cos(zlat(j+1))*ps(i,j+1)+cos(zlat(j))*ps(i,j)) 
           tmou(2)=tmou(2)+rea**2*dlon*0.5*cos(zlon(i))
     c  *(zs(i,j)-zs(i,j+1))
     c  *(cos(zlat(j+1))*ps(i,j+1)+cos(zlat(j))*ps(i,j)) 
        ENDDO
        ENDDO
           
        DO j=2,jjm 
        DO i=1,iim
           tmou(1)=tmou(1)+rea**2*dlat*0.5*sin(zlat(j))
     c  *(zs(i+1,j)-zs(i,j))
     c  *(cos(zlon(i+1))*ps(i+1,j)+cos(zlon(i))*ps(i,j))
           tmou(2)=tmou(2)+rea**2*dlat*0.5*sin(zlat(j))
     c  *(zs(i+1,j)-zs(i,j))
     c  *(sin(zlon(i+1))*ps(i+1,j)+sin(zlon(i))*ps(i,j))
           tmou(3)=tmou(3)-rea**2*dlat*0.5*
     c  cos(zlat(j))*(zs(i+1,j)-zs(i,j))*(ps(i+1,j)+ps(i,j))
        ENDDO
        ENDDO

C
C COUPLES DES DIFFERENTES FRICTION AU SOL:
C
        l=1
        DO j=2,jjm
        DO i=1,iim
        l=l+1
           tsso(1)=tsso(1)-rea**3*cos(zlat(j))*dlon*dlat*
     c     ssou(i,j)          *sin(zlat(j))*cos(zlon(i))
     c                    +rea**3*cos(zlat(j))*dlon*dlat*
     c     ssov(i,j)          *sin(zlon(i))

           tsso(2)=tsso(2)-rea**3*cos(zlat(j))*dlon*dlat*
     c     ssou(i,j)          *sin(zlat(j))*sin(zlon(i))
     c                    -rea**3*cos(zlat(j))*dlon*dlat*
     c     ssov(i,j)          *cos(zlon(i))

           tsso(3)=tsso(3)+rea**3*cos(zlat(j))*dlon*dlat*
     c     ssou(i,j)          *cos(zlat(j))

           tbls(1)=tbls(1)-rea**3*cos(zlat(j))*dlon*dlat*
     c     blsu(i,j)          *sin(zlat(j))*cos(zlon(i))
     c                    +rea**3*cos(zlat(j))*dlon*dlat*
     c     blsv(i,j)          *sin(zlon(i))

           tbls(2)=tbls(2)-rea**3*cos(zlat(j))*dlon*dlat*
     c     blsu(i,j)          *sin(zlat(j))*sin(zlon(i))
     c                    -rea**3*cos(zlat(j))*dlon*dlat*
     c     blsv(i,j)          *cos(zlon(i))

           tbls(3)=tbls(3)+rea**3*cos(zlat(j))*dlon*dlat*
     c     blsu(i,j)          *cos(zlat(j))

        ENDDO
        ENDDO
            

100   format(F12.5,15(1x,F12.5))

      aam=raam(3)
      torsfc= tmou(3)+tsso(3)+tbls(3)
c
      RETURN
      END
