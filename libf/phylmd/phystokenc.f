!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/phystokenc.F,v 1.2 2004/06/22 11:45:35 lmdzadmin Exp $
!
c
c
      SUBROUTINE phystokenc (
     I                   pdtphys,rlon,rlat,
     I                   pt,pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
     I                   pfm_therm,pentr_therm,
     I                   pcoefh,yu1,yv1,ftsol,pctsrf,
     I                   frac_impa,frac_nucl,
     I                   pphis,paire,dtime,itap)
      USE ioipsl
      use dimens_m
      use indicesol
      use dimphy
      use conf_gcm_m
      use tracstoke
      IMPLICIT none

c======================================================================
c Auteur(s) FH
c Objet: Moniteur general des tendances traceurs
c

c======================================================================
c======================================================================

c Arguments:
c
c   EN ENTREE:
c   ==========
c
c   divers:
c   -------
c
      real, intent(in):: pdtphys ! pas d'integration pour la physique (seconde)
c
      integer physid
      integer, intent(in):: itap
      save physid
      integer ndex2d(iim*(jjm+1)),ndex3d(iim*(jjm+1)*klev)

c   convection:
c   -----------
c
      REAL pmfu(klon,klev)  ! flux de masse dans le panache montant
      REAL pmfd(klon,klev)  ! flux de masse dans le panache descendant
      REAL pen_u(klon,klev) ! flux entraine dans le panache montant
      REAL pde_u(klon,klev) ! flux detraine dans le panache montant
      REAL pen_d(klon,klev) ! flux entraine dans le panache descendant
      REAL pde_d(klon,klev) ! flux detraine dans le panache descendant
	real pt(klon,klev),t(klon,klev)
c
      REAL, intent(in):: rlon(klon), rlat(klon)
      real, intent(in):: dtime
      REAL zx_tmp_3d(iim,jjm+1,klev),zx_tmp_2d(iim,jjm+1)

c   Couche limite:
c   --------------
c
      REAL pcoefh(klon,klev)    ! coeff melange CL
      REAL yv1(klon)
      REAL yu1(klon),pphis(klon),paire(klon)

c   Les Thermiques : (Abderr 25 11 02)
c   ---------------
      REAL pfm_therm(klon,klev+1)
	real fm_therm1(klon,klev)
      REAL pentr_therm(klon,klev)
      REAL entr_therm(klon,klev)
      REAL fm_therm(klon,klev)
c
c   Lessivage:
c   ----------
c
      REAL frac_impa(klon,klev)
      REAL frac_nucl(klon,klev)
c
c Arguments necessaires pour les sources et puits de traceur
C
      real ftsol(klon,nbsrf)  ! Temperature du sol (surf)(Kelvin)
      real pctsrf(klon,nbsrf) ! Pourcentage de sol f(nature du sol)
c======================================================================
c
      INTEGER i, k
c
      REAL mfu(klon,klev)  ! flux de masse dans le panache montant
      REAL mfd(klon,klev)  ! flux de masse dans le panache descendant
      REAL en_u(klon,klev) ! flux entraine dans le panache montant
      REAL de_u(klon,klev) ! flux detraine dans le panache montant
      REAL en_d(klon,klev) ! flux entraine dans le panache descendant
      REAL de_d(klon,klev) ! flux detraine dans le panache descendant
      REAL coefh(klon,klev) ! flux detraine dans le panache descendant

      REAL pyu1(klon),pyv1(klon)
      REAL pftsol(klon,nbsrf),ppsrf(klon,nbsrf)
      real pftsol1(klon),pftsol2(klon),pftsol3(klon),pftsol4(klon)
      real ppsrf1(klon),ppsrf2(klon),ppsrf3(klon),ppsrf4(klon)

      REAL dtcum

      integer iadvtr,irec
      real zmin,zmax
      logical ok_sync
 
      save t,mfu,mfd,en_u,de_u,en_d,de_d,coefh,dtcum
	save fm_therm,entr_therm
      save iadvtr,irec
      save pyu1,pyv1,pftsol,ppsrf

      data iadvtr,irec/0,1/
c
c   Couche limite:
c======================================================================

      ok_sync = .true.
	print*,'Dans phystokenc.F'
      print*,'iadvtr= ',iadvtr
      print*,'istphy= ',istphy
      print*,'istdyn= ',istdyn

      IF (iadvtr.eq.0) THEN
	
	CALL initphysto('phystoke',
     . rlon,rlat,dtime, dtime*istphy,dtime*istphy,nqmx,physid)
  	
	write(*,*) 'apres initphysto ds phystokenc' 

	
      ENDIF
c
      ndex2d = 0
      ndex3d = 0
      i=itap 
      CALL gr_fi_ecrit(1,klon,iim,jjm+1,pphis,zx_tmp_2d)
      CALL histwrite(physid,"phis",i,zx_tmp_2d)
c
      i=itap
      CALL gr_fi_ecrit(1,klon,iim,jjm+1,paire,zx_tmp_2d)
      CALL histwrite(physid,"aire",i,zx_tmp_2d)

      iadvtr=iadvtr+1
c
      if (mod(iadvtr,istphy).eq.1.or.istphy.eq.1) then
	print*,'reinitialisation des champs cumules 
     s          a iadvtr=',iadvtr
         do k=1,klev
            do i=1,klon
               mfu(i,k)=0.
               mfd(i,k)=0.
               en_u(i,k)=0.
               de_u(i,k)=0.
               en_d(i,k)=0.
               de_d(i,k)=0.
               coefh(i,k)=0.
                t(i,k)=0.
		fm_therm(i,k)=0.
               entr_therm(i,k)=0.
            enddo
         enddo
         do i=1,klon
            pyv1(i)=0.
            pyu1(i)=0.
         end do
         do k=1,nbsrf
             do i=1,klon
               pftsol(i,k)=0.
               ppsrf(i,k)=0.
            enddo
         enddo

         dtcum=0.
      endif

      do k=1,klev
         do i=1,klon
            mfu(i,k)=mfu(i,k)+pmfu(i,k)*pdtphys
            mfd(i,k)=mfd(i,k)+pmfd(i,k)*pdtphys
            en_u(i,k)=en_u(i,k)+pen_u(i,k)*pdtphys
            de_u(i,k)=de_u(i,k)+pde_u(i,k)*pdtphys
            en_d(i,k)=en_d(i,k)+pen_d(i,k)*pdtphys
            de_d(i,k)=de_d(i,k)+pde_d(i,k)*pdtphys
            coefh(i,k)=coefh(i,k)+pcoefh(i,k)*pdtphys
                t(i,k)=t(i,k)+pt(i,k)*pdtphys
       fm_therm(i,k)=fm_therm(i,k)+pfm_therm(i,k)*pdtphys
       entr_therm(i,k)=entr_therm(i,k)+pentr_therm(i,k)*pdtphys
         enddo
      enddo
         do i=1,klon
            pyv1(i)=pyv1(i)+yv1(i)*pdtphys
            pyu1(i)=pyu1(i)+yu1(i)*pdtphys
         end do
         do k=1,nbsrf
             do i=1,klon
               pftsol(i,k)=pftsol(i,k)+ftsol(i,k)*pdtphys
               ppsrf(i,k)=ppsrf(i,k)+pctsrf(i,k)*pdtphys
            enddo
         enddo

      dtcum=dtcum+pdtphys

      IF(mod(iadvtr,istphy).eq.0) THEN 
c
c   normalisation par le temps cumule
         do k=1,klev
            do i=1,klon
               mfu(i,k)=mfu(i,k)/dtcum
               mfd(i,k)=mfd(i,k)/dtcum
               en_u(i,k)=en_u(i,k)/dtcum
               de_u(i,k)=de_u(i,k)/dtcum
               en_d(i,k)=en_d(i,k)/dtcum
               de_d(i,k)=de_d(i,k)/dtcum
               coefh(i,k)=coefh(i,k)/dtcum
c Unitel a enlever
	      t(i,k)=t(i,k)/dtcum	
               fm_therm(i,k)=fm_therm(i,k)/dtcum
	       entr_therm(i,k)=entr_therm(i,k)/dtcum
            enddo
         enddo
         do i=1,klon
            pyv1(i)=pyv1(i)/dtcum
            pyu1(i)=pyu1(i)/dtcum
         end do
         do k=1,nbsrf
             do i=1,klon
               pftsol(i,k)=pftsol(i,k)/dtcum
               pftsol1(i) = pftsol(i,1)
               pftsol2(i) = pftsol(i,2)
               pftsol3(i) = pftsol(i,3)
               pftsol4(i) = pftsol(i,4)

               ppsrf(i,k)=ppsrf(i,k)/dtcum
               ppsrf1(i) = ppsrf(i,1)
               ppsrf2(i) = ppsrf(i,2)
               ppsrf3(i) = ppsrf(i,3)
               ppsrf4(i) = ppsrf(i,4)

            enddo
         enddo
c
c   ecriture des champs
c
         irec=irec+1

ccccc
         CALL gr_fi_ecrit(klev,klon,iim,jjm+1, t, zx_tmp_3d)
         CALL histwrite(physid,"t",itap,zx_tmp_3d)

         CALL gr_fi_ecrit(klev,klon,iim,jjm+1, mfu, zx_tmp_3d)
      CALL histwrite(physid,"mfu",itap,zx_tmp_3d)
	CALL gr_fi_ecrit(klev,klon,iim,jjm+1, mfd, zx_tmp_3d)
      CALL histwrite(physid,"mfd",itap,zx_tmp_3d)
        CALL gr_fi_ecrit(klev,klon,iim,jjm+1, en_u, zx_tmp_3d)
      CALL histwrite(physid,"en_u",itap,zx_tmp_3d)
        CALL gr_fi_ecrit(klev,klon,iim,jjm+1, de_u, zx_tmp_3d)
      CALL histwrite(physid,"de_u",itap,zx_tmp_3d)
        CALL gr_fi_ecrit(klev,klon,iim,jjm+1, en_d, zx_tmp_3d)
      CALL histwrite(physid,"en_d",itap,zx_tmp_3d)
        CALL gr_fi_ecrit(klev,klon,iim,jjm+1, de_d, zx_tmp_3d)       
      CALL histwrite(physid,"de_d",itap,zx_tmp_3d)
        CALL gr_fi_ecrit(klev,klon,iim,jjm+1, coefh, zx_tmp_3d)         
      CALL histwrite(physid,"coefh",itap,zx_tmp_3d)	

c ajou...
	do k=1,klev
           do i=1,klon
	 fm_therm1(i,k)=fm_therm(i,k)	
	   enddo
	enddo

      CALL gr_fi_ecrit(klev,klon,iim,jjm+1, fm_therm1, zx_tmp_3d)
      CALL histwrite(physid,"fm_th",itap,zx_tmp_3d)
c
      CALL gr_fi_ecrit(klev,klon,iim,jjm+1, entr_therm, zx_tmp_3d)
      CALL histwrite(physid,"en_th",itap,zx_tmp_3d)
cccc
       CALL gr_fi_ecrit(klev,klon,iim,jjm+1,frac_impa,zx_tmp_3d)
        CALL histwrite(physid,"frac_impa",itap,zx_tmp_3d)

        CALL gr_fi_ecrit(klev,klon,iim,jjm+1,frac_nucl,zx_tmp_3d)
        CALL histwrite(physid,"frac_nucl",itap,zx_tmp_3d)
 
        CALL gr_fi_ecrit(1, klon,iim,jjm+1, pyu1,zx_tmp_2d)
      CALL histwrite(physid,"pyu1",itap,zx_tmp_2d)
	
	CALL gr_fi_ecrit(1, klon,iim,jjm+1, pyv1,zx_tmp_2d)
      CALL histwrite(physid,"pyv1",itap,zx_tmp_2d)
	
	CALL gr_fi_ecrit(1,klon,iim,jjm+1, pftsol1, zx_tmp_2d)
      CALL histwrite(physid,"ftsol1",itap,zx_tmp_2d)
         CALL gr_fi_ecrit(1,klon,iim,jjm+1, pftsol2, zx_tmp_2d)
      CALL histwrite(physid,"ftsol2",itap,zx_tmp_2d)
          CALL gr_fi_ecrit(1,klon,iim,jjm+1, pftsol3, zx_tmp_2d)
      CALL histwrite(physid,"ftsol3",itap,zx_tmp_2d)
         CALL gr_fi_ecrit(1,klon,iim,jjm+1, pftsol4, zx_tmp_2d)
      CALL histwrite(physid,"ftsol4",itap,zx_tmp_2d)

        CALL gr_fi_ecrit(1,klon,iim,jjm+1, ppsrf1, zx_tmp_2d)
      CALL histwrite(physid,"psrf1",itap,zx_tmp_2d)
        CALL gr_fi_ecrit(1,klon,iim,jjm+1, ppsrf2, zx_tmp_2d)
      CALL histwrite(physid,"psrf2",itap,zx_tmp_2d)
        CALL gr_fi_ecrit(1,klon,iim,jjm+1, ppsrf3, zx_tmp_2d)
      CALL histwrite(physid,"psrf3",itap,zx_tmp_2d)
        CALL gr_fi_ecrit(1,klon,iim,jjm+1, ppsrf4, zx_tmp_2d)
      CALL histwrite(physid,"psrf4",itap,zx_tmp_2d)

      if (ok_sync) call histsync(physid)
c     if (ok_sync) call histsync
	
c
cAA Test sur la valeur des coefficients de lessivage 
c
         zmin=1e33
         zmax=-1e33
         do k=1,klev
            do i=1,klon
                  zmax=max(zmax,frac_nucl(i,k))
                  zmin=min(zmin,frac_nucl(i,k))
            enddo
         enddo
         Print*,'------ coefs de lessivage (min et max) --------'
         Print*,'facteur de nucleation ',zmin,zmax
         zmin=1e33
         zmax=-1e33
         do k=1,klev
            do i=1,klon
                  zmax=max(zmax,frac_impa(i,k))
                  zmin=min(zmin,frac_impa(i,k))
            enddo
         enddo
         Print*,'facteur d impaction ',zmin,zmax

      ENDIF 

c   reinitialisation des champs cumules
	go to 768
      if (mod(iadvtr,istphy).eq.1) then
         do k=1,klev
            do i=1,klon
               mfu(i,k)=0.
               mfd(i,k)=0.
               en_u(i,k)=0.
               de_u(i,k)=0.
               en_d(i,k)=0.
               de_d(i,k)=0.
               coefh(i,k)=0.
	       t(i,k)=0.
               fm_therm(i,k)=0.
	       entr_therm(i,k)=0.
            enddo
         enddo
         do i=1,klon
            pyv1(i)=0.
            pyu1(i)=0.
         end do
         do k=1,nbsrf
             do i=1,klon
               pftsol(i,k)=0.
               ppsrf(i,k)=0.
            enddo
         enddo

         dtcum=0.
      endif

      do k=1,klev
         do i=1,klon
            mfu(i,k)=mfu(i,k)+pmfu(i,k)*pdtphys
            mfd(i,k)=mfd(i,k)+pmfd(i,k)*pdtphys
            en_u(i,k)=en_u(i,k)+pen_u(i,k)*pdtphys
            de_u(i,k)=de_u(i,k)+pde_u(i,k)*pdtphys
            en_d(i,k)=en_d(i,k)+pen_d(i,k)*pdtphys
            de_d(i,k)=de_d(i,k)+pde_d(i,k)*pdtphys
            coefh(i,k)=coefh(i,k)+pcoefh(i,k)*pdtphys
		t(i,k)=t(i,k)+pt(i,k)*pdtphys
       fm_therm(i,k)=fm_therm(i,k)+pfm_therm(i,k)*pdtphys
       entr_therm(i,k)=entr_therm(i,k)+pentr_therm(i,k)*pdtphys
         enddo
      enddo
         do i=1,klon
            pyv1(i)=pyv1(i)+yv1(i)*pdtphys
            pyu1(i)=pyu1(i)+yu1(i)*pdtphys
         end do
         do k=1,nbsrf
             do i=1,klon
               pftsol(i,k)=pftsol(i,k)+ftsol(i,k)*pdtphys
               ppsrf(i,k)=ppsrf(i,k)+pctsrf(i,k)*pdtphys
            enddo
         enddo

      dtcum=dtcum+pdtphys
768   continue

      RETURN
      END
