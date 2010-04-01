!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/initphysto.F,v 1.2 2004/06/22 11:45:32 lmdzadmin Exp $
!
C
C
      subroutine initphysto
     .  (infile,
     .  rlon, rlat, tstep,t_ops,t_wrt,nq,fileid)

      USE calendar
      use histcom
      use dimens_m
      use paramet_m
      use comconst
      use indicesol
      use dimphy
      use logic
      use comgeom
      use serre
      use temps
      use ener
      implicit none

C
C   Routine d'initialisation des ecritures des fichiers histoires LMDZ
C   au format IOIPSL
C
C   Appels succesifs des routines: histbeg
C                                  histhori
C                                  histver
C                                  histdef
C                                  histend
C
C   Entree:
C
C      infile: nom du fichier histoire a creer
C      day0,anne0: date de reference
C      tstep: duree du pas de temps en seconde
C      t_ops: frequence de l'operation pour IOIPSL
C      t_wrt: frequence d'ecriture sur le fichier
C      nq: nombre de traceurs
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
C      filevid:ID du fichier netcdf pour la grille v
C
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations

C   Arguments
      character*(*) infile
      integer*4 nhoriid, i
      real, intent(in):: tstep
      real t_ops, t_wrt
      integer fileid, filevid
      integer nq,l
      real nivsigs(llm)

C   Variables locales
C
      integer tau0
      real zjulian
      character*3 str
      character*10 ctrac
      integer iq
      integer uhoriid, vhoriid, thoriid, zvertiid
      integer ii,jj
      integer zan, idayref
      logical ok_sync
      REAL zx_lon(iim,jjm+1), zx_lat(iim,jjm+1)
C
      REAL, intent(in):: rlon(klon), rlat(klon)

C  Initialisations
C
      pi = 4. * atan (1.)
      str='q  '
      ctrac = 'traceur   '
      ok_sync= .true.
C
C  Appel a histbeg: creation du fichier netcdf et initialisations diverses
C         

      zan = annee_ref
      idayref = day_ref
      CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)
      tau0 = 0
	
	CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,zx_lon)
         DO i = 1, iim
            zx_lon(i,1) = rlon(i+1)
            zx_lon(i,jjm+1) = rlon(i+1)
         ENDDO
         CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,zx_lat)


      call histbeg_totreg(infile, zx_lon(:,1), zx_lat(1,:),
     .             1, iim, 1, jjm+1,
     .             tau0, zjulian, tstep, nhoriid, fileid)

C  Appel a histvert pour la grille verticale
C
	DO l=1,llm
            nivsigs(l)=float(l)
         ENDDO
	
	write(*,*) 'avant histvert ds initphysto'

      call histvert(fileid, 'sig_s', 'Niveaux sigma',
     . 'sigma_level',
     .              llm, nivsigs, zvertiid)
C
C  Appels a histdef pour la definition des variables a sauvegarder
C
	write(*,*) 'apres histvert ds initphysto'

       CALL histdef(fileid, "phis", "Surface geop. height", "-",
     .                iim,jjm+1,nhoriid, 1,1,1, -99, 
     .                "once", t_ops, t_wrt)
c
	write(*,*) 'apres phis ds initphysto'

         CALL histdef(fileid, "aire", "Grid area", "-",
     .                iim,jjm+1,nhoriid, 1,1,1, -99, 
     .                "once", t_ops, t_wrt)
         write(*,*) 'apres aire ds initphysto'

         CALL histdef(fileid, "dtime", "tps phys ", "s",
     .                1,1,nhoriid, 1,1,1, -99, 
     .                "once", t_ops, t_wrt)
	
	 CALL histdef(fileid, "istphy", "tps stock", "s",
     .                1,1,nhoriid, 1,1,1, -99, 
     .                "once", t_ops, t_wrt)

C T 
C
      call histdef(fileid, 't', 'Temperature', 'K',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)
        write(*,*) 'apres t ds initphysto'
C mfu 
C
      call histdef(fileid, 'mfu', 'flx m. pan. mt', 'kg m/s',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)
        write(*,*) 'apres mfu ds initphysto'
C
C mfd 
C
      call histdef(fileid, 'mfd', 'flx m. pan. des', 'kg m/s',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)

C
C en_u 
C
      call histdef(fileid, 'en_u', 'flx ent pan mt', 'kg m/s',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)
               write(*,*) 'apres en_u ds initphysto'
C
C de_u 
C
      call histdef(fileid, 'de_u', 'flx det pan mt', 'kg m/s',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)

C
C en_d 
C
      call histdef(fileid, 'en_d', 'flx ent pan dt', 'kg m/s',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)
C

C
C de_d 
C
      call histdef(fileid, 'de_d', 'flx det pan dt', 'kg m/s',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)

c coefh frac_impa,frac_nucl
	
	call histdef(fileid, "coefh", " ", " ",
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              "inst(X)", t_ops, t_wrt)

c abderrahmane le 16 09 02
        call histdef(fileid, "fm_th", " ", " ",
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              "inst(X)", t_ops, t_wrt)

        call histdef(fileid, "en_th", " ", " ",
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              "inst(X)", t_ops, t_wrt)
c fin aj
	
	write(*,*) 'apres coefh ds initphysto'	

	call histdef(fileid, 'frac_impa', ' ', ' ',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)
	
	call histdef(fileid, 'frac_nucl', ' ', ' ',
     .             iim, jjm+1, nhoriid, llm, 1, llm, zvertiid,
     .              'inst(X)', t_ops, t_wrt)

c
c pyu1
c
      CALL histdef(fileid, "pyu1", " ", " ",
     .                iim,jjm+1,nhoriid, 1,1,1, -99, 
     .                "inst(X)", t_ops, t_wrt)

c
c pyv1
c
	CALL histdef(fileid, "pyv1", " ", " ",
     .                iim,jjm+1,nhoriid, 1,1,1, -99, 
     .                "inst(X)", t_ops, t_wrt)
	
	write(*,*) 'apres pyv1 ds initphysto'
c
c ftsol1
c
	call histdef(fileid, "ftsol1", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1,1, -99,
     .             "inst(X)", t_ops, t_wrt)

c
c ftsol2
c
        call histdef(fileid, "ftsol2", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1,1, -99,
     .             "inst(X)", t_ops, t_wrt)

c
c ftsol3
c
        call histdef(fileid, "ftsol3", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1,1, -99,
     .              "inst(X)", t_ops, t_wrt)

c
c ftsol4
c
        call histdef(fileid, "ftsol4", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1,1, -99,
     .              "inst(X)", t_ops, t_wrt)
	
c
c rain
c
        call histdef(fileid, "rain", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1,1, -99,
     .              "inst(X)", t_ops, t_wrt)

c
c psrf1
c
	call histdef(fileid, "psrf1", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1, 1, -99,
     .              "inst(X)", t_ops, t_wrt)
	
c
c psrf2
c
        call histdef(fileid, "psrf2", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1, 1, -99,
     .              "inst(X)", t_ops, t_wrt)

c
c psrf3
c
        call histdef(fileid, "psrf3", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1, 1, -99,
     .              "inst(X)", t_ops, t_wrt)

c
c psrf4
c
        call histdef(fileid, "psrf4", " ", " ",
     .             iim, jjm+1, nhoriid, 1, 1, 1, -99,
     .              "inst(X)", t_ops, t_wrt)
	
	write(*,*) 'avant histend ds initphysto'	

      call histend(fileid)
c     if (ok_sync) call histsync(fileid)
      if (ok_sync) call histsync

	

      return
      end
