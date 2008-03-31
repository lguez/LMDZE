!
! $Header: /home/cvsroot/LMDZ4/libf/bibio/initfluxsto.F,v 1.1.1.1 2004/05/19 12:53:05 lmdzadmin Exp $
!
      subroutine initfluxsto
     .  (infile,tstep,t_ops,t_wrt,nq,
     .                    fileid,filevid,filedid)

       USE IOIPSL

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
       use dimens_m
      use paramet_m
      use comconst
      use comvert
      use logic
      use comgeom
      use serre
      use temps, only: annee_ref, day_ref, itau_dyn
      use ener
      implicit none

C   Arguments
C
      character*(*) infile
      integer*4 itau
      real tstep, t_ops, t_wrt
      integer fileid, filevid,filedid
      integer nq,ndex(1)
      real nivd(1)

C   Variables locales
C
      integer tau0
      real zjulian
      character*3 str
      character*10 ctrac
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1),rl(1,1)
      integer uhoriid, vhoriid, thoriid, zvertiid,dhoriid,dvertiid
      integer ii,jj
      integer zan, idayref
      logical ok_sync
C
C  Initialisations
C
      pi = 4. * atan (1.)
      str='q  '
      ctrac = 'traceur   '
      ok_sync = .true.
C
C  Appel a histbeg: creation du fichier netcdf et initialisations diverses
C         

      zan = annee_ref
      idayref = day_ref
      CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)
      tau0 = itau_dyn
	
	do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonu(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo
 
      call histbeg_totreg(infile, iip1, rlong(:,1), jjp1, rlat(1,:),
     .             1, iip1, 1, jjp1,
     .             tau0, zjulian, tstep, uhoriid, fileid)
C
C  Creation du fichier histoire pour la grille en V (oblige pour l'instant,
C  IOIPSL ne permet pas de grilles avec des nombres de point differents dans 
C  un meme fichier)


      do jj = 1, jjm
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatv(jj) * 180. / pi
        enddo
      enddo

      call histbeg_totreg('fluxstokev.nc', iip1, rlong(:,1), jjm,
     .             rlat(1,:),1, iip1, 1, jjm,
     .             tau0, zjulian, tstep, vhoriid, filevid)
	
	rl(1,1) = 1.	
      call histbeg_regular('defstoke.nc', 1, rl, 1, rl,
     .             1, 1, 1, 1,
     .             tau0, zjulian, tstep, dhoriid, filedid)

C
C  Appel a histhori pour rajouter les autres grilles horizontales
C
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo

      call histhori(fileid, iip1, rlong, jjp1, rlat, 'scalar',
     .              'Grille points scalaires', thoriid)
	
C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'sig_s', 'Niveaux sigma',
     . 'sigma_level',
     .              llm, nivsigs, zvertiid)
C Pour le fichier V
      call histvert(filevid, 'sig_s', 'Niveaux sigma',
     .  'sigma_level',
     .              llm, nivsigs, zvertiid)
c pour le fichier def
      nivd(1) = 1
      call histvert(filedid, 'sig_s', 'Niveaux sigma',
     .  'sigma_level',
     .              1, nivd, dvertiid)

C
C  Appels a histdef pour la definition des variables a sauvegarder
	
	CALL histdef(fileid, "phis", "Surface geop. height", "-",
     .                iip1,jjp1,thoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)

         CALL histdef(fileid, "aire", "Grid area", "-",
     .                iip1,jjp1,thoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
	
	CALL histdef(filedid, "dtvr", "tps dyn", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
        
         CALL histdef(filedid, "istdyn", "tps stock", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
         
         CALL histdef(filedid, "istphy", "tps stock phy", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)


C
C Masse 
C
      call histdef(fileid, 'masse', 'Masse', 'kg',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Pbaru 
C
      call histdef(fileid, 'pbaru', 'flx de masse zonal', 'kg m/s',
     .             iip1, jjp1, uhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Pbarv 
C
      call histdef(filevid, 'pbarv', 'flx de masse mer', 'kg m/s',
     .             iip1, jjm, vhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  w 
C
      call histdef(fileid, 'w', 'flx de masse vert', 'kg m/s',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Temperature potentielle
C
      call histdef(fileid, 'teta', 'temperature potentielle', '-',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C

C
C Geopotentiel 
C
      call histdef(fileid, 'phi', 'geopotentiel instantane', '-',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Fin
C
      call histend(fileid)
      call histend(filevid)
      call histend(filedid)
      if (ok_sync) then
        call histsync(fileid)
        call histsync(filevid)
        call histsync(filedid)
      endif
	
      return
      end
