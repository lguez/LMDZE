SUBROUTINE initphysto(infile,rlon,rlat,tstep,t_ops,t_wrt,nq,fileid)

  ! From phylmd/initphysto.F,v 1.2 2004/06/22 11:45:32
  !   Routine d'initialisation des ecritures des fichiers histoires LMDZ
  !   au format IOIPSL

  !   Appels succesifs des routines: histbeg
  !                                  histhori
  !                                  histver
  !                                  histdef
  !                                  histend

  !   Entree:

  !      infile: nom du fichier histoire a creer
  !      day0,anne0: date de reference
  !      tstep: duree du pas de temps en seconde
  !      t_ops: frequence de l'operation pour IOIPSL
  !      t_wrt: frequence d'ecriture sur le fichier
  !      nq: nombre de traceurs

  !   Sortie:
  !      fileid: ID du fichier netcdf cree
  !      filevid:ID du fichier netcdf pour la grille v

  !   L. Fairhead, 03/99

  USE calendar
  USE histbeg_totreg_m, ONLY : histbeg_totreg
  USE histdef_m, ONLY : histdef
  USE histend_m, ONLY : histend
  use histsync_m, only: histsync
  USE histvert_m, ONLY : histvert
  USE dimens_m
  USE paramet_m
  USE comconst
  USE indicesol
  USE dimphy
  use conf_gcm_m
  USE comgeom
  USE serre
  USE temps
  USE ener
  USE nr_util, ONLY : pi

  IMPLICIT NONE

  !   Arguments
  CHARACTER*(*) infile
  INTEGER nhoriid, i
  REAL, INTENT (IN) :: tstep
  REAL t_ops, t_wrt
  INTEGER fileid, filevid
  INTEGER nq, l
  REAL nivsigs(llm)

  !   Variables locales

  INTEGER tau0
  REAL zjulian
  CHARACTER*3 str
  CHARACTER*10 ctrac
  INTEGER iq
  INTEGER uhoriid, vhoriid, thoriid, zvertiid
  INTEGER ii, jj
  INTEGER zan, idayref
  LOGICAL ok_sync
  REAL zx_lon(iim,jjm+1), zx_lat(iim,jjm+1)

  REAL, INTENT (IN) :: rlon(klon), rlat(klon)

  !-----------------------------------------------------

  !  Initialisations
  str = 'q  '
  ctrac = 'traceur   '
  ok_sync = .TRUE.

  !  Appel a histbeg: creation du fichier netcdf et initialisations
  !     diverses

  zan = annee_ref
  idayref = day_ref
  CALL ymds2ju(zan,1,idayref,0.0,zjulian)
  tau0 = 0

  CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,zx_lon)
  DO i = 1, iim
     zx_lon(i,1) = rlon(i+1)
     zx_lon(i,jjm+1) = rlon(i+1)
  END DO
  CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,zx_lat)


  CALL histbeg_totreg(infile,zx_lon(:,1),zx_lat(1,:),1,iim,1,jjm+1,tau0, &
       zjulian,tstep,nhoriid,fileid)

  !  Appel a histvert pour la grille verticale

  DO l = 1, llm
     nivsigs(l) = float(l)
  END DO

  CALL histvert(fileid,'sig_s','Niveaux sigma','sigma_level',llm,nivsigs, &
       zvertiid)

  !  Appels a histdef pour la definition des variables a sauvegarder

  CALL histdef(fileid,'phis','Surface geop. height','-',iim,jjm+1,nhoriid, &
       1,1,1,-99,'once',t_ops,t_wrt)

  CALL histdef(fileid,'aire','Grid area','-',iim,jjm+1,nhoriid,1,1,1,-99, &
       'once',t_ops,t_wrt)

  CALL histdef(fileid,'dtime','tps phys ','s',1,1,nhoriid,1,1,1,-99, &
       'once',t_ops,t_wrt)

  CALL histdef(fileid,'istphy','tps stock','s',1,1,nhoriid,1,1,1,-99, &
       'once',t_ops,t_wrt)

  ! T

  CALL histdef(fileid,'t','Temperature','K',iim,jjm+1,nhoriid,llm,1,llm, &
       zvertiid,'inst(X)',t_ops,t_wrt)

  CALL histdef(fileid,'mfu','flx m. pan. mt','kg m/s',iim,jjm+1,nhoriid, &
       llm,1,llm,zvertiid,'inst(X)',t_ops,t_wrt)

  CALL histdef(fileid,'mfd','flx m. pan. des','kg m/s',iim,jjm+1,nhoriid, &
       llm,1,llm,zvertiid,'inst(X)',t_ops,t_wrt)


  ! en_u

  CALL histdef(fileid,'en_u','flx ent pan mt','kg m/s',iim,jjm+1,nhoriid, &
       llm,1,llm,zvertiid,'inst(X)',t_ops,t_wrt)

  CALL histdef(fileid,'de_u','flx det pan mt','kg m/s',iim,jjm+1,nhoriid, &
       llm,1,llm,zvertiid,'inst(X)',t_ops,t_wrt)


  ! en_d

  CALL histdef(fileid,'en_d','flx ent pan dt','kg m/s',iim,jjm+1,nhoriid, &
       llm,1,llm,zvertiid,'inst(X)',t_ops,t_wrt)



  ! de_d

  CALL histdef(fileid,'de_d','flx det pan dt','kg m/s',iim,jjm+1,nhoriid, &
       llm,1,llm,zvertiid,'inst(X)',t_ops,t_wrt)

  ! coefh frac_impa,frac_nucl

  CALL histdef(fileid,'coefh',' ',' ',iim,jjm+1,nhoriid,llm,1,llm, &
       zvertiid,'inst(X)',t_ops,t_wrt)

  ! abderrahmane le 16 09 02
  CALL histdef(fileid,'fm_th',' ',' ',iim,jjm+1,nhoriid,llm,1,llm, &
       zvertiid,'inst(X)',t_ops,t_wrt)

  CALL histdef(fileid,'en_th',' ',' ',iim,jjm+1,nhoriid,llm,1,llm, &
       zvertiid,'inst(X)',t_ops,t_wrt)
  ! fin aj

  CALL histdef(fileid,'frac_impa',' ',' ',iim,jjm+1,nhoriid,llm,1,llm, &
       zvertiid,'inst(X)',t_ops,t_wrt)

  CALL histdef(fileid,'frac_nucl',' ',' ',iim,jjm+1,nhoriid,llm,1,llm, &
       zvertiid,'inst(X)',t_ops,t_wrt)


  ! pyu1

  CALL histdef(fileid,'pyu1',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! pyv1

  CALL histdef(fileid,'pyv1',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)

  CALL histdef(fileid,'ftsol1',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! ftsol2

  CALL histdef(fileid,'ftsol2',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! ftsol3

  CALL histdef(fileid,'ftsol3',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! ftsol4

  CALL histdef(fileid,'ftsol4',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! rain

  CALL histdef(fileid,'rain',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! psrf1

  CALL histdef(fileid,'psrf1',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! psrf2

  CALL histdef(fileid,'psrf2',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! psrf3

  CALL histdef(fileid,'psrf3',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)


  ! psrf4

  CALL histdef(fileid,'psrf4',' ',' ',iim,jjm+1,nhoriid,1,1,1,-99, &
       'inst(X)',t_ops,t_wrt)

  CALL histend(fileid)
  IF (ok_sync) CALL histsync

END SUBROUTINE initphysto
