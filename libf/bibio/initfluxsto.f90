SUBROUTINE initfluxsto(tstep, t_ops, t_wrt, nq, fileid, filevid, filedid)

  ! From bibio/initfluxsto.F, v 1.1.1.1 2004/05/19 12:53:05

  !   Routine d'initialisation des ecritures des fichiers histoires LMDZ  
  !   au format IOIPSL                                                    
  !   Appels succesifs des routines: histbeg                              
  !                                  histhori                             
  !                                  histver                              
  !                                  histdef                              
  !                                  histend                              

  !   Entree:                                                             
  !      day0, anne0: date de reference                                    
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
  use histhori_regular_m, only: histhori_regular
  use histsync_m, only: histsync
  USE histvert_m, ONLY : histvert
  USE dimens_m
  USE paramet_m
  USE comconst
  USE disvert_m
  use conf_gcm_m
  USE comgeom
  USE serre
  USE temps, ONLY : annee_ref, day_ref, itau_dyn
  USE ener
  USE nr_util, ONLY : pi

  IMPLICIT NONE

  !   Arguments                                                           
  INTEGER itau
  REAL, INTENT (IN) :: tstep
  REAL t_ops, t_wrt
  INTEGER fileid, filevid, filedid
  INTEGER nq, ndex(1)
  REAL nivd(1)

  !   Variables locales                                                   
  REAL zjulian
  CHARACTER*3 str
  CHARACTER*10 ctrac
  INTEGER iq
  REAL rlong(iip1, jjp1), rlat(iip1, jjp1)
  INTEGER uhoriid, vhoriid, thoriid, zvertiid, dhoriid, dvertiid
  INTEGER ii, jj
  INTEGER zan, idayref
  LOGICAL ok_sync

  !---------------------------------------------------------

  !  Initialisations                                                      
  str = 'q  '
  ctrac = 'traceur   '
  ok_sync = .TRUE.

  !  Appel a histbeg: creation du fichier netcdf et initialisations diverses

  zan = annee_ref
  idayref = day_ref
  CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)

  DO jj = 1, jjp1
     DO ii = 1, iip1
        rlong(ii, jj) = rlonu(ii)*180./pi
        rlat(ii, jj) = rlatu(jj)*180./pi
     END DO
  END DO

  CALL histbeg_totreg('fluxstoke', rlong(:, 1), rlat(1, :), 1, iip1, 1, jjp1, &
       itau_dyn, zjulian, tstep, uhoriid, fileid)

  !  Creation du fichier histoire pour la grille en V (oblige pour l'instant,
  !  IOIPSL ne permet pas de grilles avec des nombres de point differents dans
  !  un meme fichier)                                                     

  DO jj = 1, jjm
     DO ii = 1, iip1
        rlong(ii, jj) = rlonv(ii)*180./pi
        rlat(ii, jj) = rlatv(jj)*180./pi
     END DO
  END DO

  CALL histbeg_totreg('fluxstokev.nc', rlong(:, 1), rlat(1, :jjm), 1, iip1, &
       1, jjm, itau_dyn, zjulian, tstep, vhoriid, filevid)

  CALL histbeg_totreg('defstoke.nc', (/1./), (/1./), 1, 1, 1, 1, itau_dyn, &
       zjulian, tstep, dhoriid, filedid)

  !  Appel a histhori pour rajouter les autres grilles horizontales       

  DO jj = 1, jjp1
     DO ii = 1, iip1
        rlong(ii, jj) = rlonv(ii)*180./pi
        rlat(ii, jj) = rlatu(jj)*180./pi
     END DO
  END DO

  CALL histhori_regular(fileid, iip1, rlong, jjp1, rlat, 'scalar', &
       'Grille points scalaires', thoriid)

  !  Appel a histvert pour la grille verticale                            

  CALL histvert(fileid, 'sig_s', 'Niveaux sigma', 'sigma_level', llm, &
       nivsigs, zvertiid)
  ! Pour le fichier V                                                     
  CALL histvert(filevid, 'sig_s', 'Niveaux sigma', 'sigma_level', llm, &
       nivsigs, zvertiid)
  ! pour le fichier def                                                   
  nivd(1) = 1
  CALL histvert(filedid, 'sig_s', 'Niveaux sigma', 'sigma_level', 1, nivd, &
       dvertiid)

  !  Appels a histdef pour la definition des variables a sauvegarder      
  CALL histdef(fileid, 'phis', 'Surface geop. height', '-', iip1, jjp1, &
       thoriid, 1, 1, 1, -99, 'once', t_ops, t_wrt)
  CALL histdef(fileid, 'aire', 'Grid area', '-', iip1, jjp1, thoriid, 1, 1, &
       1, -99, 'once', t_ops, t_wrt)
  CALL histdef(filedid, 'dtvr', 'tps dyn', 's', 1, 1, dhoriid, 1, 1, 1, -99, &
       'once', t_ops, t_wrt)
  CALL histdef(filedid, 'istdyn', 'tps stock', 's', 1, 1, dhoriid, 1, 1, 1, &
       -99, 'once', t_ops, t_wrt)
  CALL histdef(filedid, 'istphy', 'tps stock phy', 's', 1, 1, dhoriid, 1, 1, &
       1, -99, 'once', t_ops, t_wrt)
  CALL histdef(fileid, 'masse', 'Masse', 'kg', iip1, jjp1, thoriid, llm, 1, &
       llm, zvertiid, 'inst(X)', t_ops, t_wrt)
  CALL histdef(fileid, 'pbaru', 'flx de masse zonal', 'kg m/s', iip1, jjp1, &
       uhoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
  CALL histdef(filevid, 'pbarv', 'flx de masse mer', 'kg m/s', iip1, jjm, &
       vhoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
  CALL histdef(fileid, 'w', 'flx de masse vert', 'kg m/s', iip1, jjp1, &
       thoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
  CALL histdef(fileid, 'teta', 'temperature potentielle', '-', iip1, jjp1, &
       thoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)
  CALL histdef(fileid, 'phi', 'geopotentiel instantane', '-', iip1, jjp1, &
       thoriid, llm, 1, llm, zvertiid, 'inst(X)', t_ops, t_wrt)

  CALL histend(fileid)
  CALL histend(filevid)
  CALL histend(filedid)
  IF (ok_sync) THEN
     CALL histsync(fileid)
     CALL histsync(filevid)
     CALL histsync(filedid)
  END IF

END SUBROUTINE initfluxsto
