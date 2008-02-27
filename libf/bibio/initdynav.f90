module initdynav_m

  ! This module is clean: no C preprocessor directive, no include line

  implicit none

contains

  subroutine initdynav(day0, anne0, tstep, nq, fileid, infile, t_ops, t_wrt)

    ! From initdynav.F,v 1.1.1.1 2004/05/19 12:53:05

    USE IOIPSL, only: ymds2ju, histbeg_totreg, histvert, histdef, histend
    use dimens_m
    use paramet_m
    use comconst, only: pi
    use comvert, only: nivsigs
    use logic
    use comgeom
    use serre
    use temps
    use ener
    use advtrac_m, only: ttext

    !   Routine d'initialisation des ecritures des fichiers histoires LMDZ
    !   au format IOIPSL. Initialisation du fichier histoire moyenne.

    !   Appels succesifs des routines: histbeg
    !                                  histhori
    !                                  histver
    !                                  histdef
    !                                  histend

    !   Entree:
    !      infile: nom du fichier histoire a creer
    !      day0,anne0: date de reference
    !      tstep : frequence d'ecriture
    !      t_ops: frequence de l'operation pour IOIPSL
    !      t_wrt: frequence d'ecriture sur le fichier
    !      nq: nombre de traceurs

    !   Sortie:
    !      fileid: ID du fichier netcdf cree

    !   L. Fairhead, LMD, 03/99

    !   Arguments
    character(len=*) infile
    integer day0, anne0
    real, intent(in):: tstep, t_ops, t_wrt
    integer fileid
    integer nq
    integer thoriid, zvertiid

    !   Variables locales

    integer tau0
    real zjulian
    integer iq
    real rlong(iip1,jjp1), rlat(iip1,jjp1)
    integer ii,jj
    integer zan, dayref

    !----------------------------------------------------

    !  Initialisations

    pi = 4. * atan (1.)
    !
    !  Appel a histbeg: creation du fichier netcdf et initialisations diverses
    !         

    zan = anne0
    dayref = day0
    CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
    tau0 = itau_dyn

    do jj = 1, jjp1
       do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj)  = rlatu(jj) * 180. / pi
       enddo
    enddo

    call histbeg_totreg(infile, iip1, rlong(:,1), jjp1, rlat(1,:), &
         1, iip1, 1, jjp1, &
         tau0, zjulian, tstep, thoriid, fileid)

    !
    !  Appel a histvert pour la grille verticale
    !
    call histvert(fileid, 'sigss', 'Niveaux sigma','Pa', &
         llm, nivsigs, zvertiid)
    !
    !  Appels a histdef pour la definition des variables a sauvegarder
    !
    !  Vents U
    !
    write(6,*)'inithistave',tstep
    call histdef(fileid, 'u', 'vents u scalaires moyennes', &
         'm/s', iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         32, 'ave(X)', t_ops, t_wrt)

    !
    !  Vents V
    !
    call histdef(fileid, 'v', 'vents v scalaires moyennes', &
         'm/s', iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         32, 'ave(X)', t_ops, t_wrt)

    !
    !  Temperature
    !
    call histdef(fileid, 'temp', 'temperature moyennee', 'K', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         32, 'ave(X)', t_ops, t_wrt)
    !
    !  Temperature potentielle
    !
    call histdef(fileid, 'theta', 'temperature potentielle', 'K', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         32, 'ave(X)', t_ops, t_wrt)


    !
    !  Geopotentiel
    !
    call histdef(fileid, 'phi', 'geopotentiel moyenne', '-', &
         iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
         32, 'ave(X)', t_ops, t_wrt)
    !
    !  Traceurs
    !
    DO iq=1,nq
       call histdef(fileid, ttext(iq), ttext(iq), '-', &
            iip1, jjp1, thoriid, llm, 1, llm, zvertiid, &
            32, 'ave(X)', t_ops, t_wrt)
    enddo
    !
    !  Masse
    !
    call histdef(fileid, 'masse', 'masse', 'kg', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         32, 'ave(X)', t_ops, t_wrt)
    !
    !  Pression au sol
    !
    call histdef(fileid, 'ps', 'pression naturelle au sol', 'Pa', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         32, 'ave(X)', t_ops, t_wrt)
    !
    !  Pression au sol
    !
    call histdef(fileid, 'phis', 'geopotentiel au sol', '-', &
         iip1, jjp1, thoriid, 1, 1, 1, -99, &
         32, 'ave(X)', t_ops, t_wrt)
    !
    !  Fin
    !
    call histend(fileid)

  end subroutine initdynav

end module initdynav_m
