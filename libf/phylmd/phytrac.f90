module phytrac_m

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT none

  private
  public phytrac

contains

  SUBROUTINE phytrac(rnpb, nstep, julien, gmtime, debutphy, lafin, nqmax, &
       pdtphys, u, v, t_seri, paprs, pplay, pmfu, pmfd, pen_u, &
       pde_u, pen_d, pde_d, coefh, fm_therm, entr_therm, yu1, yv1, ftsol, &
       pctsrf, frac_impa, frac_nucl, presnivs, pphis, &
       pphi, albsol, sh, rh, cldfra, rneb, diafra, cldliq, itop_con, &
       ibas_con, pmflxr, pmflxs, prfl, psfl, da, phi, mp, upwd, dnwd, tr_seri)

    ! From phylmd/phytrac.F, version 1.15 2006/02/21 08:08:30

    ! Authors : Frédéric Hourdin, Abderrahmane Idelkadi, Marie-Alice
    ! Foujols, Olivia
    ! Objet : moniteur général des tendances des traceurs

    ! Remarques :
    ! 1/ L'appel de "phytrac" se fait avec "nq-2" donc nous avons bien 
    ! les vrais traceurs (en nombre "nbtr", sans la vapeur d'eau ni l'eau
    ! liquide) dans "phytrac".
    ! 2/ Le choix du radon et du plomb se fait juste avec un "data" 
    ! (peu propre).
    ! Pourrait-on avoir une variable qui indiquerait le type de traceur ?

    use dimens_m, only: iim, jjm, llm
    use indicesol, only: nbsrf
    use dimphy, only: klon, nbtr
    use clesphys, only: ecrit_tra, iflag_con
    use abort_gcm_m, only: abort_gcm
    use YOMCST, only: rg
    use ctherm, only: iflag_thermals
    use read_coefoz_m, only: read_coefoz
    use phyetat0_m, only: rlat
    use o3_chem_m, only: o3_chem

    ! Arguments:

    !   EN ENTREE:

    !   divers:

    logical, intent(in):: rnpb

    integer, intent(in):: nqmax
    ! (nombre de traceurs auxquels on applique la physique)

    integer, intent(in):: nstep  ! appel physique
    integer, intent(in):: julien !jour julien, 1 <= julien <= 360
    integer itop_con(klon)
    integer ibas_con(klon)
    real, intent(in):: gmtime ! heure de la journée en fraction de jour
    real pdtphys  ! pas d'integration pour la physique (s)
    real, intent(in):: t_seri(klon, llm) ! temperature, in K

    real tr_seri(klon, llm, nbtr)
    ! (mass fractions of tracers, excluding water, at mid-layers)

    real u(klon, llm)
    real v(klon, llm)
    real sh(klon, llm)     ! humidite specifique
    real rh(klon, llm)     ! humidite relative
    real cldliq(klon, llm) ! eau liquide nuageuse
    real cldfra(klon, llm) ! fraction nuageuse (tous les nuages)

    real diafra(klon, llm)
    ! (fraction nuageuse (convection ou stratus artificiels))

    real rneb(klon, llm)   ! fraction nuageuse (grande echelle)
    real albsol(klon)  ! albedo surface

    real, intent(in):: paprs(klon, llm+1)
    ! (pression pour chaque inter-couche, en Pa)

    real pplay(klon, llm)  ! pression pour le mileu de chaque couche (en Pa)
    real pphi(klon, llm) ! geopotentiel
    real pphis(klon)
    REAL, intent(in):: presnivs(llm)
    logical, intent(in):: debutphy ! le flag de l'initialisation de la physique
    logical, intent(in):: lafin ! fin de la physique

    integer nsplit
    REAL pmflxr(klon, llm+1), pmflxs(klon, llm+1)   !--lessivage convection
    REAL prfl(klon, llm+1),   psfl(klon, llm+1)     !--lessivage large-scale

    !   convection:

    REAL pmfu(klon, llm)  ! flux de masse dans le panache montant
    REAL pmfd(klon, llm)  ! flux de masse dans le panache descendant
    REAL pen_u(klon, llm) ! flux entraine dans le panache montant

    !   thermiques:

    real fm_therm(klon, llm+1), entr_therm(klon, llm)

    REAL pde_u(klon, llm) ! flux detraine dans le panache montant
    REAL pen_d(klon, llm) ! flux entraine dans le panache descendant
    REAL pde_d(klon, llm) ! flux detraine dans le panache descendant
    ! KE
    real da(klon, llm), phi(klon, llm, llm), mp(klon, llm)
    REAL upwd(klon, llm)      ! saturated updraft mass flux
    REAL dnwd(klon, llm)      ! saturated downdraft mass flux

    !   Couche limite:

    REAL coefh(klon, llm) ! coeff melange CL
    REAL yu1(klon)        ! vents au premier niveau
    REAL yv1(klon)        ! vents au premier niveau

    !   Lessivage:

    ! pour le ON-LINE

    REAL frac_impa(klon, llm)  ! fraction d'aerosols impactes
    REAL frac_nucl(klon, llm)  ! fraction d'aerosols nuclees

    ! Arguments necessaires pour les sources et puits de traceur:

    real ftsol(klon, nbsrf)  ! Temperature du sol (surf)(Kelvin)
    real pctsrf(klon, nbsrf) ! Pourcentage de sol f(nature du sol)

    real pftsol1(klon), pftsol2(klon), pftsol3(klon), pftsol4(klon)
    real ppsrf1(klon), ppsrf2(klon), ppsrf3(klon), ppsrf4(klon)

    !  VARIABLES LOCALES TRACEURS

    ! Sources et puits des traceurs:

    ! Pour l'instant seuls les cas du rn et du pb ont ete envisages.

    REAL source(klon)       ! a voir lorsque le flux est prescrit 
    ! 
    ! Pour la source de radon et son reservoir de sol

    REAL, save:: trs(klon, nbtr)    ! Concentration de radon dans le sol

    REAL masktr(klon, nbtr) ! Masque reservoir de sol traceur
    !                            Masque de l'echange avec la surface
    !                           (1 = reservoir) ou (possible => 1 )
    SAVE masktr
    REAL fshtr(klon, nbtr)  ! Flux surfacique dans le reservoir de sol
    SAVE fshtr
    REAL hsoltr(nbtr)      ! Epaisseur equivalente du reservoir de sol
    SAVE hsoltr
    REAL tautr(nbtr)       ! Constante de decroissance radioactive
    SAVE tautr
    REAL vdeptr(nbtr)      ! Vitesse de depot sec dans la couche Brownienne
    SAVE vdeptr
    REAL scavtr(nbtr)      ! Coefficient de lessivage
    SAVE scavtr

    CHARACTER itn
    INTEGER, save:: nid_tra

    ! nature du traceur

    logical aerosol(nbtr)  ! Nature du traceur
    !                            ! aerosol(it) = true  => aerosol 
    !                            ! aerosol(it) = false => gaz 
    logical clsol(nbtr)    ! couche limite sol calculée
    logical radio(nbtr)    ! décroisssance radioactive
    save aerosol, clsol, radio

    ! convection tiedtke
    INTEGER i, k, it
    REAL delp(klon, llm)

    ! Variables liees a l'ecriture de la bande histoire physique

    ! Variables locales pour effectuer les appels en serie

    REAL d_tr(klon, llm), d_trs(klon) ! tendances de traceurs 
    REAL d_tr_cl(klon, llm, nbtr) ! tendance de traceurs  couche limite
    REAL d_tr_cv(klon, llm, nbtr) ! tendance de traceurs  conv pour chq traceur
    REAL d_tr_th(klon, llm, nbtr) ! la tendance des thermiques
    REAL d_tr_dec(klon, llm, 2) ! la tendance de la decroissance 
    !                                   ! radioactive du rn - > pb 
    REAL d_tr_lessi_impa(klon, llm, nbtr) ! la tendance du lessivage 
    !                                          ! par impaction
    REAL d_tr_lessi_nucl(klon, llm, nbtr) ! la tendance du lessivage 
    !                                          ! par nucleation 
    REAL flestottr(klon, llm, nbtr) ! flux de lessivage 
    !                                    ! dans chaque couche 

    real zmasse(klon, llm) 
    ! (column-density of mass of air in a layer, in kg m-2)

    real ztra_th(klon, llm)

    character(len=20) modname
    character(len=80) abort_message
    integer isplit

    ! Controls:
    logical:: couchelimite = .true.
    logical:: convection = .true.
    logical:: lessivage = .true.
    logical, save:: inirnpb

    !--------------------------------------

    modname='phytrac'

    if (debutphy) then
       print *, 'phytrac: pdtphys = ', pdtphys
       PRINT *, 'Fréquence de sortie des traceurs : ecrit_tra = ', ecrit_tra
       if (nbtr < nqmax) then
          abort_message='See above'
          call abort_gcm(modname, abort_message, 1)
       endif
       inirnpb=rnpb

       ! Initialisation des sorties :
       call ini_histrac(nid_tra, pdtphys, presnivs, nqmax, lessivage)

       ! Initialisation de certaines variables pour le radon et le plomb 
       ! Initialisation du traceur dans le sol (couche limite radonique)
       trs(:, :) = 0.

       open (unit=99, file='starttrac', status='old', err=999, &
            form='formatted')
       read(unit=99, fmt=*) (trs(i, 1), i=1, klon)
999    continue
       close(unit=99)

       ! Initialisation de la fraction d'aerosols lessivee

       d_tr_lessi_impa(:, :, :) = 0.
       d_tr_lessi_nucl(:, :, :) = 0. 

       ! Initialisation de la nature des traceurs

       DO it = 1, nqmax
          aerosol(it) = .FALSE.  ! Tous les traceurs sont des gaz par defaut
          radio(it) = .FALSE. ! par défaut pas de passage par "radiornpb"
          clsol(it) = .FALSE.  ! Par defaut couche limite avec flux prescrit
       ENDDO

       ! Get the parameters for ozone chemistry:
       call read_coefoz
    ENDIF

    ! Initialisation du traceur dans le sol (couche limite radonique)
    if (inirnpb) THEN

       radio(1)= .true.
       radio(2)= .true.
       clsol(1)= .true.
       clsol(2)= .true.
       aerosol(2) = .TRUE. ! le Pb est un aerosol 

       call initrrnpb(ftsol, pctsrf, masktr, fshtr, hsoltr, tautr, vdeptr, &
            scavtr)
       inirnpb=.false.
    endif

    do i=1, klon
       pftsol1(i) = ftsol(i, 1)
       pftsol2(i) = ftsol(i, 2)
       pftsol3(i) = ftsol(i, 3)
       pftsol4(i) = ftsol(i, 4)

       ppsrf1(i) = pctsrf(i, 1)
       ppsrf2(i) = pctsrf(i, 2)
       ppsrf3(i) = pctsrf(i, 3)
       ppsrf4(i) = pctsrf(i, 4)

    enddo

    ! Calcul de l'effet de la convection

    if (convection) then
       DO it=1, nqmax
          if (iflag_con.eq.2) then
             ! tiedke
             CALL nflxtr(pdtphys, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
                  pplay, paprs, tr_seri(1, 1, it), d_tr_cv(1, 1, it))
          else if (iflag_con.eq.3) then
             ! KE
             call cvltr(pdtphys, da, phi, mp, paprs, pplay, &
                  tr_seri(1, 1, it), upwd, dnwd, d_tr_cv(1, 1, it))
          endif

          DO k = 1, llm
             DO i = 1, klon
                tr_seri(i, k, it) = tr_seri(i, k, it) + d_tr_cv(i, k, it)
             ENDDO
          ENDDO
          WRITE(unit=itn, fmt='(i1)') it
          CALL minmaxqfi(tr_seri(:, :, it), 0., 1.e33, &
               'convection, tracer index = ' // itn)
       ENDDO
    endif

    forall (k=1: llm) zmasse(:, k) = (paprs(:, k)-paprs(:, k+1)) / rg

    ! Calcul de l'effet des thermiques

    do it=1, nqmax
       do k=1, llm
          do i=1, klon
             d_tr_th(i, k, it)=0.
             tr_seri(i, k, it)=max(tr_seri(i, k, it), 0.)
             tr_seri(i, k, it)=min(tr_seri(i, k, it), 1.e10)
          enddo
       enddo
    enddo

    if (iflag_thermals > 0) then
       nsplit=10
       DO it=1, nqmax
          do isplit=1, nsplit
             ! Thermiques 
             call dqthermcell(klon, llm, pdtphys/nsplit &
                  , fm_therm, entr_therm, zmasse &
                  , tr_seri(1:klon, 1:llm, it), d_tr, ztra_th)

             do k=1, llm
                do i=1, klon
                   d_tr(i, k)=pdtphys*d_tr(i, k)/nsplit
                   d_tr_th(i, k, it)=d_tr_th(i, k, it)+d_tr(i, k)
                   tr_seri(i, k, it)=max(tr_seri(i, k, it)+d_tr(i, k), 0.)
                enddo
             enddo
          enddo
       ENDDO
    endif

    !   Calcul de l'effet de la couche limite

    if (couchelimite) then

       DO k = 1, llm
          DO i = 1, klon
             delp(i, k) = paprs(i, k)-paprs(i, k+1)
          ENDDO
       ENDDO

       ! MAF modif pour tenir compte du cas rnpb + traceur
       DO it=1, nqmax
          if (clsol(it)) then 
             ! couche limite avec quantite dans le sol calculee
             CALL cltracrn(it, pdtphys, yu1, yv1, &
                  coefh, t_seri, ftsol, pctsrf, &
                  tr_seri(1, 1, it), trs(1, it), &
                  paprs, pplay, delp, &
                  masktr(1, it), fshtr(1, it), hsoltr(it), &
                  tautr(it), vdeptr(it), &
                  rlat, &
                  d_tr_cl(1, 1, it), d_trs)
             DO k = 1, llm
                DO i = 1, klon
                   tr_seri(i, k, it) = tr_seri(i, k, it) + d_tr_cl(i, k, it)
                ENDDO
             ENDDO

             ! Traceur ds sol

             DO i = 1, klon
                trs(i, it) = trs(i, it) + d_trs(i)
             END DO
          else ! couche limite avec flux prescrit
             !MAF provisoire source / traceur a creer
             DO i=1, klon
                source(i) = 0.0 ! pas de source, pour l'instant
             ENDDO

             CALL cltrac(pdtphys, coefh, t_seri, &
                  tr_seri(1, 1, it), source, &
                  paprs, pplay, delp, &
                  d_tr_cl(1, 1, it))
             DO k = 1, llm
                DO i = 1, klon
                   tr_seri(i, k, it) = tr_seri(i, k, it) + d_tr_cl(i, k, it)
                ENDDO
             ENDDO
          endif
       ENDDO

    endif ! couche limite

    !   Calcul de l'effet du puits radioactif

    ! MAF il faudrait faire une modification pour passer dans radiornpb 
    ! si radio=true mais pour l'instant radiornpb propre au cas rnpb
    if (rnpb) then
       d_tr_dec(:, :, :) = radiornpb(tr_seri, pdtphys, tautr)
       DO it=1, nqmax
          if (radio(it)) then
             tr_seri(:, :, it) = tr_seri(:, :, it) + d_tr_dec(:, :, it)
             WRITE(unit=itn, fmt='(i1)') it
             CALL minmaxqfi(tr_seri(:, :, it), 0., 1.e33, 'puits rn it='//itn)
          endif
       ENDDO
    endif ! rnpb decroissance  radioactive

    ! Ozone as a tracer:
    call o3_chem(julien, gmtime, t_seri, zmasse, pdtphys, tr_seri(:, :, 3))

    ! Calcul de l'effet de la precipitation

    IF (lessivage) THEN
       d_tr_lessi_nucl(:, :, :) = 0. 
       d_tr_lessi_impa(:, :, :) = 0. 
       flestottr(:, :, :) = 0. 

       ! tendance des aerosols nuclees et impactes 

       DO it = 1, nqmax
          IF (aerosol(it)) THEN
             DO k = 1, llm
                DO i = 1, klon
                   d_tr_lessi_nucl(i, k, it) = d_tr_lessi_nucl(i, k, it) + &
                        ( 1 - frac_nucl(i, k) )*tr_seri(i, k, it)
                   d_tr_lessi_impa(i, k, it) = d_tr_lessi_impa(i, k, it) + &
                        ( 1 - frac_impa(i, k) )*tr_seri(i, k, it)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       ! Mises a jour des traceurs + calcul des flux de lessivage 
       ! Mise a jour due a l'impaction et a la nucleation

       DO it = 1, nqmax
          IF (aerosol(it)) THEN
             DO k = 1, llm
                DO i = 1, klon
                   tr_seri(i, k, it)=tr_seri(i, k, it) &
                        *frac_impa(i, k)*frac_nucl(i, k)
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       ! Flux lessivage total 

       DO it = 1, nqmax
          DO k = 1, llm
             DO i = 1, klon
                flestottr(i, k, it) = flestottr(i, k, it) - &
                     ( d_tr_lessi_nucl(i, k, it)   + &
                     d_tr_lessi_impa(i, k, it) ) * &
                     ( paprs(i, k)-paprs(i, k+1) ) /  &
                     (RG * pdtphys)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    !   Ecriture des sorties
    call write_histrac(lessivage, nqmax, nstep, nid_tra)

    if (lafin) then
       print *, "C'est la fin de la physique."
       open (unit=99, file='restarttrac',  form='formatted')
       do i=1, klon
          write(unit=99, fmt=*) trs(i, 1)
       enddo
       PRINT *, 'Ecriture du fichier restarttrac'
       close(99)
    endif

  contains

    subroutine write_histrac(lessivage, nqmax, nstep, nid_tra)

      ! From phylmd/write_histrac.h, version 1.9 2006/02/21 08:08:30

      use dimens_m, only: iim, jjm, llm
      use ioipsl, only: histwrite, histsync
      use temps, only: itau_phy
      use advtrac_m, only: tnom
      use comgeomphy, only: airephy
      use dimphy, only: klon

      logical, intent(in):: lessivage

      integer, intent(in):: nqmax
      ! (nombre de traceurs auxquels on applique la physique)

      integer, intent(in):: nstep  ! appel physique
      integer, intent(in):: nid_tra

      ! Variables local to the procedure:
      INTEGER ndex2d(iim*(jjm+1)), ndex3d(iim*(jjm+1)*llm)
      integer it
      integer itau_w   ! pas de temps ecriture
      REAL zx_tmp_2d(iim, jjm+1), zx_tmp_3d(iim, jjm+1, llm)
      logical, parameter:: ok_sync = .true.

      !-----------------------------------------------------

      ndex2d = 0
      ndex3d = 0
      itau_w = itau_phy + nstep

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, pphis, zx_tmp_2d)
      CALL histwrite(nid_tra, "phis", itau_w, zx_tmp_2d, iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, airephy, zx_tmp_2d)      
      CALL histwrite(nid_tra, "aire", itau_w, zx_tmp_2d, iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, zmasse, zx_tmp_3d)      
      CALL histwrite(nid_tra, "zmasse", itau_w, zx_tmp_3d, iim*(jjm+1)*llm, &
           ndex3d)

      DO it=1, nqmax
         CALL gr_fi_ecrit(llm, klon, iim, jjm+1, tr_seri(1, 1, it), zx_tmp_3d)
         CALL histwrite(nid_tra, tnom(it+2), itau_w, zx_tmp_3d, &
              iim*(jjm+1)*llm, ndex3d)
         if (lessivage) THEN
            CALL gr_fi_ecrit(llm, klon, iim, jjm+1, flestottr(1, 1, it), &
                 zx_tmp_3d)
            CALL histwrite(nid_tra, "fl"//tnom(it+2), itau_w, zx_tmp_3d, &
                 iim*(jjm+1)*llm, ndex3d)
         endif

         CALL gr_fi_ecrit(llm, klon, iim, jjm+1, d_tr_th(1, 1, it), zx_tmp_3d)
         CALL histwrite(nid_tra, "d_tr_th_"//tnom(it+2), itau_w, zx_tmp_3d, &
              iim*(jjm+1)*llm, ndex3d)
         CALL gr_fi_ecrit(llm, klon, iim, jjm+1, d_tr_cv(1, 1, it), zx_tmp_3d)
         CALL histwrite(nid_tra, "d_tr_cv_"//tnom(it+2), itau_w, zx_tmp_3d, &
              iim*(jjm+1)*llm, ndex3d)
         CALL gr_fi_ecrit(llm, klon, iim, jjm+1, d_tr_cl(1, 1, it), zx_tmp_3d)
         CALL histwrite(nid_tra, "d_tr_cl_"//tnom(it+2), itau_w, zx_tmp_3d, &
              iim*(jjm+1)*llm, ndex3d)
      ENDDO

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, yu1, zx_tmp_2d)
      CALL histwrite(nid_tra, "pyu1", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, yv1, zx_tmp_2d)
      CALL histwrite(nid_tra, "pyv1", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol1, zx_tmp_2d)
      CALL histwrite(nid_tra, "ftsol1", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol2, zx_tmp_2d)
      CALL histwrite(nid_tra, "ftsol2", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol3, zx_tmp_2d)
      CALL histwrite(nid_tra, "ftsol3", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, pftsol4, zx_tmp_2d)
      CALL histwrite(nid_tra, "ftsol4", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf1, zx_tmp_2d)
      CALL histwrite(nid_tra, "psrf1", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf2, zx_tmp_2d)
      CALL histwrite(nid_tra, "psrf2", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf3, zx_tmp_2d)
      CALL histwrite(nid_tra, "psrf3", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)

      CALL gr_fi_ecrit(1, klon, iim, jjm+1, ppsrf4, zx_tmp_2d)
      CALL histwrite(nid_tra, "psrf4", itau_w, zx_tmp_2d, &
           iim*(jjm+1), ndex2d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pplay, zx_tmp_3d)
      CALL histwrite(nid_tra, "pplay", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)

      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, t_seri, zx_tmp_3d)
      CALL histwrite(nid_tra, "t", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pmfu, zx_tmp_3d)
      CALL histwrite(nid_tra, "mfu", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pmfd, zx_tmp_3d)
      CALL histwrite(nid_tra, "mfd", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pen_u, zx_tmp_3d)
      CALL histwrite(nid_tra, "en_u", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pen_d, zx_tmp_3d)
      CALL histwrite(nid_tra, "en_d", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pde_d, zx_tmp_3d)
      CALL histwrite(nid_tra, "de_d", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, pde_u, zx_tmp_3d)
      CALL histwrite(nid_tra, "de_u", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)
      CALL gr_fi_ecrit(llm, klon, iim, jjm+1, coefh, zx_tmp_3d)
      CALL histwrite(nid_tra, "coefh", itau_w, zx_tmp_3d, &
           iim*(jjm+1)*llm, ndex3d)

      ! abder

      if (ok_sync) then
         call histsync(nid_tra)
      endif

    end subroutine write_histrac

  END SUBROUTINE phytrac

  !*************************************************

  subroutine ini_histrac(nid_tra, pdtphys, presnivs, nqmax, lessivage)

    ! From phylmd/ini_histrac.h, version 1.10 2006/02/21 08:08:30

    use dimens_m, only: iim, jjm, llm
    use ioipsl, only: ymds2ju, histbeg_totreg, histvert, histdef, histend
    use temps, only: annee_ref, day_ref, itau_phy
    use advtrac_m, only: niadv, tnom, ttext
    use dimphy, only: klon
    use clesphys, only: ecrit_tra
    use grid_change, only: gr_phy_write
    use phyetat0_m, only: rlon, rlat

    INTEGER, intent(out):: nid_tra
    real, intent(in):: pdtphys  ! pas d'integration pour la physique (s)
    REAL, intent(in):: presnivs(:)

    integer, intent(in):: nqmax
    ! (nombre de traceurs auxquels on applique la physique)

    logical, intent(in):: lessivage

    ! Variables local to the procedure:

    REAL zjulian
    REAL zx_lat(iim, jjm+1)
    INTEGER nhori, nvert
    REAL zsto, zout
    integer it, iq, iiq

    !---------------------------------------------------------

    CALL ymds2ju(annee_ref, month=1, day=day_ref, sec=0.0, julian=zjulian)
    zx_lat(:, :) = gr_phy_write(rlat)
    CALL histbeg_totreg("histrac", iim, rlon(2:iim+1), jjm+1, zx_lat(1, :), &
         1, iim, 1, jjm+1, itau_phy, zjulian, pdtphys, nhori, nid_tra)
    CALL histvert(nid_tra, "presnivs", "Vertical levels", "mb", llm, &
         presnivs, nvert)

    zsto = pdtphys
    zout = pdtphys * REAL(ecrit_tra)

    CALL histdef(nid_tra, "phis", "Surface geop. height", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "once",  zsto, zout)
    CALL histdef(nid_tra, "aire", "Grid area", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "once",  zsto, zout)
    CALL histdef(nid_tra, "zmasse", "column density of air in cell", &
         "kg m-2", iim, jjm + 1, nhori, llm, 1, llm, nvert, 32, "ave(X)", &
         zsto, zout)

    DO it=1, nqmax
       ! champ 2D
       iq=it+2
       iiq=niadv(iq)
       CALL histdef(nid_tra, tnom(iq), ttext(iiq), "U/kga", iim, jjm+1, &
            nhori, llm, 1, llm, nvert, 32, "ave(X)", zsto, zout)
       if (lessivage) THEN
          CALL histdef(nid_tra, "fl"//tnom(iq), "Flux "//ttext(iiq), &
               "U/m2/s", iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
               "ave(X)", zsto, zout)
       endif

       !---Ajout Olivia
       CALL histdef(nid_tra, "d_tr_th_"//tnom(iq), &
            "tendance thermique"// ttext(iiq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
            "ave(X)", zsto, zout)
       CALL histdef(nid_tra, "d_tr_cv_"//tnom(iq), &
            "tendance convection"// ttext(iiq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
            "ave(X)", zsto, zout)
       CALL histdef(nid_tra, "d_tr_cl_"//tnom(iq), &
            "tendance couche limite"// ttext(iiq), "?", &
            iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
            "ave(X)", zsto, zout)
       !---fin Olivia	 

    ENDDO

    CALL histdef(nid_tra, "pyu1", "Vent niv 1", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)", zout, zout)

    CALL histdef(nid_tra, "pyv1", "Vent niv 1", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "psrf1", "nature sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "psrf2", "nature sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "psrf3", "nature sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "psrf4", "nature sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "ftsol1", "temper sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "ftsol2", "temper sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "ftsol3", "temper sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "ftsol4", "temper sol", "-", &
         iim, jjm+1, nhori, 1, 1, 1, -99, 32, &
         "inst(X)",  zout, zout)
    CALL histdef(nid_tra, "pplay", "flux u mont", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "inst(X)", zout, zout)
    CALL histdef(nid_tra, "t", "flux u mont", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "inst(X)", zout, zout)
    CALL histdef(nid_tra, "mfu", "flux u mont", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)
    CALL histdef(nid_tra, "mfd", "flux u decen", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)
    CALL histdef(nid_tra, "en_u", "flux u mont", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)
    CALL histdef(nid_tra, "en_d", "flux u mont", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)
    CALL histdef(nid_tra, "de_d", "flux u mont", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)
    CALL histdef(nid_tra, "de_u", "flux u decen", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)
    CALL histdef(nid_tra, "coefh", "turbulent coef", "-", &
         iim, jjm+1, nhori, llm, 1, llm, nvert, 32, &
         "ave(X)", zsto, zout)

    CALL histend(nid_tra)

  end subroutine ini_histrac

  !*************************************************

  function radiornpb(tr_seri, pdtphys, tautr)

    ! From phylmd/radiornpb.F, v 1.2 2005/05/25 13:10:09

    ! Auteurs: AA + CG (LGGE/CNRS) Date 24-06-94
    ! Objet: Decroissance radioactive d'un traceur dans l'atmosphere
    !G 24 06 94 : Pour un traceur, le radon
    !G 16 12 94 : Plus un 2eme traceur, le 210Pb. Le radon decroit en plomb.

    ! Le pas de temps "pdtphys" est supposé beaucoup plus petit que la
    ! constante de temps de décroissance.

    use dimens_m, only: llm
    use dimphy, only: klon, nbtr
    use nrutil, only: assert

    IMPLICIT none

    REAL, intent(in):: tr_seri(:, :, :), pdtphys, tautr(:)
    real radiornpb(klon, llm, 2)

    ! Variable local to the procedure:
    INTEGER it

    !-----------------------------------------------

    call assert(shape(tr_seri) == (/klon, llm, nbtr/), "radiornpb tr_seri")
    call assert(size(tautr) == nbtr, "radiornpb tautr")

    DO it = 1, 2
       IF (tautr(it) > 0.) THEN
          radiornpb(:, :, it) = - tr_seri(:, :, it) * pdtphys / tautr(it)
       ELSE
          radiornpb(:, :, it) = 0.
       END IF
    END DO

    !G161294 : Cas particulier radon 1 => plomb 2
    radiornpb(:, :, 2) = radiornpb(:, :, 2) - radiornpb(:, :, 1)

  END function radiornpb

  !*************************************************

  SUBROUTINE minmaxqfi(zq, qmin, qmax, comment)

    ! From phylmd/minmaxqfi.F, version 1.1.1.1 2004/05/19 12:53:09

    use dimens_m, only: llm
    use dimphy, only: klon
    use nrutil, only: assert

    IMPLICIT none

    real, intent(in):: zq(:, :), qmin, qmax
    CHARACTER(len=*), intent(in):: comment

    ! Variables local to the procedure:

    INTEGER jadrs(klon), jbad, k, i

    !---------------------------------

    call assert(shape(zq) == (/klon, llm/), "minmaxqfi")

    DO k = 1, llm
       jbad = 0
       DO i = 1, klon
          IF (zq(i, k) > qmax .OR. zq(i, k) < qmin) THEN
             jbad = jbad + 1
             jadrs(jbad) = i
          ENDIF
       ENDDO
       IF (jbad > 0) THEN
          PRINT *, comment
          DO i = 1, jbad
             PRINT *, "zq(", jadrs(i), ", ", k, ") = ", zq(jadrs(i), k)
          ENDDO
       ENDIF
    ENDDO

  end SUBROUTINE minmaxqfi

end module phytrac_m
