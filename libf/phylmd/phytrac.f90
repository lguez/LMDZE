module phytrac_m

  ! This module is clean: no C preprocessor directive, no include line.

  IMPLICIT none

  private
  public phytrac

contains

  SUBROUTINE phytrac(rnpb, itap, lmt_pas, julien, gmtime, firstcal, lafin, &
       nq_phys, pdtphys, u, v, t_seri, paprs, pplay, pmfu, pmfd, pen_u, &
       pde_u, pen_d, pde_d, coefh, fm_therm, entr_therm, yu1, yv1, ftsol, &
       pctsrf, frac_impa, frac_nucl, pphis, pphi, albsol, rh, cldfra, rneb, &
       diafra, cldliq, itop_con, ibas_con, pmflxr, pmflxs, prfl, psfl, da, &
       phi, mp, upwd, dnwd, tr_seri, zmasse)

    ! From phylmd/phytrac.F, version 1.15 2006/02/21 08:08:30

    ! Authors: Frédéric Hourdin, Abderrahmane Idelkadi, Marie-Alice
    ! Foujols, Olivia
    ! Objet : moniteur général des tendances des traceurs

    ! Remarques :
    ! 1/ L'appel de "phytrac" se fait avec "nq-2" donc nous avons bien 
    ! les vrais traceurs (en nombre "nbtr", sans la vapeur d'eau ni l'eau
    ! liquide) dans "phytrac".
    ! 2/ Le choix du radon et du plomb se fait juste avec un "data" 
    ! (peu propre).
    ! Pourrait-on avoir une variable qui indiquerait le type de traceur ?

    use dimens_m, only: llm
    use indicesol, only: nbsrf
    use dimphy, only: klon, nbtr
    use clesphys, only: ecrit_tra
    use clesphys2, only: iflag_con
    use abort_gcm_m, only: abort_gcm
    use YOMCST, only: rg
    use ctherm, only: iflag_thermals
    use regr_pr_comb_coefoz_m, only: regr_pr_comb_coefoz
    use phyetat0_m, only: rlat
    use o3_chem_m, only: o3_chem
    use ini_hist, only: ini_histrac
    use radiornpb_m, only: radiornpb
    use minmaxqfi_m, only: minmaxqfi
    use numer_rec, only: assert
    use press_coefoz_m, only: press_coefoz

    ! Arguments:

    !   EN ENTREE:

    !   divers:

    logical, intent(in):: rnpb

    integer, intent(in):: nq_phys
    ! (nombre de traceurs auxquels on applique la physique)

    integer, intent(in):: itap  ! number of calls to "physiq"
    integer, intent(in):: lmt_pas ! number of time steps of "physics" per day
    integer, intent(in):: julien !jour julien, 1 <= julien <= 360
    integer itop_con(klon)
    integer ibas_con(klon)
    real, intent(in):: gmtime ! heure de la journée en fraction de jour
    real, intent(in):: pdtphys  ! pas d'integration pour la physique (s)
    real, intent(in):: t_seri(klon, llm) ! temperature, in K

    real, intent(inout):: tr_seri(:, :, :) ! (klon, llm, nbtr)
    ! (mass fractions of tracers, excluding water, at mid-layers)

    real u(klon, llm)
    real v(klon, llm)
    real rh(klon, llm)     ! humidite relative
    real cldliq(klon, llm) ! eau liquide nuageuse
    real cldfra(klon, llm) ! fraction nuageuse (tous les nuages)

    real diafra(klon, llm)
    ! (fraction nuageuse (convection ou stratus artificiels))

    real rneb(klon, llm)   ! fraction nuageuse (grande echelle)
    real albsol(klon)  ! albedo surface

    real, intent(in):: paprs(klon, llm+1)
    ! (pression pour chaque inter-couche, en Pa)

    real, intent(in):: pplay(klon, llm)
    ! (pression pour le mileu de chaque couche, en Pa)

    real pphi(klon, llm) ! geopotentiel
    real pphis(klon)
    logical, intent(in):: firstcal ! first call to "calfis"
    logical, intent(in):: lafin ! fin de la physique

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

    real, intent(in):: zmasse(:, :)  ! (klon, llm)
    ! (column-density of mass of air in a cell, in kg m-2)

    ! Variables local to the procedure:

    integer nsplit

    !  TRACEURS

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

    real ztra_th(klon, llm)
    integer isplit

    ! Controls:
    logical:: couchelimite = .true.
    logical:: convection = .true.
    logical:: lessivage = .true.
    logical, save:: inirnpb

    !--------------------------------------

    call assert(shape(zmasse) == (/klon, llm/), "phytrac zmasse")
    call assert(shape(tr_seri) == (/klon, llm, nbtr/), "phytrac tr_seri")

    if (firstcal) then
       print *, 'phytrac: pdtphys = ', pdtphys
       PRINT *, 'Fréquence de sortie des traceurs : ecrit_tra = ', ecrit_tra
       if (nbtr < nq_phys) call abort_gcm('phytrac', 'nbtr < nq_phys', 1)
       inirnpb=rnpb

       ! Initialisation des sorties :
       call ini_histrac(nid_tra, pdtphys, nq_phys, lessivage)

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

       DO it = 1, nq_phys
          aerosol(it) = .FALSE.  ! Tous les traceurs sont des gaz par defaut
          radio(it) = .FALSE. ! par défaut pas de passage par "radiornpb"
          clsol(it) = .FALSE.  ! Par defaut couche limite avec flux prescrit
       ENDDO

       if (nq_phys >= 3) then
          call press_coefoz ! read input pressure levels for ozone coefficients
       end if
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

    ! Calcul de l'effet de la convection

    if (convection) then
       DO it=1, nq_phys
          if (iflag_con.eq.2) then
             ! tiedke
             CALL nflxtr(pdtphys, pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
                  paprs, tr_seri(1, 1, it), d_tr_cv(1, 1, it))
          else if (iflag_con.eq.3) then
             ! KE
             call cvltr(pdtphys, da, phi, mp, paprs, &
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

    ! Calcul de l'effet des thermiques

    do it=1, nq_phys
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
       DO it=1, nq_phys
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
       DO it=1, nq_phys
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
       DO it=1, nq_phys
          if (radio(it)) then
             tr_seri(:, :, it) = tr_seri(:, :, it) + d_tr_dec(:, :, it)
             WRITE(unit=itn, fmt='(i1)') it
             CALL minmaxqfi(tr_seri(:, :, it), 0., 1.e33, 'puits rn it='//itn)
          endif
       ENDDO
    endif ! rnpb decroissance  radioactive

    if (nq_phys >= 3) then
       ! Ozone as a tracer:
       if (mod(itap - 1, lmt_pas) == 0) then
          ! Once per day, update the coefficients for ozone chemistry:
          call regr_pr_comb_coefoz(julien)
       end if
       call o3_chem(julien, gmtime, t_seri, zmasse, pdtphys, tr_seri(:, :, 3))
    end if

    ! Calcul de l'effet de la precipitation

    IF (lessivage) THEN
       d_tr_lessi_nucl(:, :, :) = 0. 
       d_tr_lessi_impa(:, :, :) = 0. 
       flestottr(:, :, :) = 0. 

       ! tendance des aerosols nuclees et impactes 

       DO it = 1, nq_phys
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

       DO it = 1, nq_phys
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

       DO it = 1, nq_phys
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
    call write_histrac(lessivage, nq_phys, itap, nid_tra)

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

    subroutine write_histrac(lessivage, nq_phys, itap, nid_tra)

      ! From phylmd/write_histrac.h, version 1.9 2006/02/21 08:08:30

      use dimens_m, only: iim, jjm, llm
      use ioipsl, only: histwrite, histsync
      use temps, only: itau_phy
      use iniadvtrac_m, only: tnom
      use comgeomphy, only: airephy
      use dimphy, only: klon
      use grid_change, only: gr_phy_write_2d, gr_phy_write_3d

      logical, intent(in):: lessivage

      integer, intent(in):: nq_phys
      ! (nombre de traceurs auxquels on applique la physique)

      integer, intent(in):: itap  ! number of calls to "physiq"
      integer, intent(in):: nid_tra

      ! Variables local to the procedure:
      integer it
      integer itau_w   ! pas de temps ecriture
      REAL zx_tmp_2d(iim, jjm+1), zx_tmp_3d(iim, jjm+1, llm)
      logical, parameter:: ok_sync = .true.

      !-----------------------------------------------------

      itau_w = itau_phy + itap

      CALL histwrite(nid_tra, "phis", itau_w, gr_phy_write_2d(pphis))
      CALL histwrite(nid_tra, "aire", itau_w, gr_phy_write_2d(airephy))
      CALL histwrite(nid_tra, "zmasse", itau_w, gr_phy_write_3d(zmasse))

      DO it=1, nq_phys
         CALL histwrite(nid_tra, tnom(it+2), itau_w, &
              gr_phy_write_3d(tr_seri(:, :, it)))
         if (lessivage) THEN
            CALL histwrite(nid_tra, "fl"//tnom(it+2), itau_w, &
                 gr_phy_write_3d(flestottr(:, :, it)))
         endif
         CALL histwrite(nid_tra, "d_tr_th_"//tnom(it+2), itau_w, &
              gr_phy_write_3d(d_tr_th(:, :, it)))
         CALL histwrite(nid_tra, "d_tr_cv_"//tnom(it+2), itau_w, &
              gr_phy_write_3d(d_tr_cv(:, :, it)))
         CALL histwrite(nid_tra, "d_tr_cl_"//tnom(it+2), itau_w, &
              gr_phy_write_3d(d_tr_cl(:, :, it)))
      ENDDO

      CALL histwrite(nid_tra, "pplay", itau_w, gr_phy_write_3d(pplay))
      CALL histwrite(nid_tra, "t", itau_w, gr_phy_write_3d(t_seri))

      if (ok_sync) then
         call histsync(nid_tra)
      endif

    end subroutine write_histrac

  END SUBROUTINE phytrac

end module phytrac_m
