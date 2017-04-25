module phytrac_m

  IMPLICIT none

  private
  public phytrac

contains

  SUBROUTINE phytrac(julien, gmtime, firstcal, lafin, pdtphys, t_seri, paprs, &
       pplay, pmfu, pmfd, pde_u, pen_d, coefh, fm_therm, entr_therm, yu1, &
       yv1, ftsol, pctsrf, frac_impa, frac_nucl, da, phi, mp, upwd, dnwd, &
       tr_seri, zmasse, ncid_startphy)

    ! From phylmd/phytrac.F, version 1.15, 2006/02/21 08:08:30 (SVN
    ! revision 679) and phylmd/write_histrac.h, version 1.9,
    ! 2006/02/21 08:08:30

    ! Authors: Fr\'ed\'eric Hourdin, Abderrahmane Idelkadi, Marie-Alice
    ! Foujols

    ! Objet : moniteur g\'en\'eral des tendances des traceurs

    ! L'appel de "phytrac" se fait avec "nqmx - 2" donc nous avons
    ! bien les vrais traceurs, sans la vapeur d'eau ni l'eau liquide.

    ! Modifications pour les traceurs :
    ! - uniformisation des param\'etrisations dans phytrac
    ! - stockage des moyennes des champs n\'ecessaires en mode traceur off-line 

    use abort_gcm_m, only: abort_gcm
    use clesphys2, only: conv_emanuel
    use cltrac_m, only: cltrac
    use cltracrn_m, only: cltracrn
    USE conf_gcm_m, ONLY: lmt_pas
    use ctherm, only: iflag_thermals
    use cvltr_m, only: cvltr
    use dimens_m, only: llm, nqmx
    use dimphy, only: klon
    use histwrite_phy_m, only: histwrite_phy
    use indicesol, only: nbsrf
    use iniadvtrac_m, only: tname
    use initrrnpb_m, only: initrrnpb
    use minmaxqfi_m, only: minmaxqfi
    use netcdf, only: NF90_FILL_float
    use netcdf95, only: nf95_inq_varid, nf95_get_var, nf95_put_var
    use nflxtr_m, only: nflxtr
    use nr_util, only: assert
    use o3_chem_m, only: o3_chem
    use phyetat0_m, only: rlat
    use phyredem0_m, only: ncid_restartphy
    use press_coefoz_m, only: press_coefoz
    use radiornpb_m, only: radiornpb
    use regr_pr_comb_coefoz_m, only: regr_pr_comb_coefoz
    use SUPHEC_M, only: rg
    use time_phylmdz, only: itap

    integer, intent(in):: julien !jour julien, 1 <= julien <= 360
    real, intent(in):: gmtime ! heure de la journ\'ee en fraction de jour
    logical, intent(in):: firstcal ! first call to "calfis"
    logical, intent(in):: lafin ! fin de la physique
    real, intent(in):: pdtphys ! pas d'integration pour la physique (s)
    real, intent(in):: t_seri(klon, llm) ! temperature, in K

    real, intent(in):: paprs(klon, llm+1)
    ! (pression pour chaque inter-couche, en Pa)

    real, intent(in):: pplay(klon, llm)
    ! (pression pour le mileu de chaque couche, en Pa)

    ! convection:

    REAL, intent(in):: pmfu(klon, llm) ! flux de masse dans le panache montant

    REAL, intent(in):: pmfd(klon, llm)
    ! flux de masse dans le panache descendant

    REAL pde_u(klon, llm) ! flux detraine dans le panache montant
    REAL pen_d(klon, llm) ! flux entraine dans le panache descendant
    REAL coefh(klon, llm) ! coeff melange couche limite

    ! thermiques:
    real fm_therm(klon, llm+1), entr_therm(klon, llm)

    ! Couche limite:
    REAL yu1(klon) ! vents au premier niveau
    REAL yv1(klon) ! vents au premier niveau

    ! Arguments n\'ecessaires pour les sources et puits de traceur :
    real, intent(in):: ftsol(:, :) ! (klon, nbsrf) surface temperature (K)
    real, intent(in):: pctsrf(klon, nbsrf) ! Pourcentage de sol f(nature du sol)

    ! Lessivage pour le on-line
    REAL, intent(in):: frac_impa(klon, llm) ! fraction d'aerosols impactes
    REAL, intent(in):: frac_nucl(klon, llm) ! fraction d'aerosols nuclees

    ! Kerry Emanuel
    real, intent(in):: da(klon, llm), phi(klon, llm, llm), mp(klon, llm)
    REAL, intent(in):: upwd(klon, llm) ! saturated updraft mass flux
    REAL, intent(in):: dnwd(klon, llm) ! saturated downdraft mass flux

    real, intent(inout):: tr_seri(:, :, :) ! (klon, llm, nqmx - 2)
    ! (mass fractions of tracers, excluding water, at mid-layers)

    real, intent(in):: zmasse(:, :) ! (klon, llm)
    ! (column-density of mass of air in a cell, in kg m-2)

    integer, intent(in):: ncid_startphy

    ! Local:

    integer nsplit

    ! TRACEURS

    ! Sources et puits des traceurs:

    ! Pour l'instant seuls les cas du rn et du pb ont ete envisages.

    REAL source(klon) ! a voir lorsque le flux est prescrit 
    ! 
    ! Pour la source de radon et son reservoir de sol

    REAL, save:: trs(klon, nqmx - 2) ! Concentration de traceur dans le sol

    REAL masktr(klon, nqmx - 2) ! Masque reservoir de sol traceur
    ! Masque de l'echange avec la surface
    ! (1 = reservoir) ou (possible => 1)
    SAVE masktr
    REAL fshtr(klon, nqmx - 2) ! Flux surfacique dans le reservoir de sol
    SAVE fshtr
    REAL hsoltr(nqmx - 2) ! Epaisseur equivalente du reservoir de sol
    SAVE hsoltr
    REAL tautr(nqmx - 2) ! Constante de decroissance radioactive
    SAVE tautr
    REAL vdeptr(nqmx - 2) ! Vitesse de depot sec dans la couche Brownienne
    SAVE vdeptr
    REAL scavtr(nqmx - 2) ! Coefficient de lessivage
    SAVE scavtr

    CHARACTER itn

    ! nature du traceur

    logical aerosol(nqmx - 2) ! Nature du traceur
    ! ! aerosol(it) = true => aerosol 
    ! ! aerosol(it) = false => gaz 
    logical clsol(nqmx - 2) ! couche limite sol calcul\'ee
    logical radio(nqmx - 2) ! d\'ecroisssance radioactive
    save aerosol, clsol, radio

    ! convection tiedtke
    INTEGER i, k, it
    REAL delp(klon, llm)

    ! Variables liees a l'ecriture de la bande histoire physique

    ! Variables locales pour effectuer les appels en serie

    REAL d_tr(klon, llm), d_trs(klon) ! tendances de traceurs 
    REAL d_tr_cl(klon, llm, nqmx - 2) ! tendance de traceurs couche limite

    REAL d_tr_cv(klon, llm, nqmx - 2) 
    ! tendance de traceurs conv pour chq traceur

    REAL d_tr_th(klon, llm, nqmx - 2) ! la tendance des thermiques
    REAL d_tr_dec(klon, llm, 2) ! la tendance de la decroissance 
    ! ! radioactive du rn - > pb 

    REAL d_tr_lessi_impa(klon, llm, nqmx - 2)
    ! tendance du lessivage par impaction

    REAL d_tr_lessi_nucl(klon, llm, nqmx - 2)
    ! tendance du lessivage par nucleation

    REAL flestottr(klon, llm, nqmx - 2) ! flux de lessivage dans chaque couche 

    real ztra_th(klon, llm)
    integer isplit, varid

    ! Controls:
    logical:: couchelimite = .true.
    logical:: convection = .true.
    logical, save:: inirnpb

    !--------------------------------------

    call assert(shape(zmasse) == (/klon, llm/), "phytrac zmasse")
    call assert(shape(tr_seri) == (/klon, llm, nqmx - 2/), "phytrac tr_seri")

    if (firstcal) then
       inirnpb = .true.

       ! Initialisation de certaines variables pour le radon et le plomb 
       ! Initialisation du traceur dans le sol (couche limite radonique)
       trs(:, 2:) = 0.

       call nf95_inq_varid(ncid_startphy, "trs", varid)
       call nf95_get_var(ncid_startphy, varid, trs(:, 1))
       if (any(trs(:, 1) == NF90_FILL_float)) call abort_gcm("phytrac", &
            "some missing values in trs(:, 1)")

       ! Initialisation de la fraction d'aerosols lessivee

       d_tr_lessi_impa = 0.
       d_tr_lessi_nucl = 0. 

       ! Initialisation de la nature des traceurs

       DO it = 1, nqmx - 2
          aerosol(it) = .FALSE. ! Tous les traceurs sont des gaz par defaut
          radio(it) = .FALSE. ! par d\'efaut pas de passage par "radiornpb"
          clsol(it) = .FALSE. ! Par defaut couche limite avec flux prescrit
       ENDDO

       if (nqmx >= 5) then
          call press_coefoz ! read input pressure levels for ozone coefficients
       end if
    ENDIF

    if (inirnpb) THEN
       ! Initialisation du traceur dans le sol (couche limite radonique)
       radio(1)= .true.
       radio(2)= .true.
       clsol(1)= .true.
       clsol(2)= .true.
       aerosol(2) = .TRUE. ! le Pb est un aerosol 
       call initrrnpb(pctsrf, masktr, fshtr, hsoltr, tautr, vdeptr, scavtr)
       inirnpb=.false.
    endif

    if (convection) then
       ! Calcul de l'effet de la convection
       DO it=1, nqmx - 2
          if (conv_emanuel) then
             call cvltr(pdtphys, da, phi, mp, paprs, tr_seri(:, :, it), upwd, &
                  dnwd, d_tr_cv(:, :, it))
          else
             CALL nflxtr(pdtphys, pmfu, pmfd, pde_u, pen_d, paprs, &
                  tr_seri(:, :, it), d_tr_cv(:, :, it))
          endif

          DO k = 1, llm
             DO i = 1, klon
                tr_seri(i, k, it) = tr_seri(i, k, it) + d_tr_cv(i, k, it)
             ENDDO
          ENDDO
          WRITE(unit=itn, fmt='(i1)') it
          CALL minmaxqfi(tr_seri(:, :, it), 0., 1e33, &
               'convection, tracer index = ' // itn)
       ENDDO
    endif

    ! Calcul de l'effet des thermiques

    do it=1, nqmx - 2
       do k=1, llm
          do i=1, klon
             d_tr_th(i, k, it)=0.
             tr_seri(i, k, it)=max(tr_seri(i, k, it), 0.)
             tr_seri(i, k, it)=min(tr_seri(i, k, it), 1e10)
          enddo
       enddo
    enddo

    if (iflag_thermals > 0) then
       nsplit=10
       DO it=1, nqmx - 2
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

    ! Calcul de l'effet de la couche limite

    if (couchelimite) then
       DO k = 1, llm
          DO i = 1, klon
             delp(i, k) = paprs(i, k)-paprs(i, k+1)
          ENDDO
       ENDDO

       ! MAF modif pour tenir compte du cas traceur
       DO it=1, nqmx - 2
          if (clsol(it)) then 
             ! couche limite avec quantite dans le sol calculee
             CALL cltracrn(it, pdtphys, yu1, yv1, coefh, t_seri, ftsol, &
                  pctsrf, tr_seri(:, :, it), trs(:, it), paprs, pplay, delp, &
                  masktr(1, it), fshtr(1, it), hsoltr(it), tautr(it), &
                  vdeptr(it), rlat, d_tr_cl(1, 1, it), d_trs)
             DO k = 1, llm
                DO i = 1, klon
                   tr_seri(i, k, it) = tr_seri(i, k, it) + d_tr_cl(i, k, it)
                ENDDO
             ENDDO

             trs(:, it) = trs(:, it) + d_trs
          else
             ! couche limite avec flux prescrit
             !MAF provisoire source / traceur a creer
             DO i=1, klon
                source(i) = 0. ! pas de source, pour l'instant
             ENDDO

             CALL cltrac(pdtphys, coefh, t_seri, tr_seri(:, :, it), source, &
                  paprs, pplay, delp, d_tr_cl(1, 1, it))
             DO k = 1, llm
                DO i = 1, klon
                   tr_seri(i, k, it) = tr_seri(i, k, it) + d_tr_cl(i, k, it)
                ENDDO
             ENDDO
          endif
       ENDDO
    endif

    ! Calcul de l'effet du puits radioactif

    ! MAF il faudrait faire une modification pour passer dans radiornpb 
    ! si radio=true
    d_tr_dec = radiornpb(tr_seri, pdtphys, tautr)
    DO it = 1, nqmx - 2
       if (radio(it)) then
          tr_seri(:, :, it) = tr_seri(:, :, it) + d_tr_dec(:, :, it)
          WRITE(unit=itn, fmt='(i1)') it
          CALL minmaxqfi(tr_seri(:, :, it), 0., 1e33, 'puits rn it='//itn)
       endif
    ENDDO

    if (nqmx >= 5) then
       ! Ozone as a tracer:
       if (mod(itap - 1, lmt_pas) == 0) then
          ! Once per day, update the coefficients for ozone chemistry:
          call regr_pr_comb_coefoz(julien, paprs, pplay)
       end if
       call o3_chem(julien, gmtime, t_seri, zmasse, pdtphys, tr_seri(:, :, 3))
    end if

    ! Calcul de l'effet de la precipitation

    d_tr_lessi_nucl = 0. 
    d_tr_lessi_impa = 0. 
    flestottr = 0. 

    ! tendance des aerosols nuclees et impactes 

    DO it = 1, nqmx - 2
       IF (aerosol(it)) THEN
          DO k = 1, llm
             DO i = 1, klon
                d_tr_lessi_nucl(i, k, it) = d_tr_lessi_nucl(i, k, it) + &
                     (1 - frac_nucl(i, k))*tr_seri(i, k, it)
                d_tr_lessi_impa(i, k, it) = d_tr_lessi_impa(i, k, it) + &
                     (1 - frac_impa(i, k))*tr_seri(i, k, it)
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    ! Mises a jour des traceurs + calcul des flux de lessivage 
    ! Mise a jour due a l'impaction et a la nucleation

    DO it = 1, nqmx - 2
       IF (aerosol(it)) THEN
          DO k = 1, llm
             DO i = 1, klon
                tr_seri(i, k, it) = tr_seri(i, k, it) * frac_impa(i, k) &
                     * frac_nucl(i, k)
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    ! Flux lessivage total 
    DO it = 1, nqmx - 2
       DO k = 1, llm
          DO i = 1, klon
             flestottr(i, k, it) = flestottr(i, k, it) &
                  - (d_tr_lessi_nucl(i, k, it) + d_tr_lessi_impa(i, k, it)) &
                  * (paprs(i, k)-paprs(i, k+1)) / (RG * pdtphys)
          ENDDO
       ENDDO
    ENDDO

    ! Ecriture des sorties
    CALL histwrite_phy("zmasse", zmasse)
    DO it=1, nqmx - 2
       CALL histwrite_phy(tname(it+2), tr_seri(:, :, it))
       CALL histwrite_phy("fl"//tname(it+2), flestottr(:, :, it))
       CALL histwrite_phy("d_tr_th_"//tname(it+2), d_tr_th(:, :, it))
       CALL histwrite_phy("d_tr_cv_"//tname(it+2), d_tr_cv(:, :, it))
       CALL histwrite_phy("d_tr_cl_"//tname(it+2), d_tr_cl(:, :, it))
    ENDDO

    if (lafin) then
       call nf95_inq_varid(ncid_restartphy, "trs", varid)
       call nf95_put_var(ncid_restartphy, varid, trs(:, 1))
    endif

  END SUBROUTINE phytrac

end module phytrac_m
