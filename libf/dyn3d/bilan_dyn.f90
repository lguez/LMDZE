module bilan_dyn_m

  IMPLICIT NONE

contains

  SUBROUTINE bilan_dyn(ps, masse, pk, flux_u, flux_v, teta, phi, ucov, vcov, &
       trac, dt_app)

    ! From LMDZ4/libf/dyn3d/bilan_dyn.F, version 1.5 2005/03/16 10:12:17

    ! Sous-programme consacré à des diagnostics dynamiques de base.
    ! De façon générale, les moyennes des scalaires Q sont pondérées
    ! par la masse. Les flux de masse sont, eux, simplement moyennés.

    USE calendar, ONLY: ymds2ju
    USE conf_gcm_m, ONLY: day_step, iperiod, periodav
    USE comconst, ONLY: cpp
    USE comvert, ONLY: presnivs
    USE comgeom, ONLY: constang_2d, cu_2d, cv_2d, rlatv
    USE dimens_m, ONLY: iim, jjm, llm
    USE histcom, ONLY: histbeg_totreg, histdef, histend, histvert
    USE histwrite_m, ONLY: histwrite
    USE nr_util, ONLY: pi
    USE paramet_m, ONLY: iip1, jjp1
    USE temps, ONLY: annee_ref, day_ref, itau_dyn

    ! Arguments:

    real, intent(in):: ps(iip1, jjp1)
    real, intent(in):: masse(iip1, jjp1, llm), pk(iip1, jjp1, llm)
    real, intent(in):: flux_u(iip1, jjp1, llm)
    real, intent(in):: flux_v(iip1, jjm, llm)
    real, intent(in):: teta(iip1, jjp1, llm)
    real, intent(in):: phi(iip1, jjp1, llm)
    real, intent(in):: ucov(iip1, jjp1, llm)
    real, intent(in):: vcov(iip1, jjm, llm)
    real, intent(in):: trac(:, :, :) ! (iim + 1, jjm + 1, llm)
    real, intent(in):: dt_app

    ! Local:

    real dt_cum
    integer:: icum  = 0
    integer, save:: ncum
    logical:: first = .true.
    real zqy, zfactv(jjm, llm)

    integer, parameter:: nQ = 7
    character(len=4), parameter:: nom(nQ) = (/'T   ', 'gz  ', 'K   ', 'ang ', &
         'u   ', 'ovap', 'un  '/)
    character(len=5), parameter:: unites(nQ) = (/'K    ', 'm2/s2', 'm2/s2', &
         'ang  ', 'm/s  ', 'kg/kg', 'un   '/)

    integer:: itau = 0
    real ww

    ! Variables dynamiques intermédiaires
    REAL vcont(iip1, jjm, llm), ucont(iip1, jjp1, llm)
    REAL ang(iip1, jjp1, llm), unat(iip1, jjp1, llm)
    REAL massebx(iip1, jjp1, llm), masseby(iip1, jjm, llm)
    REAL w(iip1, jjp1, llm), ecin(iip1, jjp1, llm), convm(iip1, jjp1, llm)

    ! Champ contenant les scalaires advectés
    real Q(iip1, jjp1, llm, nQ)

    ! Champs cumulés
    real, save:: ps_cum(iip1, jjp1)
    real, save:: masse_cum(iip1, jjp1, llm)
    real, save:: flux_u_cum(iip1, jjp1, llm)
    real, save:: flux_v_cum(iip1, jjm, llm)
    real, save:: Q_cum(iip1, jjp1, llm, nQ)
    real, save:: flux_uQ_cum(iip1, jjp1, llm, nQ)
    real, save:: flux_vQ_cum(iip1, jjm, llm, nQ)
    real dQ(iip1, jjp1, llm, nQ)

    ! champs de tansport en moyenne zonale
    integer itr
    integer, parameter:: ntr = 5

    character(len=10), save:: znom(ntr, nQ)
    character(len=26), save:: znoml(ntr, nQ)
    character(len=12), save:: zunites(ntr, nQ)

    integer, parameter:: iave = 1, itot = 2, immc = 3, itrs = 4, istn = 5
    character(len=3), parameter:: ctrs(ntr) = (/'   ', 'TOT', 'MMC', 'TRS', &
         'STN'/)

    real zvQ(jjm, llm, ntr, nQ), zvQtmp(jjm, llm)
    real zavQ(jjm, 2: ntr, nQ), psiQ(jjm, llm + 1, nQ)
    real zmasse(jjm, llm)
    real zv(jjm, llm), psi(jjm, llm + 1)
    integer i, j, l, iQ

    ! Initialisation du fichier contenant les moyennes zonales.

    integer, save:: fileid
    integer thoriid, zvertiid

    real zjulian
    integer zan, dayref

    real rlong(jjm), rlatg(jjm)

    !-----------------------------------------------------------------

    !!print *, "Call sequence information: bilan_dyn"

    first_call: if (first) then
       ! initialisation des fichiers
       first = .false.
       ! ncum est la frequence de stokage en pas de temps
       ncum = day_step / iperiod * periodav
       dt_cum = ncum * dt_app

       ! Initialisation du fichier contenant les moyennes zonales

       zan = annee_ref
       dayref = day_ref
       CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)

       rlong = 0.
       rlatg = rlatv*180./pi

       call histbeg_totreg('dynzon', rlong(:1), rlatg, 1, 1, 1, jjm, itau_dyn, &
            zjulian, dt_cum, thoriid, fileid)

       ! Appel à histvert pour la grille verticale

       call histvert(fileid, 'presnivs', 'Niveaux sigma', 'mb', llm, presnivs, &
            zvertiid)

       ! Appels à histdef pour la définition des variables à sauvegarder
       do iQ = 1, nQ
          do itr = 1, ntr
             if (itr == 1) then
                znom(itr, iQ) = nom(iQ)
                znoml(itr, iQ) = nom(iQ)
                zunites(itr, iQ) = unites(iQ)
             else
                znom(itr, iQ) = ctrs(itr)//'v'//nom(iQ)
                znoml(itr, iQ) = 'transport : v * '//nom(iQ)//' '//ctrs(itr)
                zunites(itr, iQ) = 'm/s * '//unites(iQ)
             endif
          enddo
       enddo

       ! Déclarations des champs avec dimension verticale
       do iQ = 1, nQ
          do itr = 1, ntr
             call histdef(fileid, znom(itr, iQ), znoml(itr, iQ), &
                  zunites(itr, iQ), 1, jjm, thoriid, llm, 1, llm, zvertiid, &
                  'ave(X)', dt_cum, dt_cum)
          enddo
          ! Déclarations pour les fonctions de courant
          call histdef(fileid, 'psi'//nom(iQ), 'stream fn. '//znoml(itot, iQ), &
               zunites(itot, iQ), 1, jjm, thoriid, llm, 1, llm, zvertiid, &
               'ave(X)', dt_cum, dt_cum)
       enddo

       ! Déclarations pour les champs de transport d'air
       call histdef(fileid, 'masse', 'masse', &
            'kg', 1, jjm, thoriid, llm, 1, llm, zvertiid, &
            'ave(X)', dt_cum, dt_cum)
       call histdef(fileid, 'v', 'v', &
            'm/s', 1, jjm, thoriid, llm, 1, llm, zvertiid, &
            'ave(X)', dt_cum, dt_cum)
       ! Déclarations pour les fonctions de courant
       call histdef(fileid, 'psi', 'stream fn. MMC ', 'mega t/s', &
            1, jjm, thoriid, llm, 1, llm, zvertiid, &
            'ave(X)', dt_cum, dt_cum)

       ! Déclaration des champs 1D de transport en latitude
       do iQ = 1, nQ
          do itr = 2, ntr
             call histdef(fileid, 'a'//znom(itr, iQ), znoml(itr, iQ), &
                  zunites(itr, iQ), 1, jjm, thoriid, 1, 1, 1, -99, &
                  'ave(X)', dt_cum, dt_cum)
          enddo
       enddo

       CALL histend(fileid)
    endif first_call

    itau = itau + 1

    ! Calcul des champs dynamiques

    ! Énergie cinétique
    ucont = 0
    CALL covcont(llm, ucov, vcov, ucont, vcont)
    CALL enercin(vcov, ucov, vcont, ucont, ecin)

    ! moment cinétique
    do l = 1, llm
       ang(:, :, l) = ucov(:, :, l) + constang_2d
       unat(:, :, l) = ucont(:, :, l)*cu_2d
    enddo

    Q(:, :, :, 1) = teta * pk / cpp
    Q(:, :, :, 2) = phi
    Q(:, :, :, 3) = ecin
    Q(:, :, :, 4) = ang
    Q(:, :, :, 5) = unat
    Q(:, :, :, 6) = trac
    Q(:, :, :, 7) = 1.

    ! Cumul

    if (icum == 0) then
       ps_cum = 0.
       masse_cum = 0.
       flux_u_cum = 0.
       flux_v_cum = 0.
       Q_cum = 0.
       flux_vQ_cum = 0.
       flux_uQ_cum = 0.
    endif

    icum = icum + 1

    ! Accumulation des flux de masse horizontaux
    ps_cum = ps_cum + ps
    masse_cum = masse_cum + masse
    flux_u_cum = flux_u_cum + flux_u
    flux_v_cum = flux_v_cum + flux_v
    do iQ = 1, nQ
       Q_cum(:, :, :, iQ) = Q_cum(:, :, :, iQ) + Q(:, :, :, iQ)*masse
    enddo

    ! FLUX ET TENDANCES

    ! Flux longitudinal
    forall (iQ = 1: nQ, i = 1: iim) flux_uQ_cum(i, :, :, iQ) &
         = flux_uQ_cum(i, :, :, iQ) &
         + flux_u(i, :, :) * 0.5 * (Q(i, :, :, iQ) + Q(i + 1, :, :, iQ))
    flux_uQ_cum(iip1, :, :, :) = flux_uQ_cum(1, :, :, :)

    ! Flux méridien
    forall (iQ = 1: nQ, j = 1: jjm) flux_vQ_cum(:, j, :, iQ) &
         = flux_vQ_cum(:, j, :, iQ) &
         + flux_v(:, j, :) * 0.5 * (Q(:, j, :, iQ) + Q(:, j + 1, :, iQ))

    ! tendances

    ! convergence horizontale
    call convflu(flux_uQ_cum, flux_vQ_cum, llm*nQ, dQ)

    ! calcul de la vitesse verticale
    call convmas(flux_u_cum, flux_v_cum, convm)
    CALL vitvert(convm, w)

    do iQ = 1, nQ
       do l = 1, llm-1
          do j = 1, jjp1
             do i = 1, iip1
                ww = -0.5*w(i, j, l + 1)*(Q(i, j, l, iQ) + Q(i, j, l + 1, iQ))
                dQ(i, j, l, iQ) = dQ(i, j, l, iQ)-ww
                dQ(i, j, l + 1, iQ) = dQ(i, j, l + 1, iQ) + ww
             enddo
          enddo
       enddo
    enddo

    ! PAS DE TEMPS D'ECRITURE

    writing_step: if (icum == ncum) then
       ! Normalisation
       do iQ = 1, nQ
          Q_cum(:, :, :, iQ) = Q_cum(:, :, :, iQ)/masse_cum
       enddo
       ps_cum = ps_cum / ncum
       masse_cum = masse_cum / ncum
       flux_u_cum = flux_u_cum / ncum
       flux_v_cum = flux_v_cum / ncum
       flux_uQ_cum = flux_uQ_cum / ncum
       flux_vQ_cum = flux_vQ_cum / ncum
       dQ = dQ / ncum

       ! A retravailler eventuellement
       ! division de dQ par la masse pour revenir aux bonnes grandeurs
       do iQ = 1, nQ
          dQ(:, :, :, iQ) = dQ(:, :, :, iQ)/masse_cum
       enddo

       ! Transport méridien

       ! cumul zonal des masses des mailles

       zv = 0.
       zmasse = 0.
       call massbar(masse_cum, massebx, masseby)
       do l = 1, llm
          do j = 1, jjm
             do i = 1, iim
                zmasse(j, l) = zmasse(j, l) + masseby(i, j, l)
                zv(j, l) = zv(j, l) + flux_v_cum(i, j, l)
             enddo
             zfactv(j, l) = cv_2d(1, j)/zmasse(j, l)
          enddo
       enddo

       ! Transport dans le plan latitude-altitude

       zvQ = 0.
       psiQ = 0.
       do iQ = 1, nQ
          zvQtmp = 0.
          do l = 1, llm
             do j = 1, jjm
                ! Calcul des moyennes zonales du transort total et de zvQtmp
                do i = 1, iim
                   zvQ(j, l, itot, iQ) = zvQ(j, l, itot, iQ) &
                        + flux_vQ_cum(i, j, l, iQ)
                   zqy =  0.5 * (Q_cum(i, j, l, iQ) * masse_cum(i, j, l) &
                        + Q_cum(i, j + 1, l, iQ) * masse_cum(i, j + 1, l))
                   zvQtmp(j, l) = zvQtmp(j, l) + flux_v_cum(i, j, l) * zqy &
                        / (0.5 * (masse_cum(i, j, l) + masse_cum(i, j + 1, l)))
                   zvQ(j, l, iave, iQ) = zvQ(j, l, iave, iQ) + zqy
                enddo
                ! Decomposition
                zvQ(j, l, iave, iQ) = zvQ(j, l, iave, iQ)/zmasse(j, l)
                zvQ(j, l, itot, iQ) = zvQ(j, l, itot, iQ)*zfactv(j, l)
                zvQtmp(j, l) = zvQtmp(j, l)*zfactv(j, l)
                zvQ(j, l, immc, iQ) = zv(j, l)*zvQ(j, l, iave, iQ)*zfactv(j, l)
                zvQ(j, l, itrs, iQ) = zvQ(j, l, itot, iQ)-zvQtmp(j, l)
                zvQ(j, l, istn, iQ) = zvQtmp(j, l)-zvQ(j, l, immc, iQ)
             enddo
          enddo
          ! fonction de courant meridienne pour la quantite Q
          do l = llm, 1, -1
             do j = 1, jjm
                psiQ(j, l, iQ) = psiQ(j, l + 1, iQ) + zvQ(j, l, itot, iQ)
             enddo
          enddo
       enddo

       ! fonction de courant pour la circulation meridienne moyenne
       psi = 0.
       do l = llm, 1, -1
          do j = 1, jjm
             psi(j, l) = psi(j, l + 1) + zv(j, l)
             zv(j, l) = zv(j, l)*zfactv(j, l)
          enddo
       enddo

       ! sorties proprement dites
       do iQ = 1, nQ
          do itr = 1, ntr
             call histwrite(fileid, znom(itr, iQ), itau, zvQ(:, :, itr, iQ))
          enddo
          call histwrite(fileid, 'psi'//nom(iQ), itau, psiQ(:, :llm, iQ))
       enddo

       call histwrite(fileid, 'masse', itau, zmasse)
       call histwrite(fileid, 'v', itau, zv)
       psi = psi*1.e-9
       call histwrite(fileid, 'psi', itau, psi(:, :llm))

       ! Intégrale verticale

       forall (iQ = 1: nQ, itr = 2: ntr) zavQ(:, itr, iQ) &
            = sum(zvQ(:, :, itr, iQ) * zmasse, dim=2) / cv_2d(1, :)

       do iQ = 1, nQ
          do itr = 2, ntr
             call histwrite(fileid, 'a'//znom(itr, iQ), itau, zavQ(:, itr, iQ))
          enddo
       enddo

       ! On doit pouvoir tracer systematiquement la fonction de courant.
       icum = 0
    endif writing_step

  end SUBROUTINE bilan_dyn

end module bilan_dyn_m
