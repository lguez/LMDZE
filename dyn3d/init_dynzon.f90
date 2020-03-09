module init_dynzon_m

  IMPLICIT NONE

  integer, parameter:: ntr = 5
  integer, parameter:: nQ = 7
  integer ncum, fileid
  character(len=10) znom(ntr, nQ)
  character(len=4), parameter:: nom(nQ) = (/'T   ', 'gz  ', 'K   ', 'ang ', &
       'u   ', 'ovap', 'un  '/)

contains

  SUBROUTINE init_dynzon

    ! From LMDZ4/libf/dyn3d/bilan_dyn.F, version 1.5 2005/03/16 10:12:17

    use comconst, only: dtvr
    USE conf_gcm_m, ONLY: day_step, iperiod, periodav
    USE dimensions, ONLY: jjm, llm
    USE disvert_m, ONLY: presnivs
    use dynetat0_m, only: rlatv, itau_dyn
    use dynetat0_chosen_m, only: day_ref, annee_ref
    USE histbeg_totreg_m, ONLY: histbeg_totreg
    USE histdef_m, ONLY: histdef
    USE histend_m, ONLY: histend
    USE histvert_m, ONLY: histvert
    USE nr_util, ONLY: pi
    USE ymds2ju_m, ONLY: ymds2ju

    ! Local:

    real dt_cum
    character(len=5), parameter:: unites(nQ) = (/'K    ', 'm2/s2', 'm2/s2', &
         'ang  ', 'm/s  ', 'kg/kg', 'un   '/)

    ! Champs de tansport en moyenne zonale
    integer itr
    character(len=26) noml(ntr, nQ)
    character(len=12) zunites(ntr, nQ)
    character(len=3), parameter:: ctrs(ntr) = (/'   ', 'TOT', 'MMC', 'TRS', &
         'STN'/)
    integer iQ

    ! Initialisation du fichier contenant les moyennes zonales.

    integer horiid, vertiid
    double precision julian
    real rlong(jjm), rlatg(jjm)

    !-----------------------------------------------------------------

    print *, "Call sequence information: init_dynzon"

    ! Initialisation des fichiers
    ! ncum est la frequence de stokage en pas de temps
    ncum = day_step / iperiod * periodav
    dt_cum = day_step * periodav * dtvr

    ! Initialisation du fichier contenant les moyennes zonales

    CALL ymds2ju(annee_ref, 1, day_ref, 0.0, julian)

    rlong = 0.
    rlatg = rlatv * 180. / pi

    call histbeg_totreg('dynzon', rlong(:1), rlatg, 1, 1, 1, jjm, itau_dyn, &
         julian, dt_cum, horiid, fileid)

    ! Appel \`a histvert pour la grille verticale

    call histvert(fileid, 'presnivs', 'Niveaux sigma', 'mb', presnivs, vertiid)

    ! Appels \`a histdef pour la d\'efinition des variables \`a sauvegarder
    do iQ = 1, nQ
       do itr = 1, ntr
          if (itr == 1) then
             znom(itr, iQ) = nom(iQ)
             noml(itr, iQ) = nom(iQ)
             zunites(itr, iQ) = unites(iQ)
          else
             znom(itr, iQ) = ctrs(itr) // 'v' // nom(iQ)
             noml(itr, iQ) = 'transport: v * ' // nom(iQ) // ' ' // ctrs(itr)
             zunites(itr, iQ) = 'm/s * ' // unites(iQ)
          endif
       enddo
    enddo

    ! D\'eclarations des champs avec dimension verticale
    do iQ = 1, nQ
       do itr = 1, ntr
          call histdef(fileid, znom(itr, iQ), noml(itr, iQ), &
               zunites(itr, iQ), 1, jjm, horiid, llm, 1, llm, vertiid, &
               'ave(X)', dt_cum, dt_cum)
       enddo
       ! D\'eclarations pour les fonctions de courant
       call histdef(fileid, 'psi' // nom(iQ), 'stream fn. ' // noml(2, iQ), &
            zunites(2, iQ), 1, jjm, horiid, llm, 1, llm, vertiid, &
            'ave(X)', dt_cum, dt_cum)
    enddo

    ! D\'eclarations pour les champs de transport d'air
    call histdef(fileid, 'masse', 'masse', 'kg', 1, jjm, horiid, llm, 1, &
         llm, vertiid, 'ave(X)', dt_cum, dt_cum)
    call histdef(fileid, 'v', 'v', 'm/s', 1, jjm, horiid, llm, 1, llm, &
         vertiid, 'ave(X)', dt_cum, dt_cum)
    ! D\'eclarations pour les fonctions de courant
    call histdef(fileid, 'psi', 'stream fn. MMC ', 'mega t/s', 1, jjm, &
         horiid, llm, 1, llm, vertiid, 'ave(X)', dt_cum, dt_cum)

    ! D\'eclaration des champs 1D de transport en latitude
    do iQ = 1, nQ
       do itr = 2, ntr
          call histdef(fileid, 'a' // znom(itr, iQ), noml(itr, iQ), &
               zunites(itr, iQ), 1, jjm, horiid, 1, 1, 1, -99, 'ave(X)', &
               dt_cum, dt_cum)
       enddo
    enddo

    CALL histend(fileid)

  end SUBROUTINE init_dynzon

end module init_dynzon_m
