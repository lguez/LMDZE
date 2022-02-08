module infotrac_init_m

  ! From advtrac.h, version 1.1.1.1 2004/05/19 12:53:06

  ! iq = 1 pour l'eau vapeur
  ! iq = 2 pour l'eau liquide
  ! et \'eventuellement iq = 3, ..., nqmx pour les autres traceurs

  use dimensions, only: nqmx

  implicit none

  private nqmx

  INTEGER, protected:: iadv(nqmx)
  ! indice du sch\'ema d'advection pour l'eau et les traceurs

  character(len = 10), protected:: tname(nqmx)
  ! nom du traceur pour fichiers restart et historiques

  character(len = 13), protected:: ttext(nqmx)
  ! nom long du traceur pour sorties

contains

  subroutine infotrac_init

    ! From dyn3d/iniadvtrac.F, version 1.3 2005/04/13 08:58:34

    ! Authors: P. Le Van, L. Fairhead, F. Hourdin, F. Codron,
    ! F. Forget, M.-A. Filiberti

    ! Initialisation des traceurs
    ! Choix du sch\'ema pour l'advection dans le fichier "traceur.def"

    use jumble, only: assert, new_unit

    ! Local:
    character(len = 3) descrq(0:14)
    integer iq, iostat, nq_local, unit

    integer, parameter:: allowed_adv(5) = (/0, 10, 12, 13, 14/)
    ! Allowed values for iadv:
    ! 10: schema Van-leer (retenu pour l'eau vapeur et liquide)
    ! 12: schema Frederic Hourdin I
    ! 13: schema Frederic Hourdin II
    ! 14: schema Van-leer + humidite specifique 

    !-----------------------------------------------------------------------

    print *, "Call sequence information: infotrac_init"

    ! Initializations:
    descrq(0) = ''
    descrq(10) = 'VL1'
    descrq(12) = 'FH1'
    descrq(13) = 'FH2'
    descrq(14) = 'VLH'

    ! Choix du sch\'ema pour l'advection dans fichier "traceur.def"
    call new_unit(unit)
    open(unit, file = 'traceur.def', status = 'old', action = "read", &
         position = "rewind", iostat = iostat)
    if (iostat == 0) then
       print *, 'Ouverture de "traceur.def" ok'
       read(unit, fmt = *) nq_local
       print *, 'nombre de traceurs ', nq_local
       call assert(nq_local == nqmx, "infotrac_init nq_local")

       do iq = 1, nqmx
          read(unit, fmt = *) iadv(iq), tname(iq)
          if (.not. any(iadv(iq) == allowed_adv)) then
             print *, "bad number for advection scheme"
             stop 1
          end if
       end do
       close(unit) 
    else
       print *, 'Could not open "traceur.def".'
       print *, 'Using default values.'
       call assert(nqmx == 4, "infotrac_init nqmx")
       iadv(:4) = (/14, 10, 10, 10/)
       tname(1) = 'H2Ov'
       tname(2) = 'H2Ol'
       tname(3) = 'RN'
       tname(4) = 'PB'
       do iq = 1, nqmx
          print *, iadv(iq), tname(iq)
       end do
    ENDIF

    ! \`A partir du nom court du traceur et du sch\'ema d'advection, on
    ! d\'etermine le nom long :
    do iq = 1, nqmx
       ttext(iq) = trim(tname(iq)) // descrq(iadv(iq))
    end do

  END subroutine infotrac_init

end module infotrac_init_m
