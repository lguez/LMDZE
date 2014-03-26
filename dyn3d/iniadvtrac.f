module iniadvtrac_m

  ! From advtrac.h, version 1.1.1.1 2004/05/19 12:53:06

  ! iq = 1 pour l'eau vapeur
  ! iq = 2 pour l'eau liquide
  ! et �ventuellement iq = 3, ..., nqmx pour les autres traceurs

  use dimens_m, only: nqmx

  implicit none

  private nqmx

  INTEGER iadv(nqmx) ! indice du sch�ma d'advection pour l'eau et les traceurs
  ! 11 means Van-Leer scheme for hadv et monotonous PPM for vadv

  integer, parameter:: allowed_adv(10) = (/0, 1, 2, 10, 12, 13, 14, 16, 17, 18/)
  ! Allowed values for hadv and vadv:
  ! 1: schema transport type "humidite specifique LMD"
  ! 2: schema amont
  ! 10: schema Van-leer (retenu pour l'eau vapeur et liquide)
  ! 12: schema Frederic Hourdin I
  ! 13: schema Frederic Hourdin II
  ! 14: schema Van-leer + humidite specifique 
  ! 16: schema PPM Monotone(Collela & Woodward 1984)
  ! 17: schema PPM Semi Monotone (overshoots autoris�s)
  ! 18: schema PPM Positif Defini (overshoots undershoots autoris�s)
  ! Pour Van-Leer plus vapeur d'eau satur�e : iadv(1)=4

  INTEGER hadv(nqmx) ! indice sch�ma transport horizontal 
  INTEGER vadv(nqmx) ! indice sch�ma transport vertical 
  character(len=8) tnom(nqmx) ! nom court du traceur
  character(len=10) tname(nqmx) ! nom du traceur pour restart
  character(len=13) ttext(nqmx) ! nom long du traceur pour sorties

contains

  subroutine iniadvtrac

    ! From dyn3d/iniadvtrac.F, version 1.3 2005/04/13 08:58:34

    ! Authors: P. Le Van, L. Fairhead, F. Hourdin, F. Codron,
    ! F. Forget, M.-A. Filiberti

    use nr_util, only: assert
    use jumble, only: new_unit

    ! Variables local to the procedure:

    character(len=3) descrq(18)
    integer iq, iostat, nq_local, unit

    !-----------------------------------------------------------------------

    print *, "Call sequence information: iniadvtrac"

    ! Initializations:
    descrq(10)='VL1'
    descrq(11)='VLP'
    descrq(12)='FH1'
    descrq(13)='FH2'
    descrq(14)='VLH'
    descrq(16)='PPM'
    descrq(17)='PPS'
    descrq(18)='PPP'

    ! Choix du sch�ma pour l'advection dans fichier "traceur.def"
    call new_unit(unit)
    open(unit, file='traceur.def', status='old', action="read", &
         position="rewind", iostat=iostat)
    if (iostat == 0) then
       print *, 'Ouverture de "traceur.def" ok'
       read(unit, fmt=*) nq_local
       print *, 'nombre de traceurs ', nq_local
       call assert(nq_local == nqmx, "iniadvtrac nq_local")

       do iq=1, nqmx
          read(unit, fmt=*) hadv(iq), vadv(iq), tnom(iq)
          if (.not. any(hadv(iq) == allowed_adv) &
               .or. .not. any(vadv(iq) == allowed_adv)) then
             print *, "bad number for advection scheme"
             stop 1
          end if
       end do
       close(unit) 
    else
       print *, 'Probl�me � l''ouverture de "traceur.def"'
       print *, 'Attention : on prend des valeurs par d�faut.'
       call assert(nqmx == 4, "iniadvtrac nqmx")
       hadv(:4) = (/14, 10, 10, 10/)
       vadv(:4) = hadv(:4)
       tnom(1) = 'H2Ov'
       tnom(2) = 'H2Ol'
       tnom(3) = 'RN'
       tnom(4) = 'PB'
       do iq = 1, nqmx
          print *, hadv(iq), vadv(iq), tnom(iq)
       end do
    ENDIF

    tname = tnom

    ! � partir du nom court du traceur et du sch�ma d'advection, on
    ! d�termine le nom long :
    do iq = 1, nqmx
       if (hadv(iq) /= vadv(iq)) then
          if (hadv(iq) == 10 .and. vadv(iq) == 16) then
             iadv(iq) = 11
          else
             print *, "Bad combination for hozizontal and vertical schemes."
             stop 1
          endif
       else
          iadv(iq) = hadv(iq)
       endif

       IF (iadv(iq) == 0) THEN
          ttext(iq) = tnom(iq)
       ELSE
          ttext(iq)=trim(tnom(iq)) // descrq(iadv(iq))
       endif
    end do

  END subroutine iniadvtrac

end module iniadvtrac_m
