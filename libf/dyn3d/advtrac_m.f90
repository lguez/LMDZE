module advtrac_m

  ! From advtrac.h, v 1.1.1.1 2004/05/19 12:53:06

  use dimens_m, only: nqmx

  implicit none

  private nqmx

  INTEGER iadv(nqmx) ! indice schema de transport 
  INTEGER hadv(nqmx) ! indice schema transport horizontal 
  INTEGER vadv(nqmx) ! indice schema transport vertical 
  INTEGER niadv(nqmx) ! equivalent dyn / physique
  character(len=8) tnom(nqmx) ! nom court du traceur
  character(len=10) tname(nqmx) ! nom du traceur pour restart
  character(len=13) ttext(nqmx) ! nom long du traceur pour sorties

contains

  subroutine iniadvtrac(nq)

    ! From dyn3d/iniadvtrac.F, version 1.3 2005/04/13 08:58:34

    ! Authors :  P. Le Van, L. Fairhead, F. Hourdin
    ! Modification spéciale traceur F. Forget 05/94
    ! Modification M.-A. Filiberti 02/02 lecture de "traceur.def"
    ! Modification de l'intégration de "q" (26/04/94)

    integer, intent(out):: nq

    ! Variables local to the procedure:

    character(len=3) descrq(30)
    character:: txts(3) = (/'x', 'y', 'z'/)
    character(len=2) txtp(9)
    character(len=13) str1, str2

    integer iq, iiq, iiiq, ierr, ii

    data txtp/'x', 'y', 'z', 'xx', 'xy', 'xz', 'yy', 'yz', 'zz'/

    !-----------------------------------------------------------------------

    print *, "Call sequence information: iniadvtrac"

    ! Initializations:
    descrq(14)='VLH'
    descrq(10)='VL1'
    descrq(11)='VLP'
    descrq(12)='FH1'
    descrq(13)='FH2'
    descrq(16)='PPM'
    descrq(17)='PPS'
    descrq(18)='PPP'
    descrq(20)='SLP'
    descrq(30)='PRA'

    !        Choix  des schemas d'advection pour l'eau et les traceurs

    ! iadv = 1    schema  transport type "humidite specifique LMD"
    ! iadv = 2    schema   amont
    ! iadv = 14    schema  Van-leer + humidite specifique 
    !                        Modif F.Codron
    ! iadv = 10   schema  Van-leer (retenu pour l'eau vapeur et liquide)
    ! iadv = 11 schema Van-Leer pour hadv et version PPM (Monotone)
    !           pour vadv
    ! iadv = 12   schema  Frederic Hourdin I
    ! iadv = 13   schema  Frederic Hourdin II
    ! iadv = 16   schema  PPM Monotone(Collela & Woodward 1984)
    ! iadv = 17   schema  PPM Semi Monotone (overshoots autorisés)
    ! iadv = 18   schema  PPM Positif Defini (overshoots undershoots autorisés)
    ! iadv = 20   schema  Slopes
    ! iadv = 30   schema  Prather

    !    Dans le tableau q(ij, l, iq) : iq = 1  pour l'eau vapeur
    !                                 iq = 2  pour l'eau liquide
    !    Et éventuellement            iq = 3, nqmx pour les autres traceurs

    !    iadv(1): choix pour l'eau vap. et  iadv(2) : choix pour l'eau liq.

    ! Choix du schema d'advection
    ! choix par defaut = van leer pour tous les traceurs
    do iq=1, nqmx
       iadv(iq)=10
       str1(1:1)='q'
       if (nqmx.le.99) then
          WRITE(str1(2:3), '(i2.2)') iq
       else
          WRITE(str1(2:4), '(i3.3)') iq
       endif
       tnom(iq)=str1
       tname(iq)=tnom(iq) 
       str2=tnom(iq) 
       ttext(iq) = trim(str2) // descrq(iadv(iq))
    end do

    ! Choix du schema pour l'advection dans fichier "traceur.def"

    open(unit=90, file='traceur.def', form='formatted', status='old', &
         iostat=ierr)
    if (ierr == 0) then
       print *, 'Ouverture de "traceur.def" ok'
       read(unit=90, fmt=*) nq
       print *, 'nombre de traceurs ', nq
       if (nq > nqmx) then
          print *, 'nombre de traceurs trop important'
          print *, 'verifier traceur.def'
          stop
       endif

       do iq=1, nq
          read(90, 999) hadv(iq), vadv(iq), tnom(iq)
       end do
       close(90)  
       PRINT *, 'lecture de traceur.def :'   
       do iq=1, nq
          write(*, *) hadv(iq), vadv(iq), tnom(iq)
       end do
    else
       print *, 'problème ouverture traceur.def'
       print *, 'ATTENTION on prend des valeurs par défaut'
       nq = 4
       hadv(1) = 14
       vadv(1) = 14
       tnom(1) = 'H2Ov'
       hadv(2) = 10
       vadv(2) = 10
       tnom(2) = 'H2Ol'
       hadv(3) = 10
       vadv(3) = 10
       tnom(3) = 'RN'
       hadv(4) = 10
       vadv(4) = 10
       tnom(4) = 'PB'
    ENDIF
    PRINT *, 'Valeur de traceur.def :'
    do iq=1, nq
       write(*, *) hadv(iq), vadv(iq), tnom(iq)
    end do

    ! À partir du nom court du traceur et du schéma d'advection, on
    ! détemine le nom long :
    iiq=0
    ii=0
    do iq=1, nq 
       iiq=iiq+1
       if (hadv(iq) /= vadv(iq)) then
          if (hadv(iq) == 10.and.vadv(iq) == 16) then
             iadv(iiq)=11
          else
             print *, 'le choix des schemas d''advection H et V'
             print *, 'est non disponible actuellement'
             stop 
          endif
       else
          iadv(iiq)=hadv(iq)
       endif
       ! verification nombre de traceurs
       if (iadv(iiq).lt.20) then
          ii=ii+1
       elseif (iadv(iiq) == 20) then
          ii=ii+4
       elseif (iadv(iiq) == 30) then
          ii=ii+10
       endif

       str1=tnom(iq)
       tname(iiq)=tnom(iq)
       IF (iadv(iiq) == 0) THEN
          ttext(iiq)=trim(str1)
       ELSE
          ttext(iiq)=trim(str1)//descrq(iadv(iiq))
       endif
       str2=ttext(iiq)
       !   schemas tenant compte des moments d'ordre superieur.
       if (iadv(iiq) == 20) then
          do iiiq=1, 3
             iiq=iiq+1
             iadv(iiq)=-20
             ttext(iiq)=trim(str2)//txts(iiiq)
             tname(iiq)=trim(str1)//txts(iiiq)
          enddo
       elseif (iadv(iiq) == 30) then
          do iiiq=1, 9
             iiq=iiq+1
             iadv(iiq)=-30
             ttext(iiq)=trim(str2)//txtp(iiiq)
             tname(iiq)=trim(str1)//txtp(iiiq)
          enddo
       endif
    end do
    if (ii /= nqmx) then
       print *, 'WARNING'
       print *, 'le nombre de traceurs et de moments eventuels'
       print *, 'est inferieur a nqmx '
    endif
    if (iiq > nqmx) then
       print *, 'le choix des schemas est incompatible avec '
       print *, 'la dimension nqmx (nombre de traceurs)'
       print *, 'verifier traceur.def ou la namelist INCA'
       print *, 'ou recompiler avec plus de traceurs'
       stop
    endif
    iiq=0
    do iq=1, nqmx
       if (iadv(iq).ge.0) then
          iiq=iiq+1
          niadv(iiq)=iq
       endif
    end do

999 format (i2, 1x, i2, 1x, a8)

  END subroutine iniadvtrac

end module advtrac_m
