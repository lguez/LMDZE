SUBROUTINE dynredem0(fichnom,iday_end,phis,nq)

  ! From dyn3d/dynredem.F, v 1.2 2004/06/22 11:45:30

  ! Ecriture du fichier de redémarrage au format NetCDF (initialisation)

  USE IOIPSL, only: ymds2ju, ju2ymds
  use dimens_m, only: iim, jjm, llm
  use paramet_m, only: ip1jmp1, iip1, jjp1, llmp1
  use comconst, only: rad, cpp, daysec, dtvr, kappa, g, omeg
  use comvert, only: pa, bp, ap, nivsigs, preff, presnivs, nivsig
  use logic
  use comgeom
  use serre
  use temps, only: annee_ref, day_ref, itaufin, itau_dyn
  use ener
  use advtrac_m, only: tname, ttext

  IMPLICIT NONE

  include "netcdf.inc"

  !   Arguments:
  INTEGER, intent(in):: iday_end
  REAL, intent(in):: phis(ip1jmp1)
  CHARACTER(len=*) fichnom
  INTEGER nq

  !   Local:
  INTEGER iq,l
  INTEGER length
  PARAMETER (length = 100)
  REAL tab_cntrl(length) ! tableau des parametres du run
  INTEGER ierr

  !   Variables locales pour NetCDF:
  !
  INTEGER dims2(2), dims3(3), dims4(4)
  INTEGER idim_index
  INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
  INTEGER idim_s, idim_sig
  INTEGER idim_tim
  INTEGER nid,nvarid

  REAL zjulian,hours
  INTEGER yyears0,jjour0, mmois0
  character(len=30) unites

  !-----------------------------------------------------------------------

  print *, "Call sequence information: dynredem0"

  call ymds2ju(annee_ref, 1, iday_end, 0.0, zjulian)
  call ju2ymds(zjulian, yyears0, mmois0, jjour0, hours)

  DO l=1,length
     tab_cntrl(l) = 0.
  ENDDO
  tab_cntrl(1)  = REAL(iim)
  tab_cntrl(2)  = REAL(jjm)
  tab_cntrl(3)  = REAL(llm)
  tab_cntrl(4)  = REAL(day_ref)
  tab_cntrl(5)  = REAL(annee_ref)
  tab_cntrl(6)  = rad
  tab_cntrl(7)  = omeg
  tab_cntrl(8)  = g
  tab_cntrl(9)  = cpp
  tab_cntrl(10) = kappa
  tab_cntrl(11) = daysec
  tab_cntrl(12) = dtvr
  tab_cntrl(13) = etot0
  tab_cntrl(14) = ptot0
  tab_cntrl(15) = ztot0
  tab_cntrl(16) = stot0
  tab_cntrl(17) = ang0
  tab_cntrl(18) = pa
  tab_cntrl(19) = preff
  !
  !    .....    parametres  pour le zoom      ......   

  tab_cntrl(20)  = clon
  tab_cntrl(21)  = clat
  tab_cntrl(22)  = grossismx
  tab_cntrl(23)  = grossismy
  !
  IF ( fxyhypb )   THEN
     tab_cntrl(24) = 1.
     tab_cntrl(25) = dzoomx
     tab_cntrl(26) = dzoomy
     tab_cntrl(27) = 0.
     tab_cntrl(28) = taux
     tab_cntrl(29) = tauy
  ELSE
     tab_cntrl(24) = 0.
     tab_cntrl(25) = dzoomx
     tab_cntrl(26) = dzoomy
     tab_cntrl(27) = 0.
     tab_cntrl(28) = 0.
     tab_cntrl(29) = 0.
     IF( ysinus )  tab_cntrl(27) = 1.
  ENDIF

  tab_cntrl(30) = REAL(iday_end)
  tab_cntrl(31) = REAL(itau_dyn + itaufin)
  !
  !    .........................................................
  !
  ! Creation du fichier:
  !
  ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
  IF (ierr.NE.NF_NOERR) THEN
     WRITE(6,*)" Pb d ouverture du fichier "//fichnom
     WRITE(6,*)' ierr = ', ierr
     stop 1
  ENDIF
  !
  ! Preciser quelques attributs globaux:
  !
  ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 27, &
       "Fichier demmarage dynamique")
  !
  ! Definir les dimensions du fichiers:
  !
  ierr = NF_DEF_DIM (nid, "index", length, idim_index)
  ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
  ierr = NF_DEF_DIM (nid, "rlatu", jjp1, idim_rlatu)
  ierr = NF_DEF_DIM (nid, "rlonv", iip1, idim_rlonv)
  ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
  ierr = NF_DEF_DIM (nid, "sigs", llm, idim_s)
  ierr = NF_DEF_DIM (nid, "sig", llmp1, idim_sig)
  ierr = NF_DEF_DIM (nid, "temps", NF_UNLIMITED, idim_tim)
  !
  ierr = NF_ENDDEF(nid) ! sortir du mode de definition
  !
  ! Definir et enregistrer certains champs invariants:
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"controle",NF_FLOAT,1,idim_index,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, &
       "Parametres de controle")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23, &
       "Longitudes des points U")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, &
       "Latitudes des points U")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23, &
       "Longitudes des points V")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, &
       "Latitudes des points V")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"nivsigs",NF_FLOAT,1,idim_s,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 28, &
       "Numero naturel des couches s")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsigs)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"nivsig",NF_FLOAT,1,idim_sig,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 32, &
       "Numero naturel des couches sigma")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsig)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_sig,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26, &
       "Coefficient A pour hybride")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_sig,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26, &
       "Coefficient B pour hybride")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
  !
  ierr = NF_REDEF (nid)
  ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_s,nvarid)
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
  !
  ! Coefficients de passage cov. <-> contra. <--> naturel
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonu
  dims2(2) = idim_rlatu
  ierr = NF_DEF_VAR (nid,"cu",NF_FLOAT,2,dims2,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29, &
       "Coefficient de passage pour U")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonv
  dims2(2) = idim_rlatv
  ierr = NF_DEF_VAR (nid,"cv",NF_FLOAT,2,dims2,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29, &
       "Coefficient de passage pour V")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
  !
  ! Aire de chaque maille:
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonv
  dims2(2) = idim_rlatu
  ierr = NF_DEF_VAR (nid,"aire",NF_FLOAT,2,dims2,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22, &
       "Aires de chaque maille")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
  !
  ! Geopentiel au sol:
  !
  ierr = NF_REDEF (nid)
  dims2(1) = idim_rlonv
  dims2(2) = idim_rlatu
  ierr = NF_DEF_VAR (nid,"phisinit",NF_FLOAT,2,dims2,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19, &
       "Geopotentiel au sol")
  ierr = NF_ENDDEF(nid)
  ierr = NF_PUT_VAR_REAL (nid,nvarid,phis)
  !
  ! Definir les variables pour pouvoir les enregistrer plus tard:
  !
  ierr = NF_REDEF (nid) ! entrer dans le mode de definition
  !
  ierr = NF_DEF_VAR (nid,"temps",NF_FLOAT,1,idim_tim,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19, &
       "Temps de simulation")
  write(unites,200)yyears0,mmois0,jjour0
200 format('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", 30, &
       unites)

  !
  dims4(1) = idim_rlonu
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
  ierr = NF_DEF_VAR (nid,"ucov",NF_FLOAT,4,dims4,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9, &
       "Vitesse U")
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatv
  dims4(3) = idim_s
  dims4(4) = idim_tim
  ierr = NF_DEF_VAR (nid,"vcov",NF_FLOAT,4,dims4,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9, &
       "Vitesse V")
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
  ierr = NF_DEF_VAR (nid,"teta",NF_FLOAT,4,dims4,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 11, &
       "Temperature")
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
  IF(nq.GE.1) THEN
     DO iq=1,nq
        ierr = NF_DEF_VAR (nid,tname(iq),NF_FLOAT,4,dims4,nvarid)
        ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,ttext(iq))
     ENDDO
  ENDIF
  !
  dims4(1) = idim_rlonv
  dims4(2) = idim_rlatu
  dims4(3) = idim_s
  dims4(4) = idim_tim
  ierr = NF_DEF_VAR (nid,"masse",NF_FLOAT,4,dims4,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12, &
       "C est quoi ?")
  !
  dims3(1) = idim_rlonv
  dims3(2) = idim_rlatu
  dims3(3) = idim_tim
  ierr = NF_DEF_VAR (nid,"ps",NF_FLOAT,3,dims3,nvarid)
  ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 15, &
       "Pression au sol")
  !
  ierr = NF_ENDDEF(nid) ! sortir du mode de definition
  ierr = NF_CLOSE(nid) ! fermer le fichier

  PRINT*,'iim,jjm,llm,iday_end',iim,jjm,llm,iday_end
  PRINT*,'rad,omeg,g,cpp,kappa', &
       rad,omeg,g,cpp,kappa

END SUBROUTINE dynredem0
