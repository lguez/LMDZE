!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/screenc.F90,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      SUBROUTINE screenc(klon, knon, nsrf, zxli, &
                         speed, temp, q_zref, zref, &
                         ts, qsurf, rugos, psol, &
                         ustar, testar, qstar, okri, ri1, &
                         pref, delu, delte, delq)
      use YOMCST
      IMPLICIT NONE
!-----------------------------------------------------------------------
! 
! Objet : calcul "correcteur" des anomalies du vent, de la temperature 
!         potentielle et de l'humidite relative au niveau de reference zref et 
!         par rapport au 1er niveau (pour u) ou a la surface (pour theta et q) 
!         a partir des equations de Louis.
!
! Reference : Hess, Colman et McAvaney (1995)
!
! I. Musat, 01.07.2002
!-----------------------------------------------------------------------
!
! klon----input-I- dimension de la grille physique (= nb_pts_latitude X nb_pts_longitude)
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
! zxli----input-L- TRUE si calcul des cdrags selon Laurent Li
! speed---input-R- module du vent au 1er niveau du modele
! temp----input-R- temperature de l'air au 1er niveau du modele
! q_zref--input-R- humidite relative au 1er niveau du modele
! zref----input-R- altitude de reference
! ts------input-R- temperature de l'air a la surface
! qsurf---input-R- humidite relative a la surface
! rugos---input-R- rugosite
! psol----input-R- pression au sol
! ustar---input-R- facteur d'echelle pour le vent
! testar--input-R- facteur d'echelle pour la temperature potentielle
! qstar---input-R- facteur d'echelle pour l'humidite relative
! okri----input-L- TRUE si on veut tester le nb. Richardson entre la sfce 
!                  et zref par rapport au Ri entre la sfce et la 1ere couche
! ri1-----input-R- nb. Richardson entre la surface et la 1ere couche 
!
! pref----input-R- pression au niveau de reference
! delu----input-R- anomalie du vent par rapport au 1er niveau
! delte---input-R- anomalie de la temperature potentielle par rapport a la surface
! delq----input-R- anomalie de l'humidite relative par rapport a la surface
!
      INTEGER, intent(in) :: klon, knon, nsrf
      LOGICAL, intent(in) :: zxli, okri 
      REAL, dimension(klon), intent(in) :: speed, temp, q_zref
      REAL, intent(in) :: zref
      REAL, dimension(klon), intent(in) :: ts, qsurf, rugos, psol
      REAL, dimension(klon), intent(in) :: ustar, testar, qstar, ri1         
!
      REAL, dimension(klon), intent(out) :: pref, delu, delte, delq 
!-----------------------------------------------------------------------
!
! Variables locales  
      INTEGER :: i 
      REAL, dimension(klon) :: cdram, cdrah, cdran, zri1, gref
!
!------------------------------------------------------------------------- 
      DO i=1, knon
        gref(i) = zref*RG
      ENDDO 
!
! Richardson at reference level 
!
      CALL coefcdrag (klon, knon, nsrf, zxli, &
                    speed, temp, q_zref, gref, &
                    psol, ts, qsurf, rugos, &
                    okri, ri1, &
                    cdram, cdrah, cdran, zri1, &
                    pref)
!
      DO i = 1, knon
        delu(i) = ustar(i)/sqrt(cdram(i))
        delte(i)= (testar(i)* sqrt(cdram(i)))/ &
                   cdrah(i)
        delq(i)= (qstar(i)* sqrt(cdram(i)))/ &
                  cdrah(i)
      ENDDO 
!
      RETURN 
      END SUBROUTINE screenc
