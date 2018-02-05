module screenc_m

  IMPLICIT NONE

contains

  SUBROUTINE screenc(klon, knon, nsrf, speed, temp, q_zref, zref, ts, &
       qsurf, rugos, psol, ustar, testar, qstar, pref, delu, delte, delq)
    
    ! From LMDZ4/libf/phylmd/screenc.F90, version 1.1.1.1, 2004/05/19 12:53:09

    ! Objet : calcul "correcteur" des anomalies du vent, de la
    ! temp\'erature potentielle et de l'humidit\'e relative au niveau
    ! de r\'ef\'erence zref et par rapport au 1er niveau (pour u) ou
    ! \`a la surface (pour theta et q) \`a partir des equations de
    ! Louis.

    ! Reference: Hess, Colman et McAvaney (1995)

    ! I. Musat, 01.07.2002

    use coefcdrag_m, only: coefcdrag
    use SUPHEC_M, only: RG

    INTEGER, intent(in):: klon
    ! klon----input-I- dimension de la grille physique (=
    ! nb_pts_latitude X nb_pts_longitude)
    INTEGER, intent(in):: knon
    ! knon----input-I- nombre de points pour un type de surface
    INTEGER, intent(in):: nsrf
    ! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
    REAL, dimension(klon), intent(in):: speed, temp, q_zref
    ! speed---input-R- module du vent au 1er niveau du modele
    ! temp----input-R- temperature de l'air au 1er niveau du modele
    ! q_zref--input-R- humidite relative au 1er niveau du modele
    REAL, intent(in):: zref
    ! zref----input-R- altitude de reference
    REAL, dimension(klon), intent(in):: ts, qsurf, rugos, psol
    ! ts------input-R- temperature de l'air a la surface
    ! qsurf---input-R- humidite relative a la surface
    ! rugos---input-R- rugosite
    ! psol----input-R- pression au sol
    REAL, intent(in):: ustar(:) ! (knon) facteur d'\'echelle pour le vent
    REAL, dimension(klon), intent(in):: testar
    ! testar--input-R- facteur d'echelle pour la temperature potentielle
    REAL, dimension(klon), intent(in):: qstar
    ! qstar---input-R- facteur d'echelle pour l'humidite relative

    REAL, dimension(klon), intent(out):: pref
    ! pref----input-R- pression au niveau de reference
    REAL, dimension(klon), intent(out):: delu
    ! delu----input-R- anomalie du vent par rapport au 1er niveau
    REAL, dimension(klon), intent(out):: delte
    ! delte---input-R- anomalie de la temperature potentielle par
    ! rapport a la surface
    REAL, dimension(klon), intent(out):: delq
    ! delq----input-R- anomalie de l'humidite relative par rapport a la surface

    ! Local:
    INTEGER i
    REAL, dimension(klon):: cdram, cdrah, cdran, zri1, gref

    !------------------------------------------------------------------------- 

    DO i=1, knon
       gref(i) = zref*RG
    ENDDO

    ! Richardson at reference level

    CALL coefcdrag(nsrf, speed(:knon), temp(:knon), q_zref(:knon), &
         gref(:knon), psol(:knon), ts, qsurf, rugos, cdram, cdrah, cdran, &
         zri1, pref)

    DO i = 1, knon
       delu(i) = ustar(i) / sqrt(cdram(i))
       delte(i) = testar(i) * sqrt(cdram(i)) / cdrah(i)
       delq(i) = qstar(i) * sqrt(cdram(i)) / cdrah(i)
    ENDDO

  END SUBROUTINE screenc

end module screenc_m
