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

    use cdrag_m, only: cdrag
    use SUPHEC_M, only: RG

    INTEGER, intent(in):: klon
    ! klon----input-I- dimension de la grille physique (=
    ! nb_pts_latitude X nb_pts_longitude)
    INTEGER, intent(in):: knon
    ! knon----input-I- nombre de points pour un type de surface
    INTEGER, intent(in):: nsrf
    ! nsrf----input-I- indice pour le type de surface; voir indicesol.inc
    REAL, dimension(klon), intent(in):: speed
    ! speed---input-R- module du vent au 1er niveau du modele
    REAL, dimension(klon), intent(in):: temp
    ! temp----input-R- temperature de l'air au 1er niveau du modele
    REAL, dimension(klon), intent(in):: q_zref
    ! q_zref--input-R- humidite relative au 1er niveau du modele
    REAL, intent(in):: zref
    ! zref----input-R- altitude de reference
    REAL, intent(in):: ts(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidite relative a la surface
    REAL, dimension(klon), intent(in):: rugos
    ! rugos---input-R- rugosite
    REAL, intent(in):: psol(:) ! (knon) pression au sol
    REAL, intent(in):: ustar(:) ! (knon) facteur d'\'echelle pour le vent
    REAL, dimension(klon), intent(in):: testar
    ! testar--input-R- facteur d'echelle pour la temperature potentielle
    REAL, dimension(klon), intent(in):: qstar
    ! qstar---input-R- facteur d'echelle pour l'humidite relative

    REAL, intent(out):: pref(:) ! (knon) pression au niveau de reference
    REAL, dimension(klon), intent(out):: delu
    ! delu----input-R- anomalie du vent par rapport au 1er niveau
    REAL, dimension(klon), intent(out):: delte
    ! delte---input-R- anomalie de la temperature potentielle par
    ! rapport a la surface
    REAL, dimension(klon), intent(out):: delq
    ! delq----input-R- anomalie de l'humidite relative par rapport a la surface

    ! Local:
    INTEGER i
    REAL, dimension(knon):: cdram, cdrah, gref

    !------------------------------------------------------------------------- 

    DO i=1, knon
       gref(i) = zref*RG
    ENDDO

    ! Richardson at reference level

    CALL cdrag(nsrf, speed(:knon), temp(:knon), q_zref(:knon), gref, psol, ts, &
         qsurf, rugos(:knon), cdram, cdrah, pref)

    DO i = 1, knon
       delu(i) = ustar(i) / sqrt(cdram(i))
       delte(i) = testar(i) * sqrt(cdram(i)) / cdrah(i)
       delq(i) = qstar(i) * sqrt(cdram(i)) / cdrah(i)
    ENDDO

  END SUBROUTINE screenc

end module screenc_m
