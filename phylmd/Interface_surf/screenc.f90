module screenc_m

  IMPLICIT NONE

contains

  SUBROUTINE screenc(nsrf, temp, q_zref, zref, ts, qsurf, rugos, psol, ustar, &
       testar, qstar, pref, u_zref, delte, delq)

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

    INTEGER, intent(in):: nsrf ! indice pour le type de surface; voir indicesol

    REAL, intent(in):: temp(:) ! (knon)
    ! temperature de l'air au 1er niveau du modele

    REAL, intent(in):: q_zref(:) ! (knon)
    ! humidite relative au 1er niveau du modele

    REAL, intent(in):: zref ! altitude de reference
    REAL, intent(in):: ts(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidite relative a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosit\'e
    REAL, intent(in):: psol(:) ! (knon) pression au sol
    REAL, intent(in):: ustar(:) ! (knon) facteur d'\'echelle pour le vent

    REAL, intent(in):: testar(:) ! (knon)
    ! facteur d'echelle pour la temperature potentielle

    REAL, intent(in):: qstar(:) ! (knon)
    ! facteur d'echelle pour l'humidite relative

    REAL, intent(out):: pref(:) ! (knon) pression au niveau de reference

    REAL, intent(inout):: u_zref(:) ! (knon)
    ! anomalie du vent par rapport au 1er niveau

    REAL, intent(out):: delte(:) ! (knon)
    ! anomalie de la temperature potentielle par rapport a la surface

    REAL, intent(out):: delq(:) ! (knon)
    ! anomalie de l'humidite relative par rapport a la surface

    ! Local:
    REAL, dimension(size(temp)):: cdram, cdrah, gref ! (knon)

    !------------------------------------------------------------------------- 

    gref = zref * RG

    ! Richardson at reference level
    CALL cdrag(nsrf, u_zref, temp, q_zref, gref, psol, ts, qsurf, rugos, &
         cdram, cdrah, pref)

    u_zref = ustar / sqrt(cdram)
    delte = testar * sqrt(cdram) / cdrah
    delq = qstar * sqrt(cdram) / cdrah

  END SUBROUTINE screenc

end module screenc_m
