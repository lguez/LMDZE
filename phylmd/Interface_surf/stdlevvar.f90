module stdlevvar_m

 IMPLICIT NONE

contains

  SUBROUTINE stdlevvar(nsrf, u1, v1, t1, q1, z1, ts1, qsurf, rugos, psol, &
       pat1, t2m, q2m, wind10m, ustar)

    ! From LMDZ4/libf/phylmd/stdlevvar.F90, version 1.3, 2005 May 25th

    ! Objet : calcul de la temp\'erature et de l'humidit\'e relative
    ! \`a 2 m et du module du vent \`a 10 m \`a partir des relations
    ! de Dyer-Businger et des \'equations de Louis.

    ! Reference: Hess, Colman and McAvaney (1995) 

    ! Author: I. Musat, July 1st, 2002

    use cdrag_m, only: cdrag
    USE suphec_m, ONLY: rg, rkappa
    use screencp_m, only: screencp

    INTEGER, intent(in):: nsrf ! indice pour le type de surface
    REAL, intent(in):: u1(:) ! (knon) vent zonal au 1er niveau du modele
    REAL, intent(in):: v1(:) ! (knon) vent meridien au 1er niveau du modele

    REAL, intent(in):: t1(:) ! (knon) temperature de l'air au 1er
                             ! niveau du modele
    
    REAL, intent(in):: q1(:) ! (knon) humidite relative au 1er niveau du modele
    REAL, intent(in):: z1(:) ! (knon) geopotentiel au 1er niveau du modele
    REAL, intent(in):: ts1(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidit\'e relative \`a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosit\'e
    REAL, intent(in):: psol(:) ! (knon) pression au sol
    REAL, intent(in):: pat1(:) ! (knon) pression au 1er niveau du modele
    REAL, intent(out):: t2m(:) ! (knon) temperature de l'air a 2m
    REAL, intent(out):: q2m(:) ! (knon) humidite relative a 2m
    REAL, intent(out):: wind10m(:) ! (knon) norme du vent \`a 10m
    REAL, intent(out):: ustar(:) ! (knon) u*

    ! Local:

    REAL, PARAMETER:: RKAR = 0.4 ! constante de von Karman
    REAL tpot(size(u1)) ! (knon) temperature potentielle
    REAL cdram(size(u1)), cdrah(size(u1))

    DOUBLE PRECISION lmon(size(u1)) ! (knon)
    ! longueur de Monin-Obukhov selon Hess, Colman and McAvaney 

    REAL, dimension(size(u1)):: speed, testar, qstar, zdte ! (knon)

    !------------------------------------------------------------------------- 

    speed = SQRT(u1**2 + v1**2)
    CALL cdrag(nsrf, speed, t1, q1, z1, psol, ts1, qsurf, rugos, cdram, cdrah)

    ! Star variables:
    tpot = t1 * (psol / pat1)**RKAPPA
    ustar = sqrt(cdram) * speed
    zdte = tpot - ts1
    zdte = sign(max(abs(zdte), 1.e-10), zdte)
    testar = (cdrah * zdte * speed)/ustar
    qstar = (cdrah * (max(q1, 0.) - max(qsurf, 0.)) * speed) / ustar
    lmon = (ustar**2 * tpot) / (RKAR * RG * testar)

    call screencp(speed, tpot, q1, ts1, qsurf, rugos, lmon, ustar, testar, &
         qstar, psol, pat1, nsrf, zref = 2., temp = t2m, q_zref = q2m)
    call screencp(speed, tpot, q1, ts1, qsurf, rugos, lmon, ustar, testar, &
         qstar, psol, pat1, nsrf, zref = 10., u_zref = wind10m)

  END subroutine stdlevvar

end module stdlevvar_m
