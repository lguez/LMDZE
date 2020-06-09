module screencp_m

  implicit none

contains

  subroutine screencp(u_zref, temp, q_zref, speed, tpot, q1, ts1, qsurf, &
       rugos, lmon, ustar, testar, qstar, psol, pat1, nsrf, zref)

    use screenc_m, only: screenc
    use screenp_m, only: screenp
    USE suphec_m, ONLY: rkappa

    REAL, intent(out):: u_zref(:) ! (knon)
    ! anomalie du vent par rapport au 1er niveau

    REAL, intent(out):: temp(:), q_zref(:) ! (knon)
    REAL, intent(in):: speed(:) ! (knon) module du vent au 1er niveau du modele
    REAL, intent(in):: tpot(:) ! (knon) temperature potentielle
    REAL, intent(in):: q1(:) ! (knon) humidite relative au 1er niveau du modele
    REAL, intent(in):: ts1(:) ! (knon) temperature de l'air a la surface
    REAL, intent(in):: qsurf(:) ! (knon) humidit\'e relative \`a la surface
    REAL, intent(in):: rugos(:) ! (knon) rugosit\'e
    DOUBLE PRECISION, intent(in):: lmon(:) ! (knon) longueur de Monin-Obukov
    REAL, intent(in):: ustar(:) ! (knon) facteur d'\'echelle pour le vent

    REAL, intent(in):: testar(:) ! (knon)
    ! facteur d'echelle pour la temperature potentielle

    REAL, intent(in):: qstar(:) ! (knon)
    ! facteur d'echelle pour l'humidite relative

    REAL, intent(in):: psol(:) ! (knon) pression au sol
    REAL, intent(in):: pat1(:) ! (knon) pression au 1er niveau du modele
    INTEGER, intent(in):: nsrf ! indice pour le type de surface
    REAL, intent(in):: zref ! altitude de reference

    ! Local:

    INTEGER, parameter:: niter = 2 ! nombre iterations calcul "corrector"
    integer n
    REAL pref(size(speed)) ! (knon)

    REAL delte(size(speed)) ! (knon)
    ! anomalie de la temperature potentielle par rapport a la surface

    REAL delq(size(speed)) ! (knon)
    ! anomalie de l'humidite relative par rapport a la surface

    !---------------------------------------------------------------------

    ! First aproximation of variables at zref  
    CALL screenp(speed, tpot, q1, ts1, qsurf, rugos, lmon, ustar, testar, &
         qstar, zref, u_zref, delte, delq)
    q_zref = max(qsurf, 0.) + delq
    temp = (ts1 + delte) * (psol / pat1)**(- RKAPPA)

    ! Iteration of the variables at the reference level zref:
    ! corrector calculation ; see Hess & McAvaney, 1995
    DO n = 1, niter
       CALL screenc(nsrf, temp, q_zref, zref, ts1, qsurf, rugos, psol, ustar, &
            testar, qstar, pref, u_zref, delte, delq)
       q_zref = delq + max(qsurf, 0.)
       temp = (delte + ts1) * (psol / pref)**(- RKAPPA)
    ENDDO

  end subroutine screencp

end module screencp_m
