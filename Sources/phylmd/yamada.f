module yamada_m

  IMPLICIT NONE

contains

  SUBROUTINE yamada(ngrid, g, zlev, zlay, u, v, teta, q2, km, kn)

    ! From LMDZ4/libf/phylmd/yamada.F,v 1.1 2004/06/22 11:45:36

    USE dimens_m
    USE dimphy
    ! .......................................................................
    ! .......................................................................

    ! g  : g
    ! zlev : altitude a chaque niveau (interface inferieure de la couche
    ! de meme indice)
    ! zlay : altitude au centre de chaque couche
    ! u,v : vitesse au centre de chaque couche
    ! (en entree : la valeur au debut du pas de temps)
    ! teta : temperature potentielle au centre de chaque couche
    ! (en entree : la valeur au debut du pas de temps)
    ! q2 : $q^2$ au bas de chaque couche
    ! (en entree : la valeur au debut du pas de temps)
    ! (en sortie : la valeur a la fin du pas de temps)
    ! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
    ! couche)
    ! (en sortie : la valeur a la fin du pas de temps)
    ! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
    ! (en sortie : la valeur a la fin du pas de temps)

    ! .......................................................................
    REAL, INTENT (IN) :: g
    REAL zlev(klon, klev+1)
    REAL zlay(klon, klev)
    REAL u(klon, klev)
    REAL v(klon, klev)
    REAL teta(klon, klev)
    REAL q2(klon, klev+1)
    REAL km(klon, klev+1)
    REAL kn(klon, klev+1)
    INTEGER ngrid


    INTEGER nlay, nlev
    PARAMETER (nlay=klev)
    PARAMETER (nlev=klev+1)

    LOGICAL first
    SAVE first
    DATA first/.TRUE./


    INTEGER ig, k

    REAL ri, zrif, zalpha, zsm
    REAL rif(klon, klev+1), sm(klon, klev+1), alpha(klon, klev)

    REAL m2(klon, klev+1), dz(klon, klev+1), zq, n2(klon, klev+1)
    REAL l(klon, klev+1), l0(klon)

    REAL sq(klon), sqz(klon), zz(klon, klev+1)
    INTEGER iter

    REAL ric, rifc, b1, kap
    SAVE ric, rifc, b1, kap
    DATA ric, rifc, b1, kap/0.195, 0.191, 16.6, 0.3/

    IF (0==1 .AND. first) THEN
       DO ig = 1, 1000
          ri = (ig-800.)/500.
          IF (ri<ric) THEN
             zrif = frif(ri)
          ELSE
             zrif = rifc
          END IF
          IF (zrif<0.16) THEN
             zalpha = falpha(zrif)
             zsm = fsm(zrif)
          ELSE
             zalpha = 1.12
             zsm = 0.085
          END IF
          PRINT *, ri, rif, zalpha, zsm
       END DO
       first = .FALSE.
    END IF

    ! Correction d'un bug sauvage a verifier.
    ! do k=2,nlev
    DO k = 2, nlay
       DO ig = 1, ngrid
          dz(ig, k) = zlay(ig, k) - zlay(ig, k-1)
          m2(ig, k) = ((u(ig,k)-u(ig,k-1))**2+(v(ig,k)-v(ig, &
               k-1))**2)/(dz(ig,k)*dz(ig,k))
          n2(ig, k) = g*2.*(teta(ig,k)-teta(ig,k-1))/(teta(ig,k-1)+teta(ig,k))/ &
               dz(ig, k)
          ri = n2(ig, k)/max(m2(ig,k), 1.E-10)
          IF (ri<ric) THEN
             rif(ig, k) = frif(ri)
          ELSE
             rif(ig, k) = rifc
          END IF
          IF (rif(ig,k)<0.16) THEN
             alpha(ig, k) = falpha(rif(ig,k))
             sm(ig, k) = fsm(rif(ig,k))
          ELSE
             alpha(ig, k) = 1.12
             sm(ig, k) = 0.085
          END IF
          zz(ig, k) = b1*m2(ig, k)*(1.-rif(ig,k))*sm(ig, k)
       END DO
    END DO

    ! iterration pour determiner la longueur de melange

    DO ig = 1, ngrid
       l0(ig) = 100.
    END DO
    DO k = 2, klev - 1
       DO ig = 1, ngrid
          l(ig, k) = l0(ig)*kap*zlev(ig, k)/(kap*zlev(ig,k)+l0(ig))
       END DO
    END DO

    DO iter = 1, 10
       DO ig = 1, ngrid
          sq(ig) = 1.E-10
          sqz(ig) = 1.E-10
       END DO
       DO k = 2, klev - 1
          DO ig = 1, ngrid
             q2(ig, k) = l(ig, k)**2*zz(ig, k)
             l(ig, k) = min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig, &
                  k)+l0(ig)), 0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.E-10)))
             zq = sqrt(q2(ig,k))
             sqz(ig) = sqz(ig) + zq*zlev(ig, k)*(zlay(ig,k)-zlay(ig,k-1))
             sq(ig) = sq(ig) + zq*(zlay(ig,k)-zlay(ig,k-1))
          END DO
       END DO
       DO ig = 1, ngrid
          l0(ig) = 0.2*sqz(ig)/sq(ig)
       END DO
       ! (abd 3 5 2)         print*,'ITER=',iter,'  L0=',l0

    END DO

    DO k = 2, klev
       DO ig = 1, ngrid
          l(ig, k) = min(l0(ig)*kap*zlev(ig,k)/(kap*zlev(ig, &
               k)+l0(ig)), 0.5*sqrt(q2(ig,k))/sqrt(max(n2(ig,k),1.E-10)))
          q2(ig, k) = l(ig, k)**2*zz(ig, k)
          km(ig, k) = l(ig, k)*sqrt(q2(ig,k))*sm(ig, k)
          kn(ig, k) = km(ig, k)*alpha(ig, k)
       END DO
    END DO

  contains

    REAL function frif(ri)
      real ri
      frif = 0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))
    end function frif

    REAL function falpha(ri)
      real ri
      falpha = 1.318*(0.2231-ri)/(0.2341-ri)
    end function falpha

    REAL function fsm(ri)
      real ri
      fsm = 1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))
    end function fsm

  END SUBROUTINE yamada

end module yamada_m
