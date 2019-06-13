module cltrac_m

  IMPLICIT NONE

contains

  SUBROUTINE cltrac(dtime, coef, t, tr, flux, paprs, pplay, delp, d_tr)

    ! From LMDZ4/libf/phylmd/cltrac.F, version 1.1.1.1 2004/05/19 12:53:07

    USE dimphy, only: klon, klev
    USE suphec_m, only: rd, rg

    ! Auteur : O. Boucher (LOA/LMD), date: 1996/11/27
    ! inspir\'e de clvent

    ! Objet: diffusion verticale de traceurs avec flux fix\'e \`a la
    ! surface ou flux du type c-drag

    REAL, INTENT(IN):: dtime ! intervalle du temps (en second)

    REAL, INTENT(IN):: coef(:, 2:) ! (klon, 2:klev)
    ! coefficient d'echange (m**2/s) l>1

    REAL, INTENT(IN):: t(klon, klev) ! temperature (K)
    REAL, INTENT(IN):: tr(klon, klev) ! la q. de traceurs
    REAL, INTENT(IN):: flux(klon) ! le flux de traceurs a la surface
    REAL, INTENT(IN):: paprs(klon, klev+1) ! pression a inter-couche (Pa)
    REAL, INTENT(IN):: pplay(klon, klev) ! pression au milieu de couche (Pa)
    REAL, INTENT(IN):: delp(klon, klev) ! epaisseur de couche (Pa)
    REAL, INTENT(out):: d_tr(klon, klev) ! le changement de tr

    ! Local:
    
    real tr0(klon)
    ! tr0 traceurs a la surface ou dans l'ocean (non active)

    INTEGER i, k
    REAL zx_ctr(klon, 2:klev)
    REAL zx_dtr(klon, 2:klev)
    REAL zx_buf(klon)
    REAL zx_coef(klon, klev)
    REAL local_tr(klon, klev)
    REAL zx_alf1(klon), zx_alf2(klon), zx_flux(klon)

    !-----------------------------------------------------------------------

    DO k = 1, klev
       DO i = 1, klon
          local_tr(i, k) = tr(i, k)
       END DO
    END DO

    DO i = 1, klon
       zx_alf1(i) = (paprs(i, 1)-pplay(i, 2))/(pplay(i, 1)-pplay(i, 2))
       zx_alf2(i) = 1.0 - zx_alf1(i)
       zx_flux(i) = -flux(i)*dtime*rg
       ! pour le moment le flux est prescrit
       ! zx_coef(1) vaut 0
       tr0(i) = 0.0
       zx_coef(i, 1) = 0.
    END DO

    DO k = 2, klev
       DO i = 1, klon
          zx_coef(i, k) = coef(i, k)*rg/(pplay(i, k-1)-pplay(i, k))* &
               (paprs(i, k)*2/(t(i, k)+t(i, k-1))/rd)**2
          zx_coef(i, k) = zx_coef(i, k)*dtime*rg
       END DO
    END DO

    DO i = 1, klon
       zx_buf(i) = delp(i, 1) + zx_coef(i, 1)*zx_alf1(i) + zx_coef(i, 2)
       zx_ctr(i, 2) = (local_tr(i, 1)*delp(i, 1)+zx_coef(i, 1)*tr0(i)-zx_flux(i))/ &
            zx_buf(i)
       zx_dtr(i, 2) = (zx_coef(i, 2)-zx_alf2(i)*zx_coef(i, 1))/zx_buf(i)
    END DO

    DO k = 3, klev
       DO i = 1, klon
          zx_buf(i) = delp(i, k-1) + zx_coef(i, k) + zx_coef(i, k-1)*(1.-zx_dtr(i &
               , k-1))
          zx_ctr(i, k) = (local_tr(i, k-1)*delp(i, k-1)+zx_coef(i, k-1)*zx_ctr(i, k-1 &
               ))/zx_buf(i)
          zx_dtr(i, k) = zx_coef(i, k)/zx_buf(i)
       END DO
    END DO
    DO i = 1, klon
       local_tr(i, klev) = (local_tr(i, klev)*delp(i, klev)+zx_coef(i, klev)*zx_ctr &
            (i, klev))/(delp(i, klev)+zx_coef(i, klev)-zx_coef(i, klev)*zx_dtr(i, klev))
    END DO
    DO k = klev - 1, 1, -1
       DO i = 1, klon
          local_tr(i, k) = zx_ctr(i, k+1) + zx_dtr(i, k+1)*local_tr(i, k+1)
       END DO
    END DO

    DO k = 1, klev
       DO i = 1, klon
          d_tr(i, k) = local_tr(i, k) - tr(i, k)
       END DO
    END DO

    RETURN
  END SUBROUTINE cltrac

end module cltrac_m
