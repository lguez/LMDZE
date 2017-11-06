module cltrac_m

  IMPLICIT NONE

contains

  SUBROUTINE cltrac(dtime, coef, t, tr, flux, paprs, pplay, delp, d_tr)

    ! From LMDZ4/libf/phylmd/cltrac.F,v 1.1.1.1 2004/05/19 12:53:07

    USE dimens_m
    USE dimphy
    USE suphec_m
    ! ======================================================================
    ! Auteur(s): O. Boucher (LOA/LMD) date: 19961127
    ! inspire de clvent
    ! Objet: diffusion verticale de traceurs avec flux fixe a la surface
    ! ou/et flux du type c-drag
    ! ======================================================================
    ! Arguments:
    ! dtime----input-R- intervalle du temps (en second)
    ! coef-----input-R- le coefficient d'echange (m**2/s) l>1
    ! tr-------input-R- la q. de traceurs
    ! flux-----input-R- le flux de traceurs a la surface
    ! paprs----input-R- pression a inter-couche (Pa)
    ! pplay----input-R- pression au milieu de couche (Pa)
    ! delp-----input-R- epaisseur de couche (Pa)
    ! cdrag----input-R- cdrag pour le flux de surface (non active)
    ! tr0------input-R- traceurs a la surface ou dans l'ocean (non active)
    ! d_tr-----output-R- le changement de tr
    ! flux_tr--output-R- flux de tr
    ! ======================================================================
    REAL, INTENT (IN) :: dtime
    REAL coef(:, 2:) ! (klon, 2:klev)
    REAL, INTENT (IN) :: t(klon, klev) ! temperature (K)
    REAL tr(klon, klev)
    REAL, INTENT (IN) :: paprs(klon, klev+1)
    REAL, INTENT (IN) :: pplay(klon, klev)
    REAL delp(klon, klev)
    REAL d_tr(klon, klev)
    REAL flux(klon), cdrag(klon), tr0(klon)
    ! REAL flux_tr(klon,klev)
    ! ======================================================================
    ! ======================================================================
    INTEGER i, k
    REAL zx_ctr(klon, 2:klev)
    REAL zx_dtr(klon, 2:klev)
    REAL zx_buf(klon)
    REAL zx_coef(klon, klev)
    REAL local_tr(klon, klev)
    REAL zx_alf1(klon), zx_alf2(klon), zx_flux(klon)
    ! ======================================================================
    DO k = 1, klev
       DO i = 1, klon
          local_tr(i, k) = tr(i, k)
       END DO
    END DO


    ! ======================================================================
    DO i = 1, klon
       zx_alf1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
       zx_alf2(i) = 1.0 - zx_alf1(i)
       zx_flux(i) = -flux(i)*dtime*rg
       ! --pour le moment le flux est prescrit
       ! --cdrag et zx_coef(1) vaut 0
       cdrag(i) = 0.0
       tr0(i) = 0.0
       zx_coef(i, 1) = cdrag(i)*dtime*rg
    END DO
    ! ======================================================================
    DO k = 2, klev
       DO i = 1, klon
          zx_coef(i, k) = coef(i, k)*rg/(pplay(i,k-1)-pplay(i,k))* &
               (paprs(i,k)*2/(t(i,k)+t(i,k-1))/rd)**2
          zx_coef(i, k) = zx_coef(i, k)*dtime*rg
       END DO
    END DO
    ! ======================================================================
    DO i = 1, klon
       zx_buf(i) = delp(i, 1) + zx_coef(i, 1)*zx_alf1(i) + zx_coef(i, 2)
       zx_ctr(i, 2) = (local_tr(i,1)*delp(i,1)+zx_coef(i,1)*tr0(i)-zx_flux(i))/ &
            zx_buf(i)
       zx_dtr(i, 2) = (zx_coef(i,2)-zx_alf2(i)*zx_coef(i,1))/zx_buf(i)
    END DO

    DO k = 3, klev
       DO i = 1, klon
          zx_buf(i) = delp(i, k-1) + zx_coef(i, k) + zx_coef(i, k-1)*(1.-zx_dtr(i &
               ,k-1))
          zx_ctr(i, k) = (local_tr(i,k-1)*delp(i,k-1)+zx_coef(i,k-1)*zx_ctr(i,k-1 &
               ))/zx_buf(i)
          zx_dtr(i, k) = zx_coef(i, k)/zx_buf(i)
       END DO
    END DO
    DO i = 1, klon
       local_tr(i, klev) = (local_tr(i,klev)*delp(i,klev)+zx_coef(i,klev)*zx_ctr &
            (i,klev))/(delp(i,klev)+zx_coef(i,klev)-zx_coef(i,klev)*zx_dtr(i,klev))
    END DO
    DO k = klev - 1, 1, -1
       DO i = 1, klon
          local_tr(i, k) = zx_ctr(i, k+1) + zx_dtr(i, k+1)*local_tr(i, k+1)
       END DO
    END DO
    ! ======================================================================
    ! == flux_tr est le flux de traceur (positif vers bas)
    ! DO i = 1, klon
    ! flux_tr(i,1) = zx_coef(i,1)/(RG*dtime)
    ! ENDDO
    ! DO k = 2, klev
    ! DO i = 1, klon
    ! flux_tr(i,k) = zx_coef(i,k)/(RG*dtime)
    ! .               * (local_tr(i,k)-local_tr(i,k-1))
    ! ENDDO
    ! ENDDO
    ! ======================================================================
    DO k = 1, klev
       DO i = 1, klon
          d_tr(i, k) = local_tr(i, k) - tr(i, k)
       END DO
    END DO

    RETURN
  END SUBROUTINE cltrac

end module cltrac_m
