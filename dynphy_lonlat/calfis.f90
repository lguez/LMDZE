module calfis_m

  IMPLICIT NONE

contains

  SUBROUTINE calfis(ucov, vcov, teta, q, p3d, pk, phis, phi, w, dufi, dvfi, &
       dtetafi, dqfi, dayvrai, time, lafin)

    ! From dyn3d/calfis.F, version 1.3, 2005/05/25 13:10:09
    ! Authors: P. Le Van, F. Hourdin

    ! 1. R\'earrangement des tableaux et transformation des variables
    ! dynamiques en variables physiques

    ! 2. Calcul des tendances physiques
    ! 3. Retransformation des tendances physiques en tendances dynamiques

    ! Remarques:

    ! - Les vents sont donn\'es dans la physique par leurs composantes 
    ! naturelles.

    ! - La variable thermodynamique de la physique est une variable
    ! intensive : T.
    ! Pour la dynamique on prend T * (preff / p)**kappa

    ! - Les deux seules variables d\'ependant de la g\'eom\'etrie
    ! n\'ecessaires pour la physique sont la latitude (pour le
    ! rayonnement) et l'aire de la maille (quand on veut int\'egrer une
    ! grandeur horizontalement).

    use comgeom, only: apoln, cu_2d, cv_2d, unsaire, apols
    use dimensions, only: iim, jjm, llm, nqmx
    use dimphy, only: klon
    use disvert_m, only: preff
    use dynetat0_m, only: rlonu, rlonv
    use grid_change, only: dyn_phy, gr_fi_dyn
    use nr_util, only: pi
    use physiq_m, only: physiq
    use suphec_m, only: rcpd, rkappa, rg

    REAL, intent(in):: ucov(:, :, :) ! (iim + 1, jjm + 1, llm) 
    ! covariant zonal velocity

    REAL, intent(in):: vcov(:, :, :) ! (iim + 1, jjm, llm) 
    !covariant meridional velocity 

    REAL, intent(in):: teta(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! potential temperature

    REAL, intent(in):: q(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)
    ! mass fractions of advected fields

    REAL, intent(in):: p3d(:, :, :) ! (iim + 1, jjm + 1, llm + 1) 
    ! pressure at layer interfaces, in Pa
    ! ("p3d(i, j, l)" is at longitude "rlonv(i)", latitude "rlatu(j)",
    ! for interface "l")

    REAL, intent(in):: pk(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! Exner = cp * (p / preff)**kappa 

    REAL, intent(in):: phis(:, :) ! (iim + 1, jjm + 1)
    REAL, intent(in):: phi(:, :, :) ! (iim + 1, jjm + 1, llm)

    REAL, intent(in):: w(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! vertical mass flux, in kg / s

    REAL, intent(out):: dufi(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! tendency for the covariant zonal velocity (m2 s-2)

    REAL, intent(out):: dvfi(:, :, :) ! (iim + 1, jjm, llm)
    ! tendency for the natural meridional velocity

    REAL, intent(out):: dtetafi(:, :, :) ! (iim + 1, jjm + 1, llm)
    ! tendency for the potential temperature

    REAL, intent(out):: dqfi(:, :, :, :) ! (iim + 1, jjm + 1, llm, nqmx)

    integer, intent(in):: dayvrai
    ! current day number, based at value 1 on January 1st of annee_ref

    REAL, intent(in):: time ! time of day, as a fraction of day length
    LOGICAL, intent(in):: lafin

    ! Local:
    INTEGER i, j, l, ig0, iq
    REAL paprs(klon, llm + 1) ! aux interfaces des couches 
    REAL play(klon, llm) ! aux milieux des couches 
    REAL pphi(klon, llm), pphis(klon)
    REAL u(klon, llm), v(klon, llm)
    real zvfi(iim + 1, jjm + 1, llm)
    REAL t(klon, llm) ! temperature, in K
    real qx(klon, llm, nqmx) ! mass fractions of advected fields
    REAL omega(klon, llm)
    REAL d_u(klon, llm), d_v(klon, llm) ! tendances physiques du vent (m s-2)
    REAL d_t(klon, llm), d_qx(klon, llm, nqmx)
    REAL z1(iim)
    REAL pksurcp(iim + 1, jjm + 1)

    !-----------------------------------------------------------------------

    !!print *, "Call sequence information: calfis"

    ! 40. Transformation des variables dynamiques en variables physiques :

    ! 42. Pression intercouches :
    forall (l = 1: llm + 1) paprs(:, l) = pack(p3d(:, :, l), dyn_phy)

    ! 43. Température et pression milieu couche
    DO l = 1, llm
       pksurcp = pk(:, :, l) / rcpd
       play(:, l) = pack(preff * pksurcp**(1./ rkappa), dyn_phy)
       t(:, l) = pack(teta(:, :, l) * pksurcp, dyn_phy)
    ENDDO

    ! 43.bis Traceurs :
    forall (iq = 1: nqmx, l = 1: llm) &
         qx(:, l, iq) = pack(q(:, :, l, iq), dyn_phy)

    ! Geopotentiel calcule par rapport a la surface locale :
    forall (l = 1 :llm) pphi(:, l) = pack(phi(:, :, l), dyn_phy)
    pphis = pack(phis, dyn_phy)
    forall (l = 1: llm) pphi(:, l) = pphi(:, l) - pphis

    ! Calcul de la vitesse verticale :
    forall (l = 1: llm)
       omega(1, l) = w(1, 1, l) * rg / apoln
       omega(2: klon - 1, l) &
            = pack(w(:iim, 2: jjm, l) * rg * unsaire(:iim, 2: jjm), .true.)
       omega(klon, l) = w(1, jjm + 1, l) * rg / apols
    END forall

    ! 45. champ u:

    DO l = 1, llm
       DO j = 2, jjm
          ig0 = 1 + (j - 2) * iim
          u(ig0 + 1, l) = 0.5 &
               * (ucov(iim, j, l) / cu_2d(iim, j) + ucov(1, j, l) / cu_2d(1, j))
          DO i = 2, iim
             u(ig0 + i, l) = 0.5 * (ucov(i - 1, j, l) / cu_2d(i - 1, j) &
                  + ucov(i, j, l) / cu_2d(i, j))
          end DO
       end DO
    end DO

    ! 46.champ v:

    forall (j = 2: jjm, l = 1: llm) zvfi(:iim, j, l) = 0.5 &
         * (vcov(:iim, j - 1, l) / cv_2d(:iim, j - 1) &
         + vcov(:iim, j, l) / cv_2d(:iim, j))
    zvfi(iim + 1, 2:jjm, :) = zvfi(1, 2:jjm, :)

    ! 47. champs de vents au p\^ole nord 
    ! U = 1 / pi * integrale [ v * cos(long) * d long ]
    ! V = 1 / pi * integrale [ v * sin(long) * d long ]

    DO l = 1, llm
       z1(1) = (rlonu(1) - rlonu(iim) + 2. * pi) * vcov(1, 1, l) / cv_2d(1, 1)
       DO i = 2, iim
          z1(i) = (rlonu(i) - rlonu(i - 1)) * vcov(i, 1, l) / cv_2d(i, 1)
       ENDDO

       u(1, l) = SUM(COS(rlonv(:iim)) * z1) / pi
       zvfi(:, 1, l) = SUM(SIN(rlonv(:iim)) * z1) / pi
    ENDDO

    ! 48. champs de vents au p\^ole sud:
    ! U = 1 / pi * integrale [ v * cos(long) * d long ]
    ! V = 1 / pi * integrale [ v * sin(long) * d long ]

    DO l = 1, llm
       z1(1) = (rlonu(1) - rlonu(iim) + 2. * pi) * vcov(1, jjm, l) &
            /cv_2d(1, jjm)
       DO i = 2, iim
          z1(i) = (rlonu(i) - rlonu(i - 1)) * vcov(i, jjm, l) / cv_2d(i, jjm)
       ENDDO

       u(klon, l) = SUM(COS(rlonv(:iim)) * z1) / pi
       zvfi(:, jjm + 1, l) = SUM(SIN(rlonv(:iim)) * z1) / pi
    ENDDO

    forall(l = 1: llm) v(:, l) = pack(zvfi(:, :, l), dyn_phy)

    CALL physiq(lafin, dayvrai, time, paprs, play, pphi, pphis, u, v, t, qx, &
         omega, d_u, d_v, d_t, d_qx)

    ! transformation des tendances physiques en tendances dynamiques:

    ! 62. enthalpie potentielle
    do l = 1, llm
       dtetafi(:, :, l) = rcpd * gr_fi_dyn(d_t(:, l)) / pk(:, :, l)
    end do

    ! 63. traceurs
    DO iq = 1, nqmx
       DO l = 1, llm
          DO i = 1, iim + 1
             dqfi(i, 1, l, iq) = d_qx(1, l, iq)
             dqfi(i, jjm + 1, l, iq) = d_qx(klon, l, iq)
          ENDDO
          DO j = 2, jjm
             ig0 = 1 + (j - 2) * iim
             DO i = 1, iim
                dqfi(i, j, l, iq) = d_qx(ig0 + i, l, iq)
             ENDDO
             dqfi(iim + 1, j, l, iq) = dqfi(1, j, l, iq)
          ENDDO
       ENDDO
    ENDDO

    ! 65. champ u:
    DO l = 1, llm
       DO i = 1, iim + 1
          dufi(i, 1, l) = 0.
          dufi(i, jjm + 1, l) = 0.
       ENDDO

       DO j = 2, jjm
          ig0 = 1 + (j - 2) * iim
          DO i = 1, iim - 1
             dufi(i, j, l) = 0.5 * (d_u(ig0 + i, l) + d_u(ig0 + i + 1, l)) &
                  * cu_2d(i, j)
          ENDDO
          dufi(iim, j, l) = 0.5 * (d_u(ig0 + 1, l) + d_u(ig0 + iim, l)) &
               * cu_2d(iim, j)
          dufi(iim + 1, j, l) = dufi(1, j, l)
       ENDDO
    ENDDO

    ! 67. champ v:

    DO l = 1, llm
       DO j = 2, jjm - 1
          ig0 = 1 + (j - 2) * iim
          DO i = 1, iim
             dvfi(i, j, l) = 0.5 * (d_v(ig0 + i, l) + d_v(ig0 + i + iim, l)) &
                  * cv_2d(i, j)
          ENDDO
          dvfi(iim + 1, j, l) = dvfi(1, j, l)
       ENDDO
    ENDDO

    ! 68. champ v pr\`es des p\^oles:
    ! v = U * cos(long) + V * SIN(long)

    DO l = 1, llm
       DO i = 1, iim
          dvfi(i, 1, l) = d_u(1, l) * COS(rlonv(i)) + d_v(1, l) * SIN(rlonv(i))
          dvfi(i, jjm, l) = d_u(klon, l) * COS(rlonv(i)) &
               + d_v(klon, l) * SIN(rlonv(i))
          dvfi(i, 1, l) = 0.5 * (dvfi(i, 1, l) + d_v(i + 1, l)) * cv_2d(i, 1)
          dvfi(i, jjm, l) = 0.5 &
               * (dvfi(i, jjm, l) + d_v(klon - iim - 1 + i, l)) * cv_2d(i, jjm)
       ENDDO

       dvfi(iim + 1, 1, l) = dvfi(1, 1, l)
       dvfi(iim + 1, jjm, l) = dvfi(1, jjm, l)
    ENDDO

  END SUBROUTINE calfis

end module calfis_m
