module bilan_dyn_m

  IMPLICIT NONE

contains

  SUBROUTINE bilan_dyn(ps, masse, pk, pbaru, pbarv, teta, phi, ucov, vcov, trac)

    ! From LMDZ4/libf/dyn3d/bilan_dyn.F, version 1.5, 2005/03/16 10:12:17

    ! Sous-programme consacr\'e \`a des diagnostics dynamiques de
    ! base. De fa\c{}con g\'en\'erale, les moyennes des scalaires Q
    ! sont pond\'er\'ees par la masse. Les flux de masse sont, eux,
    ! simplement moyenn\'es.

    USE comgeom, ONLY: constang_2d, cu_2d, cv_2d
    use covcont_m, only: covcont
    USE dimensions, ONLY: iim, jjm, llm
    use enercin_m, only: enercin
    USE histwrite_m, ONLY: histwrite
    use init_dynzon_m, only: ncum, fileid, znom, ntr, nq, nom
    use massbar_m, only: massbar
    USE paramet_m, ONLY: iip1, jjp1
    use suphec_m, only: rcpd

    real, intent(in):: ps(:, :) ! (iip1, jjp1)

    real, intent(in):: masse(:, :, :) ! (iip1, jjp1, llm)
    ! mass in a grid cell, in kg

    real, intent(in):: pk(:, :, :) ! (iip1, jjp1, llm)

    ! Flux de masse :
    real, intent(in):: pbaru(:, :, :) ! (iip1, jjp1, llm) in kg s-1
    real, intent(in):: pbarv(:, :, :) ! (iip1, jjm, llm) in kg s-1

    real, intent(in):: teta(:, :, :) ! (iip1, jjp1, llm)

    real, intent(in):: phi(:, :, :) ! (iip1, jjp1, llm)
    ! geopotential at mid-layer, in m2 s-2

    real, intent(in):: ucov(:, :, :) ! (iip1, jjp1, llm)
    real, intent(in):: vcov(:, :, :) ! (iip1, jjm, llm)
    real, intent(in):: trac(:, :, :) ! (iim + 1, jjm + 1, llm)

    ! Local:

    integer:: icum  = 0
    integer:: itau = 0
    real qy, factv(jjm, llm)

    ! Variables dynamiques interm\'ediaires :
    REAL vcont(iip1, jjm, llm), ucont(iip1, jjp1, llm)
    REAL ang(iip1, jjp1, llm), unat(iip1, jjp1, llm)
    REAL massebx(iip1, jjp1, llm), masseby(iip1, jjm, llm)
    REAL ecin(iip1, jjp1, llm)

    ! Champ contenant les scalaires advect\'es :
    real Q(iip1, jjp1, llm, nQ)

    ! Champs cumul\'es :
    real, save, allocatable:: ps_cum(:, :) ! (iip1, jjp1)
    real, save, allocatable:: masse_cum(:, :, :) ! (iip1, jjp1, llm)
    real, save, allocatable:: flux_u_cum(:, :, :) ! (iip1, jjp1, llm)
    real, save, allocatable:: flux_v_cum(:, :, :) ! (iip1, jjm, llm)
    real, save, allocatable:: Q_cum(:, :, :, :) ! (iip1, jjp1, llm, nQ)
    real, save, allocatable:: flux_uQ_cum(:, :, :, :) ! (iip1, jjp1, llm, nQ)
    real, save, allocatable:: flux_vQ_cum(:, :, :, :) ! (iip1, jjm, llm, nQ)

    ! Champs de tansport en moyenne zonale :
    integer itr
    integer, parameter:: iave = 1, itot = 2, immc = 3, itrs = 4, istn = 5

    real vq(jjm, llm, ntr, nQ), vqtmp(jjm, llm)
    real avq(jjm, 2: ntr, nQ), psiQ(jjm, llm + 1, nQ)
    real zmasse(jjm, llm)
    real v(jjm, llm), psi(jjm, llm + 1)
    integer i, j, l, iQ

    !-----------------------------------------------------------------

    if (itau == 0) then
       allocate(ps_cum(iip1, jjp1))
       allocate(masse_cum(iip1, jjp1, llm))
       allocate(flux_u_cum(iip1, jjp1, llm))
       allocate(flux_v_cum(iip1, jjm, llm))
       allocate(Q_cum(iip1, jjp1, llm, nQ))
       allocate(flux_uQ_cum(iip1, jjp1, llm, nQ))
       allocate(flux_vQ_cum(iip1, jjm, llm, nQ))
    end if

    ! Calcul des champs dynamiques

    ! \'Energie cin\'etique
    CALL covcont(ucov, vcov, ucont, vcont)
    CALL enercin(vcov, ucov, vcont, ucont, ecin)

    ! moment cin\'etique
    forall (l = 1: llm)
       ang(:, :, l) = ucov(:, :, l) + constang_2d
       unat(:, :, l) = ucont(:, :, l) * cu_2d
    end forall

    Q(:, :, :, 1) = teta * pk / rcpd
    Q(:, :, :, 2) = phi
    Q(:, :, :, 3) = ecin
    Q(:, :, :, 4) = ang
    Q(:, :, :, 5) = unat
    Q(:, :, :, 6) = trac
    Q(:, :, :, 7) = 1.

    ! Cumul

    if (icum == 0) then
       ps_cum = 0.
       masse_cum = 0.
       flux_u_cum = 0.
       flux_v_cum = 0.
       Q_cum = 0.
       flux_vQ_cum = 0.
       flux_uQ_cum = 0.
    endif

    itau = itau + 1
    icum = icum + 1

    ! Accumulation des flux de masse horizontaux
    ps_cum = ps_cum + ps
    masse_cum = masse_cum + masse
    flux_u_cum = flux_u_cum + pbaru
    flux_v_cum = flux_v_cum + pbarv
    forall (iQ = 1: nQ) Q_cum(:, :, :, iQ) = Q_cum(:, :, :, iQ) &
         + Q(:, :, :, iQ) * masse

    ! Flux longitudinal
    forall (iQ = 1: nQ, i = 1: iim) flux_uQ_cum(i, :, :, iQ) &
         = flux_uQ_cum(i, :, :, iQ) &
         + pbaru(i, :, :) * 0.5 * (Q(i, :, :, iQ) + Q(i + 1, :, :, iQ))
    flux_uQ_cum(iip1, :, :, :) = flux_uQ_cum(1, :, :, :)

    ! Flux m\'eridien
    forall (iQ = 1: nQ, j = 1: jjm) flux_vQ_cum(:, j, :, iQ) &
         = flux_vQ_cum(:, j, :, iQ) &
         + pbarv(:, j, :) * 0.5 * (Q(:, j, :, iQ) + Q(:, j + 1, :, iQ))

    writing_step: if (icum == ncum) then
       ! Normalisation
       forall (iQ = 1: nQ) Q_cum(:, :, :, iQ) = Q_cum(:, :, :, iQ) / masse_cum
       ps_cum = ps_cum / ncum
       masse_cum = masse_cum / ncum
       flux_u_cum = flux_u_cum / ncum
       flux_v_cum = flux_v_cum / ncum
       flux_uQ_cum = flux_uQ_cum / ncum
       flux_vQ_cum = flux_vQ_cum / ncum

       ! Transport m\'eridien

       ! Cumul zonal des masses des mailles

       v = 0.
       zmasse = 0.
       call massbar(masse_cum, massebx, masseby)
       do l = 1, llm
          do j = 1, jjm
             do i = 1, iim
                zmasse(j, l) = zmasse(j, l) + masseby(i, j, l)
                v(j, l) = v(j, l) + flux_v_cum(i, j, l)
             enddo
             factv(j, l) = cv_2d(1, j) / zmasse(j, l)
          enddo
       enddo

       ! Transport dans le plan latitude-altitude

       vq = 0.
       psiQ = 0.
       do iQ = 1, nQ
          vqtmp = 0.
          do l = 1, llm
             do j = 1, jjm
                ! Calcul des moyennes zonales du transport total et de vqtmp
                do i = 1, iim
                   vq(j, l, itot, iQ) = vq(j, l, itot, iQ) &
                        + flux_vQ_cum(i, j, l, iQ)
                   qy =  0.5 * (Q_cum(i, j, l, iQ) * masse_cum(i, j, l) &
                        + Q_cum(i, j + 1, l, iQ) * masse_cum(i, j + 1, l))
                   vqtmp(j, l) = vqtmp(j, l) + flux_v_cum(i, j, l) * qy &
                        / (0.5 * (masse_cum(i, j, l) + masse_cum(i, j + 1, l)))
                   vq(j, l, iave, iQ) = vq(j, l, iave, iQ) + qy
                enddo
                ! Decomposition
                vq(j, l, iave, iQ) = vq(j, l, iave, iQ) / zmasse(j, l)
                vq(j, l, itot, iQ) = vq(j, l, itot, iQ) * factv(j, l)
                vqtmp(j, l) = vqtmp(j, l) * factv(j, l)
                vq(j, l, immc, iQ) = v(j, l) * vq(j, l, iave, iQ) * factv(j, l)
                vq(j, l, itrs, iQ) = vq(j, l, itot, iQ) - vqtmp(j, l)
                vq(j, l, istn, iQ) = vqtmp(j, l) - vq(j, l, immc, iQ)
             enddo
          enddo
          ! Fonction de courant m\'eridienne pour la quantit\'e Q
          do l = llm, 1, -1
             do j = 1, jjm
                psiQ(j, l, iQ) = psiQ(j, l + 1, iQ) + vq(j, l, itot, iQ)
             enddo
          enddo
       enddo

       ! Fonction de courant pour la circulation m\'eridienne moyenne
       psi = 0.
       do l = llm, 1, -1
          do j = 1, jjm
             psi(j, l) = psi(j, l + 1) + v(j, l)
             v(j, l) = v(j, l) * factv(j, l)
          enddo
       enddo

       ! Sorties proprement dites
       do iQ = 1, nQ
          do itr = 1, ntr
             call histwrite(fileid, znom(itr, iQ), itau, vq(:, :, itr, iQ))
          enddo
          call histwrite(fileid, 'psi' // nom(iQ), itau, psiQ(:, :llm, iQ))
       enddo

       call histwrite(fileid, 'masse', itau, zmasse)
       call histwrite(fileid, 'v', itau, v)
       psi = psi * 1e-9
       call histwrite(fileid, 'psi', itau, psi(:, :llm))

       ! Int\'egrale verticale

       forall (iQ = 1: nQ, itr = 2: ntr) avq(:, itr, iQ) &
            = sum(vq(:, :, itr, iQ) * zmasse, dim=2) / cv_2d(1, :)

       do iQ = 1, nQ
          do itr = 2, ntr
             call histwrite(fileid, 'a' // znom(itr, iQ), itau, avq(:, itr, iQ))
          enddo
       enddo

       icum = 0
    endif writing_step

  end SUBROUTINE bilan_dyn

end module bilan_dyn_m
