module inter_bary_m

  implicit none

contains


  !******************************

  function inter_bary(yjdat, fdat, yjmod)

    ! From dyn3d/inter_bary.F, version 1.1.1.1 2004/05/19 12:53:06
    ! Authors: R. Sadourny, P. Le Van

    ! Interpolation barycentrique basée sur les aires.
    ! Version unidimensionnelle, en latitude.
    ! L'indice 1 correspond à l'interface maille 1 -- maille 2.

    use nr_util, only: assert


    REAL, intent(in):: yjdat(:)
    ! (angles, ordonnées des interfaces des mailles des données, in
    ! degrees, in increasing order)

    REAL, intent(in):: fdat(:) ! champ de données

    REAL, intent(in):: yjmod(:)
    ! (ordonnées des interfaces des mailles du modèle)
    ! (in degrees, in strictly increasing order)

    REAL inter_bary(size(yjmod)) ! champ du modèle

    ! Variables local to the procedure:

    REAL y0, dy, dym
    INTEGER jdat ! indice du champ de données
    integer jmod ! indice du champ du modèle

    !------------------------------------

    call assert(size(yjdat) == size(fdat), "inter_bary")

    ! Initialisation des variables
    inter_bary(:) = 0.
    y0    = -90.
    dym   = 0.
    jmod  = 1
    jdat  = 1

    do while (jmod <= size(yjmod))
       do while (yjmod(jmod) > yjdat(jdat))
          dy         = yjdat(jdat) - y0
          dym        = dym + dy
          inter_bary(jmod) = inter_bary(jmod) + dy * fdat(jdat)
          y0         = yjdat(jdat)
          jdat       = jdat + 1
       end do
       IF (yjmod(jmod) < yjdat(jdat)) THEN
          dy         = yjmod(jmod) - y0
          dym        = dym + dy
          inter_bary(jmod) = (inter_bary(jmod) + dy * fdat(jdat)) / dym
          y0         = yjmod(jmod)
          dym        = 0.
          jmod       = jmod + 1
       ELSE
          ! {yjmod(jmod) == yjdat(jdat)}
          dy         = yjmod(jmod) - y0
          dym        = dym + dy
          inter_bary(jmod) = (inter_bary(jmod) + dy * fdat(jdat)) / dym
          y0         = yjmod(jmod)
          dym        = 0.
          jmod       = jmod + 1
          jdat       = jdat + 1
       END IF
    end do
    ! Le test de fin suppose que l'interface 0 est commune aux deux
    ! grilles "yjdat" et "yjmod".

  END function inter_bary

end module inter_bary_m
