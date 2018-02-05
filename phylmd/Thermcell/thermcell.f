module thermcell_m

  IMPLICIT NONE

contains

  SUBROUTINE thermcell(ngrid, nlay, ptimestep, pplay, pplev, pphi, pu, pv, pt, &
       po, pduadj, pdvadj, pdtadj, pdoadj, fm0, entr0, r_aspect, l_mix, w2di, &
       tho)

    ! Calcul du transport vertical dans la couche limite en pr\'esence
    ! de "thermiques" explicitement repr\'esent\'es. R\'ecriture \`a partir
    ! d'un listing papier \`a Habas, le 14/02/00. Le thermique est
    ! suppos\'e homog\`ene et dissip\'e par m\'elange avec son
    ! environnement. La longueur "l_mix" contr\^ole l'efficacit\'e du
    ! m\'elange. Le calcul du transport des diff\'erentes esp\`eces se fait
    ! en prenant en compte :
    ! 1. un flux de masse montant
    ! 2. un flux de masse descendant
    ! 3. un entra\^inement
    ! 4. un d\'etra\^inement

    USE dimphy, ONLY : klev, klon
    USE suphec_m, ONLY : rd, rg, rkappa

    ! arguments:

    INTEGER ngrid, nlay, w2di
    real tho
    real ptimestep, l_mix, r_aspect
    REAL, intent(in):: pt(ngrid, nlay)
    real pdtadj(ngrid, nlay)
    REAL, intent(in):: pu(ngrid, nlay)
    real pduadj(ngrid, nlay)
    REAL, intent(in):: pv(ngrid, nlay)
    real pdvadj(ngrid, nlay)
    REAL po(ngrid, nlay), pdoadj(ngrid, nlay)
    REAL, intent(in):: pplay(ngrid, nlay)
    real, intent(in):: pplev(ngrid, nlay+1)
    real, intent(in):: pphi(ngrid, nlay)

    integer idetr
    save idetr
    data idetr/3/

    ! local:

    INTEGER ig, k, l, lmaxa(klon), lmix(klon)
    ! CR: on remplace lmax(klon, klev+1)
    INTEGER lmax(klon), lmin(klon), lentr(klon)
    real linter(klon)
    real zmix(klon), fracazmix(klon)

    real zmax(klon), zw, zw2(klon, klev+1), ztva(klon, klev)

    real zlev(klon, klev+1)
    REAL zh(klon, klev), zdhadj(klon, klev)
    REAL ztv(klon, klev)
    real zu(klon, klev), zv(klon, klev), zo(klon, klev)
    real zva(klon, klev)
    real zua(klon, klev)
    real zoa(klon, klev)

    real zha(klon, klev)
    real wa_moy(klon, klev+1)
    real fraca(klon, klev+1)
    real fracc(klon, klev+1)
    real zf, zf2
    real thetath2(klon, klev), wth2(klon, klev)
    common/comtherm/thetath2, wth2

    real rho(klon, klev), rhobarz(klon, klev+1), masse(klon, klev)
    real zpspsk(klon, klev)

    real wmax(klon), wmaxa(klon)
    real fracd(klon, klev+1)
    real xxx(klon, klev+1)
    real larg_cons(klon, klev+1)
    real larg_detr(klon, klev+1)
    real fm0(klon, klev+1), entr0(klon, klev), detr(klon, klev)
    real fm(klon, klev+1), entr(klon, klev)
    real fmc(klon, klev+1)

    !CR:nouvelles variables
    real f_star(klon, klev+1), entr_star(klon, klev)
    real entr_star_tot(klon), entr_star2(klon)
    real f(klon)
    real zlevinter(klon)

    EXTERNAL SCOPY

    !-----------------------------------------------------------------------

    ! initialisation:

    IF(ngrid.NE.klon) THEN
       PRINT *
       PRINT *, 'STOP dans convadj'
       PRINT *, 'ngrid =', ngrid
       PRINT *, 'klon =', klon
    ENDIF

    ! incrementation eventuelle de tendances precedentes:

    print *, '0 OK convect8'

    DO l=1, nlay
       DO ig=1, ngrid
          zpspsk(ig, l)=(pplay(ig, l)/pplev(ig, 1))**RKAPPA
          zh(ig, l)=pt(ig, l)/zpspsk(ig, l)
          zu(ig, l)=pu(ig, l)
          zv(ig, l)=pv(ig, l)
          zo(ig, l)=po(ig, l)
          ztv(ig, l)=zh(ig, l)*(1.+0.61*zo(ig, l))
       end DO
    end DO

    print *, '1 OK convect8'

    ! See notes, "thermcell.txt"
    ! Calcul des altitudes des couches

    do l=2, nlay
       do ig=1, ngrid
          zlev(ig, l)=0.5*(pphi(ig, l)+pphi(ig, l-1))/RG
       enddo
    enddo
    do ig=1, ngrid
       zlev(ig, 1)=0.
       zlev(ig, nlay+1)=(2.*pphi(ig, klev)-pphi(ig, klev-1))/RG
    enddo

    ! Calcul des densites

    do l=1, nlay
       do ig=1, ngrid
          rho(ig, l)=pplay(ig, l)/(zpspsk(ig, l)*RD*zh(ig, l))
       enddo
    enddo

    do l=2, nlay
       do ig=1, ngrid
          rhobarz(ig, l)=0.5*(rho(ig, l)+rho(ig, l-1))
       enddo
    enddo

    ! Calcul de w2, quarre de w a partir de la cape
    ! a partir de w2, on calcule wa, vitesse de l'ascendance

    ! ATTENTION: Dans cette version, pour cause d'economie de memoire,
    ! w2 est stoke dans wa

    ! ATTENTION: dans convect8, on n'utilise le calcule des wa
    ! independants par couches que pour calculer l'entrainement
    ! a la base et la hauteur max de l'ascendance.

    ! Indicages:
    ! l'ascendance provenant du niveau k traverse l'interface l avec
    ! une vitesse wa(k, l).
    ! See notes, "thermcell.txt".

    !CR: ponderation entrainement des couches instables
    !def des entr_star tels que entr=f*entr_star
    do l=1, klev
       do ig=1, ngrid
          entr_star(ig, l)=0.
       enddo
    enddo
    ! determination de la longueur de la couche d entrainement
    do ig=1, ngrid
       lentr(ig)=1
    enddo

    !on ne considere que les premieres couches instables
    do k=nlay-2, 1, -1
       do ig=1, ngrid
          if (ztv(ig, k).gt.ztv(ig, k+1).and. &
               ztv(ig, k+1).le.ztv(ig, k+2)) then
             lentr(ig)=k
          endif
       enddo
    enddo

    ! determination du lmin: couche d ou provient le thermique
    do ig=1, ngrid
       lmin(ig)=1
    enddo
    do ig=1, ngrid
       do l=nlay, 2, -1
          if (ztv(ig, l-1).gt.ztv(ig, l)) then
             lmin(ig)=l-1
          endif
       enddo
    enddo

    ! definition de l'entrainement des couches
    do l=1, klev-1
       do ig=1, ngrid
          if (ztv(ig, l).gt.ztv(ig, l+1).and. &
               l.ge.lmin(ig).and.l.le.lentr(ig)) then
             entr_star(ig, l)=(ztv(ig, l)-ztv(ig, l+1))* &
                  (zlev(ig, l+1)-zlev(ig, l))
          endif
       enddo
    enddo
    ! pas de thermique si couches 1->5 stables
    do ig=1, ngrid
       if (lmin(ig).gt.5) then
          do l=1, klev
             entr_star(ig, l)=0.
          enddo
       endif
    enddo
    ! calcul de l entrainement total
    do ig=1, ngrid
       entr_star_tot(ig)=0.
    enddo
    do ig=1, ngrid
       do k=1, klev
          entr_star_tot(ig)=entr_star_tot(ig)+entr_star(ig, k)
       enddo
    enddo

    print *, 'fin calcul entr_star'
    do k=1, klev
       do ig=1, ngrid
          ztva(ig, k)=ztv(ig, k)
       enddo
    enddo

    do k=1, klev+1
       do ig=1, ngrid
          zw2(ig, k)=0.
          fmc(ig, k)=0.

          f_star(ig, k)=0.

          larg_cons(ig, k)=0.
          larg_detr(ig, k)=0.
          wa_moy(ig, k)=0.
       enddo
    enddo

    do ig=1, ngrid
       linter(ig)=1.
       lmaxa(ig)=1
       lmix(ig)=1
       wmaxa(ig)=0.
    enddo

    do l=1, nlay-2
       do ig=1, ngrid
          if (ztv(ig, l).gt.ztv(ig, l+1) &
               .and.entr_star(ig, l).gt.1.e-10 &
               .and.zw2(ig, l).lt.1e-10) then
             f_star(ig, l+1)=entr_star(ig, l)
             !test:calcul de dteta
             zw2(ig, l+1)=2.*RG*(ztv(ig, l)-ztv(ig, l+1))/ztv(ig, l+1) &
                  *(zlev(ig, l+1)-zlev(ig, l)) &
                  *0.4*pphi(ig, l)/(pphi(ig, l+1)-pphi(ig, l))
             larg_detr(ig, l)=0.
          else if ((zw2(ig, l).ge.1e-10).and. &
               (f_star(ig, l)+entr_star(ig, l).gt.1.e-10)) then
             f_star(ig, l+1)=f_star(ig, l)+entr_star(ig, l)
             ztva(ig, l)=(f_star(ig, l)*ztva(ig, l-1)+entr_star(ig, l) &
                  *ztv(ig, l))/f_star(ig, l+1)
             zw2(ig, l+1)=zw2(ig, l)*(f_star(ig, l)/f_star(ig, l+1))**2+ &
                  2.*RG*(ztva(ig, l)-ztv(ig, l))/ztv(ig, l) &
                  *(zlev(ig, l+1)-zlev(ig, l))
          endif
          ! determination de zmax continu par interpolation lineaire
          if (zw2(ig, l+1).lt.0.) then
             if (abs(zw2(ig, l+1)-zw2(ig, l)).lt.1e-10) then
                print *, 'pb linter'
             endif
             linter(ig)=(l*(zw2(ig, l+1)-zw2(ig, l)) &
                  -zw2(ig, l))/(zw2(ig, l+1)-zw2(ig, l))
             zw2(ig, l+1)=0.
             lmaxa(ig)=l
          else
             if (zw2(ig, l+1).lt.0.) then
                print *, 'pb1 zw2<0'
             endif
             wa_moy(ig, l+1)=sqrt(zw2(ig, l+1))
          endif
          if (wa_moy(ig, l+1).gt.wmaxa(ig)) then
             ! lmix est le niveau de la couche ou w (wa_moy) est maximum
             lmix(ig)=l+1
             wmaxa(ig)=wa_moy(ig, l+1)
          endif
       enddo
    enddo
    print *, 'fin calcul zw2'

    ! Calcul de la couche correspondant a la hauteur du thermique
    do ig=1, ngrid
       lmax(ig)=lentr(ig)
    enddo
    do ig=1, ngrid
       do l=nlay, lentr(ig)+1, -1
          if (zw2(ig, l).le.1.e-10) then
             lmax(ig)=l-1
          endif
       enddo
    enddo
    ! pas de thermique si couches 1->5 stables
    do ig=1, ngrid
       if (lmin(ig).gt.5) then
          lmax(ig)=1
          lmin(ig)=1
       endif
    enddo

    ! Determination de zw2 max
    do ig=1, ngrid
       wmax(ig)=0.
    enddo

    do l=1, nlay
       do ig=1, ngrid
          if (l.le.lmax(ig)) then
             if (zw2(ig, l).lt.0.)then
                print *, 'pb2 zw2<0'
             endif
             zw2(ig, l)=sqrt(zw2(ig, l))
             wmax(ig)=max(wmax(ig), zw2(ig, l))
          else
             zw2(ig, l)=0.
          endif
       enddo
    enddo

    ! Longueur caracteristique correspondant a la hauteur des thermiques.
    do ig=1, ngrid
       zmax(ig)=0.
       zlevinter(ig)=zlev(ig, 1)
    enddo
    do ig=1, ngrid
       ! calcul de zlevinter
       zlevinter(ig)=(zlev(ig, lmax(ig)+1)-zlev(ig, lmax(ig)))* &
            linter(ig)+zlev(ig, lmax(ig))-lmax(ig)*(zlev(ig, lmax(ig)+1) &
            -zlev(ig, lmax(ig)))
       zmax(ig)=max(zmax(ig), zlevinter(ig)-zlev(ig, lmin(ig)))
    enddo

    print *, 'avant fermeture'
    ! Fermeture, determination de f
    do ig=1, ngrid
       entr_star2(ig)=0.
    enddo
    do ig=1, ngrid
       if (entr_star_tot(ig).LT.1.e-10) then
          f(ig)=0.
       else
          do k=lmin(ig), lentr(ig)
             entr_star2(ig)=entr_star2(ig)+entr_star(ig, k)**2 &
                  /(rho(ig, k)*(zlev(ig, k+1)-zlev(ig, k)))
          enddo
          ! Nouvelle fermeture
          f(ig)=wmax(ig)/(max(500., zmax(ig))*r_aspect &
               *entr_star2(ig))*entr_star_tot(ig)
       endif
    enddo
    print *, 'apres fermeture'

    ! Calcul de l'entrainement
    do k=1, klev
       do ig=1, ngrid
          entr(ig, k)=f(ig)*entr_star(ig, k)
       enddo
    enddo
    ! Calcul des flux
    do ig=1, ngrid
       do l=1, lmax(ig)-1
          fmc(ig, l+1)=fmc(ig, l)+entr(ig, l)
       enddo
    enddo

    ! determination de l'indice du debut de la mixed layer ou w decroit

    ! calcul de la largeur de chaque ascendance dans le cas conservatif.
    ! dans ce cas simple, on suppose que la largeur de l'ascendance provenant
    ! d'une couche est \'egale \`a la hauteur de la couche alimentante.
    ! La vitesse maximale dans l'ascendance est aussi prise comme estimation
    ! de la vitesse d'entrainement horizontal dans la couche alimentante.

    do l=2, nlay
       do ig=1, ngrid
          if (l.le.lmaxa(ig)) then
             zw=max(wa_moy(ig, l), 1.e-10)
             larg_cons(ig, l)=zmax(ig)*r_aspect &
                  *fmc(ig, l)/(rhobarz(ig, l)*zw)
          endif
       enddo
    enddo

    do l=2, nlay
       do ig=1, ngrid
          if (l.le.lmaxa(ig)) then
             if ((l_mix*zlev(ig, l)).lt.0.)then
                print *, 'pb l_mix*zlev<0'
             endif
             larg_detr(ig, l)=sqrt(l_mix*zlev(ig, l))
          endif
       enddo
    enddo

    ! calcul de la fraction de la maille concern\'ee par l'ascendance en tenant
    ! compte de l'epluchage du thermique.

    !CR def de zmix continu (profil parabolique des vitesses)
    do ig=1, ngrid
       if (lmix(ig).gt.1.) then
          if (((zw2(ig, lmix(ig)-1)-zw2(ig, lmix(ig))) &
               *((zlev(ig, lmix(ig)))-(zlev(ig, lmix(ig)+1))) &
               -(zw2(ig, lmix(ig))-zw2(ig, lmix(ig)+1)) &
               *((zlev(ig, lmix(ig)-1))-(zlev(ig, lmix(ig))))).gt.1e-10) &
               then

             zmix(ig)=((zw2(ig, lmix(ig)-1)-zw2(ig, lmix(ig))) &
                  *((zlev(ig, lmix(ig)))**2-(zlev(ig, lmix(ig)+1))**2) &
                  -(zw2(ig, lmix(ig))-zw2(ig, lmix(ig)+1)) &
                  *((zlev(ig, lmix(ig)-1))**2-(zlev(ig, lmix(ig)))**2)) &
                  /(2.*((zw2(ig, lmix(ig)-1)-zw2(ig, lmix(ig))) &
                  *((zlev(ig, lmix(ig)))-(zlev(ig, lmix(ig)+1))) &
                  -(zw2(ig, lmix(ig))-zw2(ig, lmix(ig)+1)) &
                  *((zlev(ig, lmix(ig)-1))-(zlev(ig, lmix(ig))))))
          else
             zmix(ig)=zlev(ig, lmix(ig))
             print *, 'pb zmix'
          endif
       else
          zmix(ig)=0.
       endif

       if ((zmax(ig)-zmix(ig)).lt.0.) then
          zmix(ig)=0.99*zmax(ig)
       endif
    enddo

    ! calcul du nouveau lmix correspondant
    do ig=1, ngrid
       do l=1, klev
          if (zmix(ig).ge.zlev(ig, l).and. &
               zmix(ig).lt.zlev(ig, l+1)) then
             lmix(ig)=l
          endif
       enddo
    enddo

    do l=2, nlay
       do ig=1, ngrid
          if(larg_cons(ig, l).gt.1.) then
             fraca(ig, l)=(larg_cons(ig, l)-larg_detr(ig, l)) &
                  /(r_aspect*zmax(ig))
             fraca(ig, l)=max(fraca(ig, l), 0.)
             fraca(ig, l)=min(fraca(ig, l), 0.5)
             fracd(ig, l)=1.-fraca(ig, l)
             fracc(ig, l)=larg_cons(ig, l)/(r_aspect*zmax(ig))
          else
             fraca(ig, l)=0.
             fracc(ig, l)=0.
             fracd(ig, l)=1.
          endif
       enddo
    enddo
    !CR: calcul de fracazmix
    do ig=1, ngrid
       fracazmix(ig)=(fraca(ig, lmix(ig)+1)-fraca(ig, lmix(ig)))/ &
            (zlev(ig, lmix(ig)+1)-zlev(ig, lmix(ig)))*zmix(ig) &
            +fraca(ig, lmix(ig))-zlev(ig, lmix(ig))*(fraca(ig, lmix(ig)+1) &
            -fraca(ig, lmix(ig)))/(zlev(ig, lmix(ig)+1)-zlev(ig, lmix(ig)))
    enddo

    do l=2, nlay
       do ig=1, ngrid
          if(larg_cons(ig, l).gt.1.) then
             if (l.gt.lmix(ig)) then
                if (zmax(ig)-zmix(ig).lt.1.e-10) then
                   xxx(ig, l)=(lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
                else
                   xxx(ig, l)=(zmax(ig)-zlev(ig, l))/(zmax(ig)-zmix(ig))
                endif
                if (idetr.eq.0) then
                   fraca(ig, l)=fracazmix(ig)
                else if (idetr.eq.1) then
                   fraca(ig, l)=fracazmix(ig)*xxx(ig, l)
                else if (idetr.eq.2) then
                   fraca(ig, l)=fracazmix(ig)*(1.-(1.-xxx(ig, l))**2)
                else
                   fraca(ig, l)=fracazmix(ig)*xxx(ig, l)**2
                endif
                fraca(ig, l)=max(fraca(ig, l), 0.)
                fraca(ig, l)=min(fraca(ig, l), 0.5)
                fracd(ig, l)=1.-fraca(ig, l)
                fracc(ig, l)=larg_cons(ig, l)/(r_aspect*zmax(ig))
             endif
          endif
       enddo
    enddo

    print *, 'fin calcul fraca'

    ! Calcul de fracd, wd
    ! somme wa - wd = 0

    do ig=1, ngrid
       fm(ig, 1)=0.
       fm(ig, nlay+1)=0.
    enddo

    do l=2, nlay
       do ig=1, ngrid
          fm(ig, l)=fraca(ig, l)*wa_moy(ig, l)*rhobarz(ig, l)
          if (entr(ig, l-1).lt.1e-10.and.fm(ig, l).gt.fm(ig, l-1) &
               .and.l.gt.lmix(ig)) then
             fm(ig, l)=fm(ig, l-1)
          endif
       enddo
       do ig=1, ngrid
          if(fracd(ig, l).lt.0.1) then
             stop'fracd trop petit'
          endif
       enddo
    enddo

    do l=1, nlay
       do ig=1, ngrid
          masse(ig, l)=(pplev(ig, l)-pplev(ig, l+1))/RG
       enddo
    enddo

    print *, '12 OK convect8'

    ! calcul du transport vertical

    !CR:redefinition du entr
    do l=1, nlay
       do ig=1, ngrid
          detr(ig, l)=fm(ig, l)+entr(ig, l)-fm(ig, l+1)
          if (detr(ig, l).lt.0.) then
             entr(ig, l)=entr(ig, l)-detr(ig, l)
             detr(ig, l)=0.
          endif
       enddo
    enddo

    if (w2di.eq.1) then
       fm0=fm0+ptimestep*(fm-fm0)/tho
       entr0=entr0+ptimestep*(entr-entr0)/tho
    else
       fm0=fm
       entr0=entr
    endif

    if (1.eq.1) then
       call dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse &
            , zh, zdhadj, zha)
       call dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse &
            , zo, pdoadj, zoa)
    else
       call dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca &
            , zh, zdhadj, zha)
       call dqthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse, fraca &
            , zo, pdoadj, zoa)
    endif

    if (1.eq.0) then
       call dvthermcell2(ngrid, nlay, ptimestep, fm0, entr0, masse &
            , fraca, zmax &
            , zu, zv, pduadj, pdvadj, zua, zva)
    else
       call dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse &
            , zu, pduadj, zua)
       call dqthermcell(ngrid, nlay, ptimestep, fm0, entr0, masse &
            , zv, pdvadj, zva)
    endif

    do l=1, nlay
       do ig=1, ngrid
          zf=0.5*(fracc(ig, l)+fracc(ig, l+1))
          zf2=zf/(1.-zf)
          thetath2(ig, l)=zf2*(zha(ig, l)-zh(ig, l))**2
          wth2(ig, l)=zf2*(0.5*(wa_moy(ig, l)+wa_moy(ig, l+1)))**2
       enddo
    enddo

    do l=1, nlay
       do ig=1, ngrid
          pdtadj(ig, l)=zdhadj(ig, l)*zpspsk(ig, l)
       enddo
    enddo

  end SUBROUTINE thermcell

end module thermcell_m
