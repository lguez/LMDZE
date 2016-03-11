module yamada4_m

  IMPLICIT NONE

  private
  public yamada4
  real, parameter:: kap = 0.4

contains

  SUBROUTINE yamada4(ngrid, dt, g, zlev, zlay, u, v, teta, cd, q2, km, kn, kq, &
       ustar, iflag_pbl)

    ! From LMDZ4/libf/phylmd/yamada4.F, version 1.1 2004/06/22 11:45:36

    use nr_util, only: assert
    USE dimphy, ONLY: klev

    integer, intent(in):: ngrid
    REAL, intent(in):: dt ! pas de temps
    real, intent(in):: g

    REAL zlev(ngrid, klev+1)
    ! altitude \`a chaque niveau (interface inf\'erieure de la couche de
    ! m\^eme indice)

    REAL zlay(ngrid, klev) ! altitude au centre de chaque couche

    REAL u(ngrid, klev), v(ngrid, klev)
    ! vitesse au centre de chaque couche (en entr\'ee : la valeur au
    ! d\'ebut du pas de temps)

    REAL, intent(in):: teta(ngrid, klev)
    ! temp\'erature potentielle au centre de chaque couche (en entr\'ee :
    ! la valeur au d\'ebut du pas de temps)

    REAL, intent(in):: cd(:) ! (ngrid) cdrag, valeur au d\'ebut du pas de temps

    REAL, intent(inout):: q2(ngrid, klev+1)
    ! $q^2$ au bas de chaque couche 
    ! En entr\'ee : la valeur au d\'ebut du pas de temps ; en sortie : la
    ! valeur \`a la fin du pas de temps.

    REAL km(ngrid, klev+1)
    ! diffusivit\'e turbulente de quantit\'e de mouvement (au bas de
    ! chaque couche) (en sortie : la valeur \`a la fin du pas de temps)

    REAL kn(ngrid, klev+1)
    ! diffusivit\'e turbulente des scalaires (au bas de chaque couche)
    ! (en sortie : la valeur \`a la fin du pas de temps)

    REAL kq(ngrid, klev+1)
    real ustar(ngrid)

    integer, intent(in):: iflag_pbl
    ! iflag_pbl doit valoir entre 6 et 9
    ! l = 6, on prend syst\'ematiquement une longueur d'\'equilibre
    ! iflag_pbl = 6 : MY 2.0
    ! iflag_pbl = 7 : MY 2.0.Fournier
    ! iflag_pbl = 8 : MY 2.5
    ! iflag_pbl = 9 : un test ?

    ! Local:

    real kmin, qmin, pblhmin(ngrid), coriol(ngrid)
    real qpre
    REAL unsdz(ngrid, klev)
    REAL unsdzdec(ngrid, klev+1)
    REAL kmpre(ngrid, klev+1), tmp2
    REAL mpre(ngrid, klev+1)
    real delta(ngrid, klev+1)
    real aa(ngrid, klev+1), aa1
    integer, PARAMETER:: nlev = klev+1
    logical:: first = .true.
    integer:: ipas = 0
    integer ig, k
    real ri
    real rif(ngrid, klev+1), sm(ngrid, klev+1), alpha(ngrid, klev)
    real m2(ngrid, klev+1), dz(ngrid, klev+1), zq, n2(ngrid, klev+1)
    real dtetadz(ngrid, klev+1)
    real m2cstat, mcstat, kmcstat
    real l(ngrid, klev+1)
    real l0(ngrid)
    real sq(ngrid), sqz(ngrid), zz(ngrid, klev+1)
    integer iter
    real:: ric = 0.195, rifc = 0.191, b1 = 16.6, kap = 0.4

    !-----------------------------------------------------------------------

    call assert(iflag_pbl >= 6 .and. iflag_pbl <= 9, "yamada4")

    ipas = ipas+1

    ! les increments verticaux
    DO ig = 1, ngrid
       ! alerte: zlev n'est pas declare a nlev
       zlev(ig, nlev) = zlay(ig, klev) +(zlay(ig, klev) - zlev(ig, nlev-1))
    ENDDO

    DO k = 1, klev
       DO ig = 1, ngrid
          unsdz(ig, k) = 1.E+0/(zlev(ig, k+1)-zlev(ig, k))
       ENDDO
    ENDDO

    DO ig = 1, ngrid
       unsdzdec(ig, 1) = 1.E+0/(zlay(ig, 1)-zlev(ig, 1))
    ENDDO

    DO k = 2, klev
       DO ig = 1, ngrid
          unsdzdec(ig, k) = 1.E+0/(zlay(ig, k)-zlay(ig, k-1))
       ENDDO
    ENDDO

    DO ig = 1, ngrid
       unsdzdec(ig, klev+1) = 1.E+0/(zlev(ig, klev+1)-zlay(ig, klev))
    ENDDO

    do k = 2, klev
       do ig = 1, ngrid
          dz(ig, k) = zlay(ig, k)-zlay(ig, k-1)
          m2(ig, k) = ((u(ig, k)-u(ig, k-1))**2+(v(ig, k)-v(ig, k-1))**2) &
               /(dz(ig, k)*dz(ig, k))
          dtetadz(ig, k) = (teta(ig, k)-teta(ig, k-1))/dz(ig, k)
          n2(ig, k) = g*2.*dtetadz(ig, k)/(teta(ig, k-1)+teta(ig, k))
          ri = n2(ig, k)/max(m2(ig, k), 1.e-10)
          if (ri.lt.ric) then
             rif(ig, k) = frif(ri)
          else
             rif(ig, k) = rifc
          endif
          if (rif(ig, k).lt.0.16) then
             alpha(ig, k) = falpha(rif(ig, k))
             sm(ig, k) = fsm(rif(ig, k))
          else
             alpha(ig, k) = 1.12
             sm(ig, k) = 0.085
          endif
          zz(ig, k) = b1*m2(ig, k)*(1.-rif(ig, k))*sm(ig, k)
       enddo
    enddo

    ! Au premier appel, on d\'etermine l et q2 de façon it\'erative.
    ! It\'eration pour d\'eterminer la longueur de m\'elange

    if (first .or. iflag_pbl == 6) then
       do ig = 1, ngrid
          l0(ig) = 10.
       enddo
       do k = 2, klev-1
          do ig = 1, ngrid
             l(ig, k) = l0(ig) * kap * zlev(ig, k) &
                  / (kap * zlev(ig, k) + l0(ig))
          enddo
       enddo

       do iter = 1, 10
          do ig = 1, ngrid
             sq(ig) = 1e-10
             sqz(ig) = 1e-10
          enddo
          do k = 2, klev-1
             do ig = 1, ngrid
                q2(ig, k) = l(ig, k)**2 * zz(ig, k)
                l(ig, k) = fl(zlev(ig, k), l0(ig), q2(ig, k), n2(ig, k))
                zq = sqrt(q2(ig, k))
                sqz(ig) = sqz(ig) + zq * zlev(ig, k) &
                     * (zlay(ig, k) - zlay(ig, k-1))
                sq(ig) = sq(ig) + zq * (zlay(ig, k) - zlay(ig, k-1))
             enddo
          enddo
          do ig = 1, ngrid
             l0(ig) = 0.2 * sqz(ig) / sq(ig)
          enddo
       enddo
    endif

    ! Calcul de la longueur de melange.

    ! Mise a jour de l0
    do ig = 1, ngrid
       sq(ig) = 1.e-10
       sqz(ig) = 1.e-10
    enddo
    do k = 2, klev-1
       do ig = 1, ngrid
          zq = sqrt(q2(ig, k))
          sqz(ig) = sqz(ig)+zq*zlev(ig, k)*(zlay(ig, k)-zlay(ig, k-1))
          sq(ig) = sq(ig)+zq*(zlay(ig, k)-zlay(ig, k-1))
       enddo
    enddo
    do ig = 1, ngrid
       l0(ig) = 0.2*sqz(ig)/sq(ig)
    enddo
    ! calcul de l(z)
    do k = 2, klev
       do ig = 1, ngrid
          l(ig, k) = fl(zlev(ig, k), l0(ig), q2(ig, k), n2(ig, k))
          if (first) then
             q2(ig, k) = l(ig, k)**2 * zz(ig, k)
          endif
       enddo
    enddo

    ! Yamada 2.0
    if (iflag_pbl == 6) then
       do k = 2, klev
          do ig = 1, ngrid
             q2(ig, k) = l(ig, k)**2 * zz(ig, k)
          enddo
       enddo
    else if (iflag_pbl == 7) then
       ! Yamada 2.Fournier

       ! Calcul de l, km, au pas precedent
       do k = 2, klev
          do ig = 1, ngrid
             delta(ig, k) = q2(ig, k) / (l(ig, k)**2 * sm(ig, k))
             kmpre(ig, k) = l(ig, k) * sqrt(q2(ig, k)) * sm(ig, k)
             mpre(ig, k) = sqrt(m2(ig, k))
          enddo
       enddo

       do k = 2, klev-1
          do ig = 1, ngrid
             m2cstat = max(alpha(ig, k)*n2(ig, k)+delta(ig, k)/b1, 1.e-12)
             mcstat = sqrt(m2cstat)

             ! puis on ecrit la valeur de q qui annule l'equation de m
             ! supposee en q3

             IF (k == 2) THEN
                kmcstat = 1.E+0 / mcstat &
                     *(unsdz(ig, k)*kmpre(ig, k+1) &
                     *mpre(ig, k+1) &
                     +unsdz(ig, k-1) &
                     *cd(ig) &
                     *(sqrt(u(ig, 3)**2+v(ig, 3)**2) &
                     -mcstat/unsdzdec(ig, k) &
                     -mpre(ig, k+1)/unsdzdec(ig, k+1))**2) &
                     /(unsdz(ig, k)+unsdz(ig, k-1))
             ELSE
                kmcstat = 1.E+0 / mcstat &
                     *(unsdz(ig, k)*kmpre(ig, k+1) &
                     *mpre(ig, k+1) &
                     +unsdz(ig, k-1)*kmpre(ig, k-1) &
                     *mpre(ig, k-1)) &
                     /(unsdz(ig, k)+unsdz(ig, k-1))
             ENDIF
             tmp2 = kmcstat / (sm(ig, k) / q2(ig, k)) /l(ig, k)
             q2(ig, k) = max(tmp2, 1.e-12)**(2./3.)
          enddo
       enddo
    else if (iflag_pbl >= 8) then
       ! Yamada 2.5 a la Didi

       ! Calcul de l, km, au pas precedent
       do k = 2, klev
          do ig = 1, ngrid
             delta(ig, k) = q2(ig, k)/(l(ig, k)**2*sm(ig, k))
             if (delta(ig, k).lt.1.e-20) then
                delta(ig, k) = 1.e-20
             endif
             km(ig, k) = l(ig, k)*sqrt(q2(ig, k))*sm(ig, k)
             aa1 = (m2(ig, k)*(1.-rif(ig, k))-delta(ig, k)/b1)
             aa(ig, k) = aa1*dt/(delta(ig, k)*l(ig, k))
             qpre = sqrt(q2(ig, k))
             if (iflag_pbl == 8) then
                if (aa(ig, k).gt.0.) then
                   q2(ig, k) = (qpre+aa(ig, k)*qpre*qpre)**2
                else
                   q2(ig, k) = (qpre/(1.-aa(ig, k)*qpre))**2
                endif
             else
                ! iflag_pbl = 9
                if (aa(ig, k)*qpre.gt.0.9) then
                   q2(ig, k) = (qpre*10.)**2
                else
                   q2(ig, k) = (qpre/(1.-aa(ig, k)*qpre))**2
                endif
             endif
             q2(ig, k) = min(max(q2(ig, k), 1.e-10), 1.e4)
          enddo
       enddo
    endif

    ! Calcul des coefficients de m\'elange
    do k = 2, klev
       do ig = 1, ngrid
          zq = sqrt(q2(ig, k))
          km(ig, k) = l(ig, k)*zq*sm(ig, k)
          kn(ig, k) = km(ig, k)*alpha(ig, k)
          kq(ig, k) = l(ig, k)*zq*0.2
       enddo
    enddo

    ! Traitement des cas noctrunes avec l'introduction d'une longueur
    ! minilale.

    ! Traitement particulier pour les cas tres stables.
    ! D'apres Holtslag Boville.

    do ig = 1, ngrid
       coriol(ig) = 1.e-4
       pblhmin(ig) = 0.07*ustar(ig)/max(abs(coriol(ig)), 2.546e-5)
    enddo

    print *, 'pblhmin ', pblhmin
    do k = 2, klev
       do ig = 1, ngrid
          if (teta(ig, 2).gt.teta(ig, 1)) then
             qmin = ustar(ig)*(max(1.-zlev(ig, k)/pblhmin(ig), 0.))**2
             kmin = kap*zlev(ig, k)*qmin
          else
             kmin = -1. ! kmin n'est utilise que pour les SL stables.
          endif
          if (kn(ig, k).lt.kmin.or.km(ig, k).lt.kmin) then
             kn(ig, k) = kmin
             km(ig, k) = kmin
             kq(ig, k) = kmin
             ! la longueur de melange est suposee etre l = kap z
             ! K = l q Sm d'ou q2 = (K/l Sm)**2
             q2(ig, k) = (qmin/sm(ig, k))**2
          endif
       enddo
    enddo

    first = .false.

  end SUBROUTINE yamada4

  !*******************************************************************

  real function frif(ri)

    real, intent(in):: ri

    frif = 0.6588*(ri+0.1776-sqrt(ri*ri-0.3221*ri+0.03156))

  end function frif

  !*******************************************************************

  real function falpha(ri)

    real, intent(in):: ri

    falpha = 1.318*(0.2231-ri)/(0.2341-ri)

  end function falpha

  !*******************************************************************

  real function fsm(ri)

    real, intent(in):: ri

    fsm = 1.96*(0.1912-ri)*(0.2341-ri)/((1.-ri)*(0.2231-ri))

  end function fsm

  !*******************************************************************

  real function fl(zzz, zl0, zq2, zn2)

    real, intent(in):: zzz, zl0, zq2, zn2

    fl = max(min(zl0 * kap * zzz / (kap * zzz + zl0), &
         0.5 * sqrt(zq2) / sqrt(max(zn2, 1e-10))), 1.)

  end function fl

end module yamada4_m
