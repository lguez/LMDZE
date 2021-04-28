module dqthermcell_m

  implicit none

contains

  subroutine dqthermcell(ngrid, nlay, ptimestep, fm, entr, masse, q, dq, qa)

    use dimensions
    use dimphy

    ! Calcul du transport verticale dans la couche limite en presence
    ! de "thermiques" explicitement representes
    ! calcul du dq/dt une fois qu'on connait les ascendances

    integer ngrid, nlay

    real ptimestep
    real, intent(in):: masse(ngrid, nlay)
    real fm(ngrid, nlay+1)
    real entr(ngrid, nlay)
    real q(ngrid, nlay)
    real dq(ngrid, nlay)

    real qa(klon, klev), detr(klon, klev), wqd(klon, klev+1)

    integer ig, k

    ! calcul du detrainement

    do k=1, nlay
       do ig=1, ngrid
          detr(ig, k)=fm(ig, k)-fm(ig, k+1)+entr(ig, k)
       enddo
    enddo

    ! calcul de la valeur dans les ascendances
    do ig=1, ngrid
       qa(ig, 1)=q(ig, 1)
    enddo

    do k=2, nlay
       do ig=1, ngrid
          if ((fm(ig, k+1)+detr(ig, k))*ptimestep > &
               1.e-5*masse(ig, k)) then
             qa(ig, k)=(fm(ig, k)*qa(ig, k-1)+entr(ig, k)*q(ig, k)) &
                  /(fm(ig, k+1)+detr(ig, k))
          else
             qa(ig, k)=q(ig, k)
          endif
       enddo
    enddo

    do k=2, nlay
       do ig=1, ngrid
          ! wqd(ig, k)=fm(ig, k)*0.5*(q(ig, k-1)+q(ig, k))
          wqd(ig, k)=fm(ig, k)*q(ig, k)
       enddo
    enddo
    do ig=1, ngrid
       wqd(ig, 1)=0.
       wqd(ig, nlay+1)=0.
    enddo

    do k=1, nlay
       do ig=1, ngrid
          dq(ig, k)=(detr(ig, k)*qa(ig, k)-entr(ig, k)*q(ig, k) &
               -wqd(ig, k)+wqd(ig, k+1)) &
               /masse(ig, k)
       enddo
    enddo

  end subroutine dqthermcell

end module dqthermcell_m
