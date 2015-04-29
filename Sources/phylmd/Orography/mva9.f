module MVA9_m

  implicit none

contains

  SUBROUTINE MVA9(X)

    ! From dyn3d/grid_noro.F, version 1.1.1.1 2004/05/19 12:53:06

    ! Make a moving average over 9 gridpoints of the x fields.

    REAL, intent(inout):: X(:, :) ! (IMAR, JMAR)

    ! Local:
    integer imar, jmar
    real XF(size(x, 1), size(x, 2))
    real WEIGHTpb(-1:1, -1:1)
    real my_sum
    integer i, is, js, j

    !---------------------------------------------------------------

    imar = size(x, 1)
    jmar = size(x, 2)

    MY_SUM=0.
    DO IS=-1, 1
       DO JS=-1, 1
          WEIGHTpb(IS, JS) = 1. / REAL((1 + IS**2) * (1 + JS**2))
          MY_SUM=MY_SUM+WEIGHTpb(IS, JS)
       ENDDO
    ENDDO

    DO IS=-1, 1
       DO JS=-1, 1
          WEIGHTpb(IS, JS)=WEIGHTpb(IS, JS)/MY_SUM
       ENDDO
    ENDDO

    DO J=2, JMAR-1
       DO I=2, IMAR-1
          XF(I, J)=0.
          DO IS=-1, 1
             DO JS=-1, 1
                XF(I, J)=XF(I, J)+X(I+IS, J+JS)*WEIGHTpb(IS, JS)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO J=2, JMAR-1
       XF(1, J)=0.
       IS=IMAR-1
       DO JS=-1, 1 
          XF(1, J)=XF(1, J)+X(IS, J+JS)*WEIGHTpb(-1, JS)
       ENDDO
       DO IS=0, 1 
          DO JS=-1, 1 
             XF(1, J)=XF(1, J)+X(1+IS, J+JS)*WEIGHTpb(IS, JS)
          ENDDO
       ENDDO
       XF(IMAR, J)=XF(1, J)
    ENDDO

    DO I=1, IMAR
       XF(I, 1)=XF(I, 2)
       XF(I, JMAR)=XF(I, JMAR-1)
    ENDDO

    DO I=1, IMAR
       DO J=1, JMAR
          X(I, J)=XF(I, J)
       ENDDO
    ENDDO

  END SUBROUTINE MVA9

end module MVA9_m
