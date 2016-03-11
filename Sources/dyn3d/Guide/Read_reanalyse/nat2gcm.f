module nat2gcm_m

  implicit none

contains

  subroutine nat2gcm(pk, u, v, t)

    ! Passage aux variables du mod\`ele (vents covariants,
    ! temp\'erature potentielle et humidit\'e sp\'ecifique).

    use comconst, only: cpp
    use comgeom, only: cu_2d, cv_2d
    use dimens_m, only: iim, jjm, llm
    use paramet_m, only: iip1, jjp1

    real, intent(in):: pk(iip1, jjp1, llm)
    real, intent(inout):: u(iip1, jjp1, llm), v(iip1, jjm, llm)
    real, intent(inout):: t(iip1, jjp1, llm)

    ! Local:
    integer i, j, l

    !----------------------------------------------------------------------

    print *, "Call sequence information: nat2gcm"

    ! calcul de ucov et de la temperature potentielle
    do l = 1, llm
       do j = 1, jjp1
          do i = 1, iim
             u(i, j, l) = u(i, j, l) * cu_2d(i, j)
             t(i, j, l) = t(i, j, l) * cpp / pk(i, j, l)
          enddo
          u(iip1, j, l) = u(1, j, l)
          t(iip1, j, l) = t(1, j, l)
       enddo
       do i = 1, iip1
          u(i, 1, l) = 0.
          u(i, jjp1, l) = 0.
          t(i, 1, l) = t(1, 1, l)
          t(i, jjp1, l) = t(1, jjp1, l)
       enddo
    enddo

    do l = 1, llm
       do j = 1, jjm
          do i = 1, iim
             v(i, j, l) = v(i, j, l) * cv_2d(i, j)
          enddo
          v(iip1, j, l) = v(1, j, l)
       enddo
    enddo

  end subroutine nat2gcm

end module nat2gcm_m
