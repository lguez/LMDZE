module minmaxqfi_m

    IMPLICIT none

contains

  SUBROUTINE minmaxqfi(zq, qmin, qmax, comment)

    ! From phylmd/minmaxqfi.F, version 1.1.1.1 2004/05/19 12:53:09

    use jumble, only: assert

    use dimensions, only: llm
    use dimphy, only: klon

    real, intent(in):: zq(:, :), qmin, qmax
    CHARACTER(len=*), intent(in):: comment

    ! Variables local to the procedure:

    INTEGER jadrs(klon), jbad, k, i

    !---------------------------------

    call assert(shape(zq) == (/klon, llm/), "minmaxqfi")

    DO k = 1, llm
       jbad = 0
       DO i = 1, klon
          IF (zq(i, k) > qmax .OR. zq(i, k) < qmin) THEN
             jbad = jbad + 1
             jadrs(jbad) = i
          ENDIF
       ENDDO
       IF (jbad > 0) THEN
          PRINT *, comment
          DO i = 1, jbad
             PRINT *, "zq(", jadrs(i), ", ", k, ") = ", zq(jadrs(i), k)
          ENDDO
       ENDIF
    ENDDO

  end SUBROUTINE minmaxqfi

end module minmaxqfi_m
