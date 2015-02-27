module heavyside_m

  IMPLICIT NONE

contains

  real FUNCTION heavyside(a)

    ! From LMDZ4/libf/dyn3d/heavyside.F,v 1.1.1.1 2004/05/19 12:53:06

    ! ...   P. Le Van  ....

    real, intent(in):: a

    !-------------------------------------------------

    IF (a<=0.) THEN
       heavyside = 0.
    ELSE
       heavyside = 1.
    END IF

  END FUNCTION heavyside

end module heavyside_m
