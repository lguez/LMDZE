module flxdtdq_m

  IMPLICIT none

contains

  SUBROUTINE flxdtdq(ktopm2, paph, ldcum, pten,  pmfus, pmfds, pmfuq, pmfdq, &
       pmful, pdmfup, pdmfdp,  pdpmel, dt_con, dq_con)

    USE dimphy, ONLY: klev, klon
    USE suphec_m, ONLY: rcpd, rg, rlmlt, rlstt, rlvtt, rtt

    !----------------------------------------------------------------------
    ! calculer les tendances T et Q
    !----------------------------------------------------------------------
    !     -----------------------------------------------------------------
    LOGICAL  llo1
    !
    REAL, intent(in):: pten(klon,klev), paph(klon,klev+1)
    REAL, intent(in):: pmfus(klon,klev), pmfuq(klon,klev), pmful(klon,klev)
    REAL pmfds(klon,klev), pmfdq(klon,klev)
    REAL pdmfup(klon,klev)
    REAL pdmfdp(klon,klev)
    REAL pdpmel(klon,klev)
    LOGICAL ldcum(klon)
    REAL dt_con(klon,klev), dq_con(klon,klev)
    !
    INTEGER ktopm2
    !
    INTEGER i, k
    REAL zalv, zdtdt, zdqdt
    !
    DO  k=ktopm2,klev-1
       DO  i = 1, klon
          IF (ldcum(i)) THEN
             llo1 = (pten(i,k)-RTT).GT.0.
             zalv = RLSTT
             IF (llo1) zalv = RLVTT
             zdtdt=RG/(paph(i,k+1)-paph(i,k))/RCPD &
                  *(pmfus(i,k+1)-pmfus(i,k) &
                  +pmfds(i,k+1)-pmfds(i,k) &
                  -RLMLT*pdpmel(i,k) &
                  -zalv*(pmful(i,k+1)-pmful(i,k)-pdmfup(i,k)-pdmfdp(i,k)) &
                  )
             dt_con(i,k)=zdtdt
             zdqdt=RG/(paph(i,k+1)-paph(i,k)) &
                  *(pmfuq(i,k+1)-pmfuq(i,k) &
                  +pmfdq(i,k+1)-pmfdq(i,k) &
                  +pmful(i,k+1)-pmful(i,k)-pdmfup(i,k)-pdmfdp(i,k))
             dq_con(i,k)=zdqdt
          ENDIF
       end DO
    end DO
    !
    k = klev
    DO  i = 1, klon
       IF (ldcum(i)) THEN
          llo1 = (pten(i,k)-RTT).GT.0.
          zalv = RLSTT
          IF (llo1) zalv = RLVTT
          zdtdt=-RG/(paph(i,k+1)-paph(i,k))/RCPD &
               *(pmfus(i,k)+pmfds(i,k)+RLMLT*pdpmel(i,k) &
               -zalv*(pmful(i,k)+pdmfup(i,k)+pdmfdp(i,k)))
          dt_con(i,k)=zdtdt
          zdqdt=-RG/(paph(i,k+1)-paph(i,k)) &
               *(pmfuq(i,k)+pmfdq(i,k)+pmful(i,k) &
               +pdmfup(i,k)+pdmfdp(i,k))
          dq_con(i,k)=zdqdt
       ENDIF
    end DO

  END SUBROUTINE flxdtdq

end module flxdtdq_m
