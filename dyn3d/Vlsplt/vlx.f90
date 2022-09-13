module vlx_m

  IMPLICIT NONE

contains

  SUBROUTINE vlx(q, pente_max, masse, u_m)

    ! Authors: P. Le Van, F. Hourdin, F. Forget

    ! Sch\'ema d'advection "pseudo-amont".

    use dimensions, only: iim, llm
    use paramet_m, only: ip1jmp1, iip1, iip2, ip1jm

    REAL, intent(inout):: q(ip1jmp1, llm)
    REAL, intent(in):: pente_max
    real, intent(inout):: masse(ip1jmp1, llm)
    REAL, intent(in):: u_m(ip1jmp1, llm)

    ! Local:
    INTEGER ij, l, j, i, iju, ijq, indu(ip1jmp1), niju
    INTEGER n0, iadvplus(ip1jmp1, llm), nl(llm)
    REAL new_m, zu_m, zdum(ip1jmp1, llm)
    REAL dxq(ip1jmp1, llm), dxqu(iip2:ip1jm)
    REAL zz(ip1jmp1)
    REAL adxqu(iip2:ip1jm), dxqmax(ip1jmp1, llm)
    REAL u_mq(ip1jmp1, llm)

    !-----------------------------------------------------

    ! calcul de la pente a droite et a gauche de la maille

    IF (pente_max > - 1e-5) THEN
       ! calcul des pentes avec limitation, Van Leer scheme I:

       ! calcul de la pente aux points u
       DO l = 1, llm
          DO ij = iip2, ip1jm - 1
             dxqu(ij) = q(ij + 1, l) - q(ij, l)
          ENDDO
          DO ij = iip1 + iip1, ip1jm, iip1
             dxqu(ij) = dxqu(ij - iim)
          ENDDO

          DO ij = iip2, ip1jm
             adxqu(ij) = abs(dxqu(ij))
          ENDDO

          ! calcul de la pente maximum dans la maille en valeur absolue

          DO ij = iip2 + 1, ip1jm
             dxqmax(ij, l) = pente_max * min(adxqu(ij - 1), adxqu(ij))
          ENDDO

          DO ij = iip1 + iip1, ip1jm, iip1
             dxqmax(ij - iim, l) = dxqmax(ij, l)
          ENDDO

          DO ij = iip2 + 1, ip1jm
             IF (dxqu(ij - 1) * dxqu(ij) > 0.) THEN
                dxq(ij, l) = dxqu(ij - 1) + dxqu(ij)
             ELSE
                ! extremum local
                dxq(ij, l) = 0.
             ENDIF
             dxq(ij, l) = 0.5 * dxq(ij, l)
             dxq(ij, l) =  sign(min(abs(dxq(ij, l)), dxqmax(ij, l)), dxq(ij, l))
          ENDDO
       ENDDO
    ELSE
       ! Pentes produits:

       DO l = 1, llm
          DO ij = iip2, ip1jm - 1
             dxqu(ij) = q(ij + 1, l) - q(ij, l)
          ENDDO
          DO ij = iip1 + iip1, ip1jm, iip1
             dxqu(ij) = dxqu(ij - iim)
          ENDDO

          DO ij = iip2 + 1, ip1jm
             zz(ij) = dxqu(ij - 1) * dxqu(ij)
             zz(ij) = zz(ij) + zz(ij)
             IF (zz(ij) > 0) THEN
                dxq(ij, l) = zz(ij) / (dxqu(ij - 1) + dxqu(ij))
             ELSE
                ! extremum local
                dxq(ij, l) = 0.
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! bouclage de la pente en iip1:

    DO l = 1, llm
       DO ij = iip1 + iip1, ip1jm, iip1
          dxq(ij - iim, l) = dxq(ij, l)
       ENDDO
       DO ij = 1, ip1jmp1
          iadvplus(ij, l) = 0
       ENDDO
    ENDDO

    ! calcul des flux a gauche et a droite

    ! on cumule le flux correspondant a toutes les mailles dont la masse
    ! au travers de la paroi pENDant le pas de temps.

    DO l = 1, llm
       DO ij = iip2, ip1jm - 1
          IF (u_m(ij, l) > 0.) THEN
             zdum(ij, l) = 1. - u_m(ij, l) / masse(ij, l)
             u_mq(ij, l) = u_m(ij, l) &
                  * (q(ij, l) + 0.5 * zdum(ij, l) * dxq(ij, l))
          ELSE
             zdum(ij, l) = 1. + u_m(ij, l) / masse(ij + 1, l)
             u_mq(ij, l) = u_m(ij, l) &
                  * (q(ij + 1, l) - 0.5 * zdum(ij, l) * dxq(ij + 1, l))
          ENDIF
       ENDDO
    ENDDO

    ! detection des points ou on advecte plus que la masse de la
    ! maille
    DO l = 1, llm
       DO ij = iip2, ip1jm - 1
          IF (zdum(ij, l) < 0) THEN
             iadvplus(ij, l) = 1
             u_mq(ij, l) = 0.
          ENDIF
       ENDDO
    ENDDO

    DO l = 1, llm
       DO ij = iip1 + iip1, ip1jm, iip1
          iadvplus(ij, l) = iadvplus(ij - iim, l)
       ENDDO
    ENDDO

    ! traitement special pour le cas ou on advecte en longitude plus
    ! que le contenu de la maille.  cette partie est mal vectorisee.

    ! calcul du nombre de maille sur lequel on advecte plus que la maille.

    n0 = 0
    DO l = 1, llm
       nl(l) = 0
       DO ij = iip2, ip1jm
          nl(l) = nl(l) + iadvplus(ij, l)
       ENDDO
       n0 = n0 + nl(l)
    ENDDO

    IF (n0 > 0) THEN
       DO l = 1, llm
          IF (nl(l) > 0) THEN
             iju = 0
             ! indicage des mailles concernees par le traitement special
             DO ij = iip2, ip1jm
                IF (iadvplus(ij, l) == 1.and.mod(ij, iip1) /= 0) THEN
                   iju = iju + 1
                   indu(iju) = ij
                ENDIF
             ENDDO
             niju = iju

             ! traitement des mailles
             DO iju = 1, niju
                ij = indu(iju)
                j = (ij - 1) / iip1 + 1
                zu_m = u_m(ij, l)
                u_mq(ij, l) = 0.
                IF (zu_m > 0.) THEN
                   ijq = ij
                   i = ijq - (j - 1) * iip1
                   ! accumulation pour les mailles completements advectees
                   do while(zu_m > masse(ijq, l))
                      u_mq(ij, l) = u_mq(ij, l) + q(ijq, l) * masse(ijq, l)
                      zu_m = zu_m - masse(ijq, l)
                      i = mod(i - 2 + iim, iim) + 1
                      ijq = (j - 1) * iip1 + i
                   ENDDO
                   ! ajout de la maille non completement advectee
                   u_mq(ij, l) = u_mq(ij, l) + zu_m * (q(ijq, l) + 0.5 * (1. &
                        - zu_m / masse(ijq, l)) * dxq(ijq, l))
                ELSE
                   ijq = ij + 1
                   i = ijq - (j - 1) * iip1
                   ! accumulation pour les mailles completements advectees
                   do while(- zu_m > masse(ijq, l))
                      u_mq(ij, l) = u_mq(ij, l) - q(ijq, l) * masse(ijq, l)
                      zu_m = zu_m + masse(ijq, l)
                      i = mod(i, iim) + 1
                      ijq = (j - 1) * iip1 + i
                   ENDDO
                   ! ajout de la maille non completement advectee
                   u_mq(ij, l) = u_mq(ij, l) + zu_m * (q(ijq, l) - 0.5 * (1. &
                        + zu_m / masse(ijq, l)) * dxq(ijq, l))
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDIF

    ! bouclage en latitude
    DO l = 1, llm
       DO ij = iip1 + iip1, ip1jm, iip1
          u_mq(ij, l) = u_mq(ij - iim, l)
       ENDDO
    ENDDO

    ! calcul des tENDances

    DO l = 1, llm
       DO ij = iip2 + 1, ip1jm
          new_m = masse(ij, l) + u_m(ij - 1, l) - u_m(ij, l)
          q(ij, l) = (q(ij, l) * masse(ij, l) + u_mq(ij - 1, l) &
               - u_mq(ij, l)) / new_m
          masse(ij, l) = new_m
       ENDDO

       DO ij = iip1 + iip1, ip1jm, iip1
          q(ij - iim, l) = q(ij, l)
          masse(ij - iim, l) = masse(ij, l)
       ENDDO
    ENDDO

  END SUBROUTINE vlx

end module vlx_m
