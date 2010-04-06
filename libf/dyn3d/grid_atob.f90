module grid_atob

  ! From grid_atob.F,v 1.1.1.1 2004/05/19 12:53:05

  IMPLICIT none

contains

  function grille_m(xdata, ydata, entree, x, y)

    !=======================================================================
    ! Z. X. Li (le 1 avril 1994) (voir aussi A. Harzallah et L. Fairhead)

    ! Méthode naïve pour transformer un champ d'une grille fine à une
    ! grille grossière. Je considère que les nouveaux points occupent
    ! une zone adjacente qui comprend un ou plusieurs anciens points

    ! Aucune pondération n'est considérée (voir grille_p)

    !           (c)
    !        ----d-----
    !        | . . . .|
    !        |        |
    !     (b)a . * . .b(a)
    !        |        |
    !        | . . . .|
    !        ----c-----
    !           (d)
    !=======================================================================
    ! INPUT:
    !        imdep, jmdep: dimensions X et Y pour depart
    !        xdata, ydata: coordonnees X et Y pour depart
    !        entree: champ d'entree a transformer
    ! OUTPUT:
    !        imar, jmar: dimensions X et Y d'arrivee
    !        x, y: coordonnees X et Y d'arrivee
    !        grille_m: champ de sortie deja transforme
    !=======================================================================

    use numer_rec, only: assert_eq

    REAL, intent(in):: xdata(:),ydata(:) 
    REAL, intent(in):: entree(:, :)
    REAL, intent(in):: x(:), y(:)

    real grille_m(size(x), size(y))

    ! Variables local to the procedure:
    INTEGER imdep, jmdep, imar, jmar
    INTEGER i, j, ii, jj
    REAL a(2200),b(2200),c(1100),d(1100)
    REAL number(2200,1100)
    REAL distans(2200*1100)
    INTEGER i_proche, j_proche, ij_proche
    REAL zzmin

    !-------------------------

    print *, "Call sequence information: grille_m"

    imdep = assert_eq(size(xdata), size(entree, 1), "grille_m")
    jmdep = assert_eq(size(ydata), size(entree, 2), "grille_m")
    imar = size(x)
    jmar = size(y)

    IF (imar.GT.2200 .OR. jmar.GT.1100) THEN
       PRINT*, 'imar ou jmar trop grand', imar, jmar
       STOP 1
    ENDIF

    ! Calculer les limites des zones des nouveaux points

    a(1) = x(1) - (x(2)-x(1))/2.0
    b(1) = (x(1)+x(2))/2.0
    DO i = 2, imar-1
       a(i) = b(i-1)
       b(i) = (x(i)+x(i+1))/2.0
    ENDDO
    a(imar) = b(imar-1)
    b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

    c(1) = y(1) - (y(2)-y(1))/2.0
    d(1) = (y(1)+y(2))/2.0
    DO j = 2, jmar-1
       c(j) = d(j-1)
       d(j) = (y(j)+y(j+1))/2.0
    ENDDO
    c(jmar) = d(jmar-1)
    d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

    DO i = 1, imar
       DO j = 1, jmar
          number(i,j) = 0.0
          grille_m(i,j) = 0.0
       ENDDO
    ENDDO

    ! Determiner la zone sur laquelle chaque ancien point se trouve

    DO ii = 1, imar
       DO jj = 1, jmar
          DO i = 1, imdep
             IF( ( xdata(i)-a(ii) >= 1.e-5.AND.xdata(i)-b(ii) <= 1.e-5 ).OR. &
                  (xdata(i)-a(ii) <= 1.e-5.AND.xdata(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmdep
                   IF((ydata(j)-c(jj) >= 1.e-5.AND.ydata(j)-d(jj) <= 1.e-5 ) &
                        .OR. (ydata(j)-c(jj) <= 1.e-5 .AND. &
                        ydata(j)-d(jj) >= 1.e-5)) THEN
                      number(ii,jj) = number(ii,jj) + 1.0
                      grille_m(ii,jj) = grille_m(ii,jj) + entree(i,j)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Si aucun ancien point ne tombe sur une zone, c'est un probleme

    DO i = 1, imar
       DO j = 1, jmar
          IF (number(i,j) .GT. 0.001) THEN
             grille_m(i,j) = grille_m(i,j) / number(i,j)
          ELSE
             PRINT*, 'probleme,i,j=', i,j
             CALL dist_sphe(x(i),y(j),xdata,ydata,imdep,jmdep,distans)
             ij_proche = 1
             zzmin = distans(ij_proche)
             DO ii = 2, imdep*jmdep
                IF (distans(ii).LT.zzmin) THEN
                   zzmin = distans(ii)
                   ij_proche = ii
                ENDIF
             ENDDO
             j_proche = (ij_proche-1)/imdep + 1
             i_proche = ij_proche - (j_proche-1)*imdep
             PRINT*, "solution:", ij_proche, i_proche, j_proche
             grille_m(i,j) = entree(i_proche,j_proche)
          ENDIF
       ENDDO
    ENDDO

  END function grille_m

  !**************************************************

  SUBROUTINE dist_sphe(rf_lon,rf_lat,rlon,rlat,im,jm,distance)

    ! Auteur: Laurent Li (le 30 decembre 1996)

    ! Ce programme calcule la distance minimale (selon le grand cercle)
    ! entre deux points sur la terre

    INTEGER, intent(in):: im, jm ! dimensions
    REAL, intent(in):: rf_lon ! longitude du point de reference (degres)
    REAL, intent(in):: rf_lat ! latitude du point de reference (degres)
    REAL, intent(in):: rlon(im), rlat(jm) ! longitude et latitude des points

    REAL, intent(out):: distance(im,jm) ! distances en metre

    REAL rlon1, rlat1
    REAL rlon2, rlat2
    REAL dist
    REAL pa, pb, p, pi

    REAL radius
    PARAMETER (radius=6371229.)
    integer i, j

    pi = 4.0 * ATAN(1.0)

    DO j = 1, jm
       DO  i = 1, im
      
          rlon1=rf_lon
          rlat1=rf_lat
          rlon2=rlon(i)
          rlat2=rlat(j)
          pa = pi/2.0 - rlat1*pi/180.0 ! dist. entre pole n et point a
          pb = pi/2.0 - rlat2*pi/180.0 ! dist. entre pole n et point b
          p = (rlon1-rlon2)*pi/180.0 ! angle entre a et b (leurs meridiens)
      
          dist = ACOS( COS(pa)*COS(pb) + SIN(pa)*SIN(pb)*COS(p))
          dist = radius * dist
          distance(i,j) = dist
      
       end DO
    end DO

  END SUBROUTINE dist_sphe

end module grid_atob
