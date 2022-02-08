module grille_m_m

  IMPLICIT none

contains

  function grille_m(xdata, ydata, entree, x, y)

    ! From grid_atob.F, v 1.1.1.1, 2004/05/19 12:53:05
    
    ! Z. X. Li (1er avril 1994) (voir aussi A. Harzallah et L. Fairhead)

    ! M\'ethode na\"ive pour transformer un champ d'une grille fine \`a une
    ! grille grossi\`ere. Je consid\`ere que les nouveaux points occupent
    ! une zone adjacente qui comprend un ou plusieurs anciens points.

    ! Aucune pond\'eration n'est consid\'er\'ee. Cf. grille_m.txt.

    use jumble, only: assert_eq

    use  dist_sphe_m, only: dist_sphe

    ! Coordonn\'ees :
    REAL, intent(in):: xdata(:) ! (imdep)
    REAL, intent(in):: ydata(:) ! (jmdep)

    REAL, intent(in):: entree(:, :) ! (imdep, jmdep) champ \`a transformer

    ! Coordonn\'ees :
    REAL, intent(in):: x(:) ! (imar)
    REAL, intent(in):: y(:) ! (jmar)

    real grille_m(size(x), size(y)) ! (imar, jmar) champ transform\'e

    ! Local:
    INTEGER imdep, jmdep, imar, jmar
    INTEGER i, j, ii, jj
    REAL a(size(x)), b(size(x)) ! (imar)
    real c(size(y)), d(size(y)) ! (jmar)
    REAL number(size(x), size(y)) ! (imar, jmar) 
    REAL distans(size(xdata) * size(ydata)) ! (imdep * jmdep)
    INTEGER i_proche, j_proche, ij_proche
    REAL zzmin

    !-------------------------

    imdep = assert_eq(size(xdata), size(entree, 1), "grille_m")
    jmdep = assert_eq(size(ydata), size(entree, 2), "grille_m")
    imar = size(x)
    jmar = size(y)

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
          number(i, j) = 0.0
          grille_m(i, j) = 0.0
       ENDDO
    ENDDO

    ! Determiner la zone sur laquelle chaque ancien point se trouve

    DO ii = 1, imar
       DO jj = 1, jmar
          DO i = 1, imdep
             IF((xdata(i)-a(ii) >= 1.e-5.AND.xdata(i)-b(ii) <= 1.e-5).OR. &
                  (xdata(i)-a(ii) <= 1.e-5.AND.xdata(i)-b(ii) >= 1.e-5)) &
                  THEN
                DO j = 1, jmdep
                   IF((ydata(j)-c(jj) >= 1.e-5.AND.ydata(j)-d(jj) <= 1.e-5) &
                        .OR. (ydata(j)-c(jj) <= 1.e-5 .AND. &
                        ydata(j)-d(jj) >= 1.e-5)) THEN
                      number(ii, jj) = number(ii, jj) + 1.0
                      grille_m(ii, jj) = grille_m(ii, jj) + entree(i, j)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, imar
       DO j = 1, jmar
          IF (number(i, j) > 0.001) THEN
             grille_m(i, j) = grille_m(i, j) / number(i, j)
          ELSE
             ! Si aucun ancien point ne tombe sur une zone, c'est un probl\`eme
             CALL dist_sphe(x(i), y(j), xdata, ydata, imdep, jmdep, distans)
             ij_proche = 1
             zzmin = distans(ij_proche)
             DO ii = 2, imdep*jmdep
                IF (distans(ii) < zzmin) THEN
                   zzmin = distans(ii)
                   ij_proche = ii
                ENDIF
             ENDDO
             j_proche = (ij_proche-1)/imdep + 1
             i_proche = ij_proche - (j_proche-1)*imdep
             grille_m(i, j) = entree(i_proche, j_proche)
          ENDIF
       ENDDO
    ENDDO

    if (any(number <= 0.001)) print *, "problem in grille_m"

  END function grille_m

end module grille_m_m
