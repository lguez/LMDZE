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

  SUBROUTINE grille_p(imdep, jmdep, xdata, ydata, entree, &
       imar, jmar, x, y, sortie)
    !=======================================================================
    ! z.x.li (le 1 avril 1994) (voir aussi A. Harzallah et L. Fairhead)

    ! Methode naive pour transformer un champ d'une grille fine a une
    ! grille grossiere. Je considere que les nouveaux points occupent
    ! une zone adjacente qui comprend un ou plusieurs anciens points

    ! Consideration de la distance des points (voir grille_m)

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
    !        sortie: champ de sortie deja transforme
    !=======================================================================

    INTEGER imdep, jmdep
    REAL xdata(imdep),ydata(jmdep) 
    REAL entree(imdep,jmdep)

    INTEGER imar, jmar
    REAL x(imar),y(jmar)
    REAL sortie(imar,jmar)

    INTEGER i, j, ii, jj
    REAL a(400),b(400),c(200),d(200)
    REAL number(400,200)
    INTEGER indx(400,200), indy(400,200)
    REAL dist(400,200), distsom(400,200)

    IF (imar.GT.400 .OR. jmar.GT.200) THEN
       PRINT*, 'imar ou jmar trop grand', imar, jmar
       STOP 1
    ENDIF

    IF (imdep.GT.400 .OR. jmdep.GT.200) THEN
       PRINT*, 'imdep ou jmdep trop grand', imdep, jmdep
       STOP 1
    ENDIF

    ! calculer les bords a et b de la nouvelle grille

    a(1) = x(1) - (x(2)-x(1))/2.0
    b(1) = (x(1)+x(2))/2.0
    DO i = 2, imar-1
       a(i) = b(i-1)
       b(i) = (x(i)+x(i+1))/2.0
    ENDDO
    a(imar) = b(imar-1)
    b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

    ! calculer les bords c et d de la nouvelle grille

    c(1) = y(1) - (y(2)-y(1))/2.0
    d(1) = (y(1)+y(2))/2.0
    DO j = 2, jmar-1
       c(j) = d(j-1)
       d(j) = (y(j)+y(j+1))/2.0
    ENDDO
    c(jmar) = d(jmar-1)
    d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0

    ! trouver les indices (indx,indy) de la nouvelle grille sur laquelle
    ! un point de l'ancienne grille est tombe.

    !  .....  Modif  P. Le Van ( 23/08/95 )  ....

    DO ii = 1, imar
       DO jj = 1, jmar
          DO i = 1, imdep
             IF( ( xdata(i)-a(ii) >= 1.e-5.AND.xdata(i)-b(ii) <= 1.e-5 ).OR. &
                  (   xdata(i)-a(ii) <= 1.e-5.AND.xdata(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmdep
                   IF( (ydata(j)-c(jj) >= 1.e-5.AND.ydata(j)-d(jj) <= 1.e-5 ).OR. &
                        (  ydata(j)-c(jj) <= 1.e-5.AND.ydata(j)-d(jj) >= 1.e-5 )   ) &
                        THEN
                      indx(i,j) = ii
                      indy(i,j) = jj
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! faire une verification

    DO i = 1, imdep
       DO j = 1, jmdep
          IF (indx(i,j).GT.imar .OR. indy(i,j).GT.jmar) THEN
             PRINT*, 'Probleme grave,i,j,indx,indy=', &
                  i,j,indx(i,j),indy(i,j)
             stop 1
          ENDIF
       ENDDO
    ENDDO

    ! calculer la distance des anciens points avec le nouveau point,
    ! on prend ensuite une sorte d'inverse pour ponderation.

    DO i = 1, imar
       DO j = 1, jmar
          number(i,j) = 0.0
          distsom(i,j) = 0.0
       ENDDO
    ENDDO
    DO i = 1, imdep
       DO j = 1, jmdep
          dist(i,j) = SQRT ( (xdata(i)-x(indx(i,j)))**2 &
               +(ydata(j)-y(indy(i,j)))**2 )
          distsom(indx(i,j),indy(i,j)) = distsom(indx(i,j),indy(i,j)) &
               + dist(i,j)
          number(indx(i,j),indy(i,j)) = number(indx(i,j),indy(i,j)) +1.
       ENDDO
    ENDDO
    DO i = 1, imdep
       DO j = 1, jmdep
          dist(i,j) = 1.0 - dist(i,j)/distsom(indx(i,j),indy(i,j))
       ENDDO
    ENDDO

    DO i = 1, imar
       DO j = 1, jmar
          number(i,j) = 0.0
          sortie(i,j) = 0.0
       ENDDO
    ENDDO
    DO i = 1, imdep
       DO j = 1, jmdep
          sortie(indx(i,j),indy(i,j)) = sortie(indx(i,j),indy(i,j)) &
               + entree(i,j) * dist(i,j)
          number(indx(i,j),indy(i,j)) = number(indx(i,j),indy(i,j)) &
               + dist(i,j)
       ENDDO
    ENDDO
    DO i = 1, imar
       DO j = 1, jmar
          IF (number(i,j) .GT. 0.001) THEN
             sortie(i,j) = sortie(i,j) / number(i,j)
          ELSE
             PRINT*, 'probleme,i,j=', i,j
             STOP 1
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE grille_p

  !******************************************************************

  SUBROUTINE mask_c_o(imdep, jmdep, xdata, ydata, relief, &
       imar, jmar, x, y, mask)
    !=======================================================================
    ! z.x.li (le 1 avril 1994): A partir du champ de relief, on fabrique
    !                           un champ indicateur (masque) terre/ocean
    !                           terre:1; ocean:0

    ! Methode naive (voir grille_m)
    !=======================================================================

    INTEGER imdep, jmdep
    REAL xdata(imdep),ydata(jmdep) 
    REAL relief(imdep,jmdep)

    INTEGER imar, jmar
    REAL x(imar),y(jmar)
    REAL mask(imar,jmar)

    INTEGER i, j, ii, jj
    REAL a(2200),b(2200),c(1100),d(1100)
    REAL num_tot(2200,1100), num_oce(2200,1100)

    IF (imar.GT.2200 .OR. jmar.GT.1100) THEN
       PRINT*, 'imar ou jmar trop grand', imar, jmar
       STOP 1
    ENDIF

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
          num_oce(i,j) = 0.0
          num_tot(i,j) = 0.0
       ENDDO
    ENDDO

    !  .....  Modif  P. Le Van ( 23/08/95 )  ....

    DO ii = 1, imar
       DO jj = 1, jmar
          DO i = 1, imdep
             IF( ( xdata(i)-a(ii) >= 1.e-5.AND.xdata(i)-b(ii) <= 1.e-5 ).OR. &
                  (   xdata(i)-a(ii) <= 1.e-5.AND.xdata(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmdep
                   IF( (ydata(j)-c(jj) >= 1.e-5.AND.ydata(j)-d(jj) <= 1.e-5 ).OR. &
                        (  ydata(j)-c(jj) <= 1.e-5.AND.ydata(j)-d(jj) >= 1.e-5 )   ) &
                        THEN
                      num_tot(ii,jj) = num_tot(ii,jj) + 1.0
                      IF (.NOT. ( relief(i,j) - 0.9>= 1.e-5 ) ) &
                           num_oce(ii,jj) = num_oce(ii,jj) + 1.0
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, imar
       DO j = 1, jmar
          IF (num_tot(i,j) .GT. 0.001) THEN
             IF ( num_oce(i,j)/num_tot(i,j) - 0.5  >=  1.e-5 ) THEN
                mask(i,j) = 0.
             ELSE
                mask(i,j) = 1.
             ENDIF
          ELSE
             PRINT*, 'probleme,i,j=', i,j
             STOP 1
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mask_c_o

  ! *************************************

  real function rugosite(xdata, ydata, entree, x, y, mask)

    ! Z. X. Li (le 1 avril 1994): Transformer la longueur de rugosite d'une
    ! grille fine a une grille grossiere. Sur l'ocean, on impose une valeur
    ! fixe (0.001m).

    ! Methode naive (voir grille_m)

    use numer_rec, only: assert_eq

    REAL, intent(in):: xdata(:), ydata(:), entree(:,:), x(:), y(:), mask(:,:)

    dimension rugosite(size(mask, 1), size(mask, 2))
    
    ! Variables local to the procedure:
    INTEGER imdep, jmdep
    INTEGER imar, jmar
    INTEGER i, j, ii, jj
    REAL a(400),b(400),c(400),d(400)
    REAL num_tot(400,400)
    REAL distans(400*400)
    INTEGER i_proche, j_proche, ij_proche
    REAL zzmin

    ! --------------------

    imdep = assert_eq(size(xdata), size(entree, 1), "rugosite")
    jmdep = assert_eq(size(ydata), size(entree, 2), "rugosite")
    imar = assert_eq(size(x), size(mask, 1), "rugosite")
    jmar = assert_eq(size(y), size(mask, 2), "rugosite")

    IF (imar.GT.400 .OR. jmar.GT.400) THEN
       PRINT*, 'imar ou jmar trop grand', imar, jmar
       STOP 1
    ENDIF

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
          num_tot(i,j) = 0.0
          rugosite(i,j) = 0.0
       ENDDO
    ENDDO

    !  .....  Modif  P. Le Van ( 23/08/95 )  ....

    DO ii = 1, imar
       DO jj = 1, jmar
          DO i = 1, imdep
             IF( ( xdata(i)-a(ii) >= 1.e-5.AND.xdata(i)-b(ii) <= 1.e-5 ).OR. &
                  (   xdata(i)-a(ii) <= 1.e-5.AND.xdata(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmdep
                   IF( (ydata(j)-c(jj) >= 1.e-5.AND.ydata(j)-d(jj) <= 1.e-5 ).OR. &
                        (  ydata(j)-c(jj) <= 1.e-5.AND.ydata(j)-d(jj) >= 1.e-5 )   ) &
                        THEN
                      rugosite(ii,jj)  = rugosite(ii,jj) + LOG(entree(i,j))
                      num_tot(ii,jj) = num_tot(ii,jj) + 1.0
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, imar
       DO j = 1, jmar
          IF (NINT(mask(i,j)).EQ.1) THEN
             IF (num_tot(i,j) .GT. 0.0) THEN
                rugosite(i,j) = rugosite(i,j) / num_tot(i,j)
                rugosite(i,j) = EXP(rugosite(i,j))
             ELSE
                PRINT*, 'probleme,i,j=', i,j
                !cc            STOP 1
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
                rugosite(i,j) = entree(i_proche,j_proche)
             ENDIF
          ELSE
             rugosite(i,j) = 0.001
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END function rugosite

  !************************************

  real function sea_ice(xdata, ydata, glace01, x, y)

    !=======================================================================
    ! z.x.li (le 1 avril 1994): Transformer un champ d'indicateur de la
    ! glace (1, sinon 0) d'une grille fine a un champ de fraction de glace
    ! (entre 0 et 1) dans une grille plus grossiere.

    ! Methode naive (voir grille_m)
    !=======================================================================

    use numer_rec, only: assert_eq

    REAL, intent(in):: xdata(:),ydata(:) 
    REAL, intent(in):: glace01(:,:)
    REAL, intent(in):: x(:),y(:)
    dimension sea_ice(size(x), size(y))

    ! Variables local to the procedure:
    INTEGER imdep, jmdep
    INTEGER imar, jmar
    INTEGER i, j, ii, jj
    REAL a(400),b(400),c(400),d(400)
    REAL num_tot(400,400), num_ice(400,400)
    REAL distans(400*400)
    INTEGER i_proche, j_proche, ij_proche
    REAL zzmin

    !------------------------------

    imdep = assert_eq(size(xdata), size(glace01, 1), "sea_ice")
    jmdep = assert_eq(size(ydata), size(glace01, 2), "sea_ice")
    imar = size(x)
    jmar = size(y)

    IF (imar.GT.400 .OR. jmar.GT.400) THEN
       PRINT*, 'imar ou jmar trop grand', imar, jmar
       STOP 1
    ENDIF

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
          num_ice(i,j) = 0.0
          num_tot(i,j) = 0.0
       ENDDO
    ENDDO

    !  .....  Modif  P. Le Van ( 23/08/95 )  ....

    DO ii = 1, imar
       DO jj = 1, jmar
          DO i = 1, imdep
             IF( ( xdata(i)-a(ii) >= 1.e-5.AND.xdata(i)-b(ii) <= 1.e-5 ).OR. &
                  (   xdata(i)-a(ii) <= 1.e-5.AND.xdata(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmdep
                   IF( (ydata(j)-c(jj) >= 1.e-5.AND.ydata(j)-d(jj) <= 1.e-5 ).OR. &
                        (  ydata(j)-c(jj) <= 1.e-5.AND.ydata(j)-d(jj) >= 1.e-5 )   ) &
                        THEN
                      num_tot(ii,jj) = num_tot(ii,jj) + 1.0
                      IF (NINT(glace01(i,j)).EQ.1 )  &
                           num_ice(ii,jj) = num_ice(ii,jj) + 1.0
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, imar
       DO j = 1, jmar
          IF (num_tot(i,j) .GT. 0.001) THEN
             IF (num_ice(i,j).GT.0.001) THEN
                sea_ice(i,j) = num_ice(i,j) / num_tot(i,j)
             ELSE
                sea_ice(i,j) = 0.0
             ENDIF
          ELSE
             PRINT*, 'probleme,i,j=', i,j
             !cc           STOP 1
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
             IF (NINT(glace01(i_proche,j_proche)).EQ.1 ) THEN
                sea_ice(i,j) = 1.0
             ELSE
                sea_ice(i,j) = 0.0
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END function sea_ice

  !*************************************

  SUBROUTINE rugsoro(imrel, jmrel, xrel, yrel, relief, immod, jmmod, xmod, &
       ymod, rugs)

    ! Calcule la longueur de rugosite liee au relief en utilisant
    ! l'ecart-type dans une maille de 1x1.

    INTEGER, intent(in):: imrel, jmrel
    REAL, intent(in):: xrel(imrel),yrel(jmrel)
    REAL, intent(in):: relief(imrel,jmrel)

    INTEGER, intent(in):: immod, jmmod
    REAL, intent(in):: xmod(immod),ymod(jmmod)
    REAL, intent(out):: rugs(immod,jmmod)

    REAL zzmin
    REAL amin, AMAX
    INTEGER imtmp, jmtmp
    PARAMETER (imtmp=360,jmtmp=180)
    REAL xtmp(imtmp), ytmp(jmtmp)
    double precision cham1tmp(imtmp,jmtmp), cham2tmp(imtmp,jmtmp)
    REAL zzzz

    INTEGER i, j, ii, jj
    REAL a(2200),b(2200),c(1100),d(1100)
    REAL number(2200,1100)

    REAL distans(400*400)
    INTEGER i_proche, j_proche, ij_proche

    !---------------------------------------------------------

    IF (immod.GT.2200 .OR. jmmod.GT.1100) THEN
       PRINT*, 'immod ou jmmod trop grand', immod, jmmod
       STOP 1
    ENDIF

    ! Calculs intermediares:

    xtmp(1) = -180.0 + 360.0/FLOAT(imtmp) / 2.0
    DO i = 2, imtmp
       xtmp(i) = xtmp(i-1) + 360.0/FLOAT(imtmp)
    ENDDO
    DO i = 1, imtmp
       xtmp(i) = xtmp(i) /180.0 * 4.0*ATAN(1.0)
    ENDDO
    ytmp(1) = -90.0 + 180.0/FLOAT(jmtmp) / 2.0
    DO j = 2, jmtmp
       ytmp(j) = ytmp(j-1) + 180.0/FLOAT(jmtmp)
    ENDDO
    DO j = 1, jmtmp
       ytmp(j) = ytmp(j) /180.0 * 4.0*ATAN(1.0)
    ENDDO

    a(1) = xtmp(1) - (xtmp(2)-xtmp(1))/2.0
    b(1) = (xtmp(1)+xtmp(2))/2.0
    DO i = 2, imtmp-1
       a(i) = b(i-1)
       b(i) = (xtmp(i)+xtmp(i+1))/2.0
    ENDDO
    a(imtmp) = b(imtmp-1)
    b(imtmp) = xtmp(imtmp) + (xtmp(imtmp)-xtmp(imtmp-1))/2.0

    c(1) = ytmp(1) - (ytmp(2)-ytmp(1))/2.0
    d(1) = (ytmp(1)+ytmp(2))/2.0
    DO j = 2, jmtmp-1
       c(j) = d(j-1)
       d(j) = (ytmp(j)+ytmp(j+1))/2.0
    ENDDO
    c(jmtmp) = d(jmtmp-1)
    d(jmtmp) = ytmp(jmtmp) + (ytmp(jmtmp)-ytmp(jmtmp-1))/2.0

    DO i = 1, imtmp
       DO j = 1, jmtmp
          number(i,j) = 0.0
          cham1tmp(i,j) = 0d0
          cham2tmp(i,j) = 0d0
       ENDDO
    ENDDO

    !  .....  Modif  P. Le Van ( 23/08/95 )  ....

    DO ii = 1, imtmp
       DO jj = 1, jmtmp
          DO i = 1, imrel
             IF( ( xrel(i)-a(ii) >= 1.e-5.AND.xrel(i)-b(ii) <= 1.e-5 ).OR. &
                  (   xrel(i)-a(ii) <= 1.e-5.AND.xrel(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmrel
                   IF ((yrel(j)-c(jj) >= 1.e-5.AND.yrel(j)-d(jj) <= 1.e-5 ) &
                        .OR. (yrel(j)-c(jj) <= 1.e-5 .AND. &
                        yrel(j)-d(jj) >= 1.e-5 )   ) &
                        THEN
                      number(ii,jj) = number(ii,jj) + 1.0
                      cham1tmp(ii,jj) = cham1tmp(ii,jj) + relief(i,j)
                      cham2tmp(ii,jj) = cham2tmp(ii,jj)  &
                           + relief(i,j) * relief(i,j)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, imtmp
       DO j = 1, jmtmp
          IF (number(i,j) .GT. 0.001) THEN
             cham1tmp(i,j) = cham1tmp(i,j) / number(i,j)
             cham2tmp(i,j) = cham2tmp(i,j) / number(i,j)
             zzzz = cham2tmp(i,j) - cham1tmp(i,j)**2
             if (zzzz .lt. 0.0) then
                if (zzzz .gt. -7.5) then
                   zzzz = 0.0
                   print*,'Pb rugsoro, -7.5 < zzzz < 0, => zzz = 0.0'
                else
                   stop 'Pb rugsoro, zzzz <-7.5'
                endif
             endif
             cham2tmp(i,j) = SQRT(zzzz)
          ELSE
             PRINT*, 'probleme,i,j=', i,j
             STOP 1
          ENDIF
       ENDDO
    ENDDO

    amin = cham2tmp(1,1)
    AMAX = cham2tmp(1,1)
    DO j = 1, jmtmp
       DO i = 1, imtmp
          IF (cham2tmp(i,j).GT.AMAX) AMAX = cham2tmp(i,j)
          IF (cham2tmp(i,j).LT.amin) amin = cham2tmp(i,j)
       ENDDO
    ENDDO
    PRINT*, 'Ecart-type 1x1:', amin, AMAX

    a(1) = xmod(1) - (xmod(2)-xmod(1))/2.0
    b(1) = (xmod(1)+xmod(2))/2.0
    DO i = 2, immod-1
       a(i) = b(i-1)
       b(i) = (xmod(i)+xmod(i+1))/2.0
    ENDDO
    a(immod) = b(immod-1)
    b(immod) = xmod(immod) + (xmod(immod)-xmod(immod-1))/2.0

    c(1) = ymod(1) - (ymod(2)-ymod(1))/2.0
    d(1) = (ymod(1)+ymod(2))/2.0
    DO j = 2, jmmod-1
       c(j) = d(j-1)
       d(j) = (ymod(j)+ymod(j+1))/2.0
    ENDDO
    c(jmmod) = d(jmmod-1)
    d(jmmod) = ymod(jmmod) + (ymod(jmmod)-ymod(jmmod-1))/2.0

    DO i = 1, immod
       DO j = 1, jmmod
          number(i,j) = 0.0
          rugs(i,j) = 0.0
       ENDDO
    ENDDO

    DO ii = 1, immod
       DO jj = 1, jmmod
          DO i = 1, imtmp
             IF( ( xtmp(i)-a(ii) >= 1.e-5.AND.xtmp(i)-b(ii) <= 1.e-5 ).OR. &
                  (   xtmp(i)-a(ii) <= 1.e-5.AND.xtmp(i)-b(ii) >= 1.e-5 )   ) &
                  THEN
                DO j = 1, jmtmp
                   IF ((ytmp(j) - c(jj) >= 1.e-5 &
                        .AND. ytmp(j) - d(jj) <= 1.e-5) .OR. &
                        (ytmp(j) - c(jj) <= 1.e-5 &
                        .AND. ytmp(j) - d(jj) >= 1.e-5)) &
                        THEN
                      number(ii,jj) = number(ii,jj) + 1.0
                      rugs(ii,jj) = rugs(ii,jj) &
                           + LOG(MAX(0.001d0,cham2tmp(i,j)))
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, immod
       DO j = 1, jmmod
          IF (number(i,j) .GT. 0.001) THEN
             rugs(i,j) = rugs(i,j) / number(i,j)
             rugs(i,j) = EXP(rugs(i,j))
          ELSE
             PRINT*, 'probleme,i,j=', i,j
             CALL dist_sphe(xmod(i),ymod(j),xtmp,ytmp,imtmp,jmtmp,distans)
             ij_proche = 1
             zzmin = distans(ij_proche)
             DO ii = 2, imtmp*jmtmp
                IF (distans(ii).LT.zzmin) THEN
                   zzmin = distans(ii)
                   ij_proche = ii
                ENDIF
             ENDDO
             j_proche = (ij_proche-1)/imtmp + 1
             i_proche = ij_proche - (j_proche-1)*imtmp
             PRINT*, "solution:", ij_proche, i_proche, j_proche
             rugs(i,j) = LOG(MAX(0.001d0,cham2tmp(i_proche,j_proche)))
          ENDIF
       ENDDO
    ENDDO

    amin = rugs(1,1)
    AMAX = rugs(1,1)
    DO j = 1, jmmod
       DO i = 1, immod
          IF (rugs(i,j).GT.AMAX) AMAX = rugs(i,j)
          IF (rugs(i,j).LT.amin) amin = rugs(i,j)
       ENDDO
    ENDDO
    PRINT*, 'Ecart-type du modele:', amin, AMAX

    DO j = 1, jmmod
       DO i = 1, immod
          rugs(i,j) = rugs(i,j) / AMAX * 20.0
       ENDDO
    ENDDO

    amin = rugs(1,1)
    AMAX = rugs(1,1)
    DO j = 1, jmmod
       DO i = 1, immod
          IF (rugs(i,j).GT.AMAX) AMAX = rugs(i,j)
          IF (rugs(i,j).LT.amin) amin = rugs(i,j)
       ENDDO
    ENDDO
    PRINT*, 'Longueur de rugosite du modele:', amin, AMAX

  END SUBROUTINE rugsoro
  !
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
