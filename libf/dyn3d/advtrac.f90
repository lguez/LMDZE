SUBROUTINE advtrac(pbaru,pbarv,p,masse,q,iapptrac,teta,pk)

  ! From dyn3d/advtrac.F,v 1.4 2005/04/13 08:58:34
  !     Auteur :  F. Hourdin

  !     Modif. P. Le Van     (20/12/97)
  !            F. Codron     (10/99)
  !            D. Le Croller (07/2001)
  !            M.A Filiberti (04/2002)

  USE dimens_m
  USE paramet_m
  USE comconst
  USE comvert
  USE conf_gcm_m
  USE logic
  USE comgeom
  USE temps
  USE ener
  USE advtrac_m
  USE comdissip
  IMPLICIT NONE


  !-------------------------------------------------------------------
  !     Arguments
  !-------------------------------------------------------------------
  !     Ajout PPM
  !--------------------------------------------------------
  REAL massebx(ip1jmp1,llm), masseby(ip1jm,llm)
  !--------------------------------------------------------
  INTEGER iapptrac
  REAL pbaru(ip1jmp1,llm), pbarv(ip1jm,llm)
  REAL q(ip1jmp1,llm,nqmx), masse(ip1jmp1,llm)
  REAL, intent(in):: p(ip1jmp1,llmp1)
  real teta(ip1jmp1,llm)
  REAL pk(ip1jmp1,llm)

  !-------------------------------------------------------------
  !     Variables locales
  !-------------------------------------------------------------

  REAL pbaruc(ip1jmp1,llm), pbarvc(ip1jm,llm)
  REAL massem(ip1jmp1,llm), zdp(ip1jmp1)
  REAL pbarug(ip1jmp1,llm), pbarvg(ip1jm,llm), wg(ip1jmp1,llm)
  REAL cpuadv(nqmx)
  COMMON /cpuadv/cpuadv

  INTEGER iadvtr
  INTEGER ij, l, iq
  REAL zdpmin, zdpmax
  EXTERNAL minmax
  SAVE iadvtr, massem, pbaruc, pbarvc
  DATA iadvtr/0/
  !----------------------------------------------------------
  !     Rajouts pour PPM
  !----------------------------------------------------------
  INTEGER indice, n
  ! Pas de temps adaptatif pour que CFL<1                
  REAL dtbon
  REAL cflmaxz, aaa, & ! CFL maximum                                
       bbb
  REAL psppm(iim,jjp1) ! pression  au sol                           
  REAL unatppm(iim,jjp1,llm), vnatppm(iim,jjp1,llm)
  REAL qppm(iim*jjp1,llm,nqmx)
  REAL fluxwppm(iim,jjp1,llm)
  REAL apppm(llmp1), bpppm(llmp1)
  LOGICAL dum, fill
  DATA fill/ .TRUE./
  DATA dum/ .TRUE./


  IF (iadvtr==0) THEN
     CALL initial0(ijp1llm,pbaruc)
     CALL initial0(ijmllm,pbarvc)
  END IF

  !   accumulation des flux de masse horizontaux
  DO l = 1, llm
     DO ij = 1, ip1jmp1
        pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
     END DO
     DO ij = 1, ip1jm
        pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
     END DO
  END DO

  !   selection de la masse instantannee des mailles avant le transport.
  IF (iadvtr==0) THEN

     CALL scopy(ip1jmp1*llm,masse,1,massem,1)
     !cc         CALL filtreg ( massem ,jjp1, llm,-2, 2, .TRUE., 1 )

  END IF

  iadvtr = iadvtr + 1
  iapptrac = iadvtr


  !   Test pour savoir si on advecte a ce pas de temps
  IF (iadvtr==iapp_tracvl) THEN

     !c   ..  Modif P.Le Van  ( 20/12/97 )  ....
     !c

     !   traitement des flux de masse avant advection.
     !     1. calcul de w
     !     2. groupement des mailles pres du pole.

     CALL groupe(massem,pbaruc,pbarvc,pbarug,pbarvg,wg)


     !  test sur l'eventuelle creation de valeurs negatives de la masse
     DO l = 1, llm - 1
        DO ij = iip2 + 1, ip1jm
           zdp(ij) = pbarug(ij-1,l) - pbarug(ij,l) - pbarvg(ij-iip1,l) + &
                pbarvg(ij,l) + wg(ij,l+1) - wg(ij,l)
        END DO
        CALL scopy(jjm-1,zdp(iip1+iip1),iip1,zdp(iip2),iip1)
        DO ij = iip2, ip1jm
           zdp(ij) = zdp(ij)*dtvr/massem(ij,l)
        END DO


        CALL minmax(ip1jm-iip1,zdp(iip2),zdpmin,zdpmax)

        IF (max(abs(zdpmin),abs(zdpmax))>0.5) THEN
           PRINT *, 'WARNING DP/P l=', l, '  MIN:', zdpmin, '   MAX:', zdpmax
        END IF

     END DO

     !-------------------------------------------------------------------
     !   Advection proprement dite (Modification Le Croller (07/2001)
     !-------------------------------------------------------------------

     !----------------------------------------------------
     !        Calcul des moyennes basées sur la masse
     !----------------------------------------------------
     CALL massbar(massem,massebx,masseby)

     !-----------------------------------------------------------
     !     Appel des sous programmes d'advection
     !-----------------------------------------------------------
     DO iq = 1, nqmx
        !        call clock(t_initial)
        IF (iadv(iq)==0) CYCLE
        !   ----------------------------------------------------------------
        !   Schema de Van Leer I MUSCL
        !   ----------------------------------------------------------------
        IF (iadv(iq)==10) THEN
           CALL vlsplt(q(1,1,iq),2.,massem,wg,pbarug,pbarvg,dtvr)


           !   ----------------------------------------------------------------
           !   Schema "pseudo amont" + test sur humidite specifique
           !    pour la vapeur d'eau. F. Codron
           !   ----------------------------------------------------------------
        ELSE IF (iadv(iq)==14) THEN

           CALL vlspltqs(q(1,1,1),2.,massem,wg,pbarug,pbarvg,dtvr,p,pk,teta)
           !   ----------------------------------------------------------------
           !   Schema de Frederic Hourdin
           !   ----------------------------------------------------------------
        ELSE IF (iadv(iq)==12) THEN
           !            Pas de temps adaptatif
           CALL adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           IF (n>1) THEN
              WRITE (*,*) 'WARNING horizontal dt=', dtbon, 'dtvr=', dtvr, &
                   'n=', n
           END IF
           DO indice = 1, n
              CALL advn(q(1,1,iq),massem,wg,pbarug,pbarvg,dtbon,1)
           END DO
        ELSE IF (iadv(iq)==13) THEN
           !            Pas de temps adaptatif
           CALL adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           IF (n>1) THEN
              WRITE (*,*) 'WARNING horizontal dt=', dtbon, 'dtvr=', dtvr, &
                   'n=', n
           END IF
           DO indice = 1, n
              CALL advn(q(1,1,iq),massem,wg,pbarug,pbarvg,dtbon,2)
           END DO
           !   ----------------------------------------------------------------
           !   Schema de pente SLOPES
           !   ----------------------------------------------------------------
        ELSE IF (iadv(iq)==20) THEN
           CALL pentes_ini(q(1,1,iq),wg,massem,pbarug,pbarvg,0)

           !   ----------------------------------------------------------------
           !   Schema de Prather
           !   ----------------------------------------------------------------
        ELSE IF (iadv(iq)==30) THEN
           !            Pas de temps adaptatif
           CALL adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           IF (n>1) THEN
              WRITE (*,*) 'WARNING horizontal dt=', dtbon, 'dtvr=', dtvr, &
                   'n=', n
           END IF
           CALL prather(q(1,1,iq),wg,massem,pbarug,pbarvg,n,dtbon)
           !   ----------------------------------------------------------------
           !   Schemas PPM Lin et Rood
           !   ----------------------------------------------------------------
        ELSE IF (iadv(iq)==11 .OR. (iadv(iq)>=16 .AND. iadv(iq)<=18)) THEN

           !        Test sur le flux horizontal
           !        Pas de temps adaptatif
           CALL adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           IF (n>1) THEN
              WRITE (*,*) 'WARNING horizontal dt=', dtbon, 'dtvr=', dtvr, &
                   'n=', n
           END IF
           !        Test sur le flux vertical
           cflmaxz = 0.
           DO l = 2, llm
              DO ij = iip2, ip1jm
                 aaa = wg(ij,l)*dtvr/massem(ij,l)
                 cflmaxz = max(cflmaxz,aaa)
                 bbb = -wg(ij,l)*dtvr/massem(ij,l-1)
                 cflmaxz = max(cflmaxz,bbb)
              END DO
           END DO
           IF (cflmaxz>=1) THEN
              WRITE (*,*) 'WARNING vertical', 'CFLmaxz=', cflmaxz
           END IF

           !-----------------------------------------------------------
           !        Ss-prg interface LMDZ.4->PPM3d
           !-----------------------------------------------------------

           CALL interpre(q(1,1,iq),qppm(1,1,iq),wg,fluxwppm,massem,apppm, &
                bpppm,massebx,masseby,pbarug,pbarvg,unatppm,vnatppm,psppm)

           DO indice = 1, n
              !----------------------------------------------------------------
              !                         VL (version PPM) horiz. et PPM vert.
              !---------------------------------------------------------------
              IF (iadv(iq)==11) THEN
                 !                  Ss-prg PPM3d de Lin
                 CALL ppm3d(1,qppm(1,1,iq),psppm,psppm,unatppm,vnatppm, &
                      fluxwppm,dtbon,2,2,2,1,iim,jjp1,2,llm,apppm,bpppm,0.01, &
                      6400000,fill,dum,220.)

                 !-----------------------------------------------------------
                 !                           Monotonic PPM
                 !-------------------------------------------------------
              ELSE IF (iadv(iq)==16) THEN
                 !                  Ss-prg PPM3d de Lin
                 CALL ppm3d(1,qppm(1,1,iq),psppm,psppm,unatppm,vnatppm, &
                      fluxwppm,dtbon,3,3,3,1,iim,jjp1,2,llm,apppm,bpppm,0.01, &
                      6400000,fill,dum,220.)
                 !                           Semi Monotonic PPM
              ELSE IF (iadv(iq)==17) THEN
                 !                  Ss-prg PPM3d de Lin
                 CALL ppm3d(1,qppm(1,1,iq),psppm,psppm,unatppm,vnatppm, &
                      fluxwppm,dtbon,4,4,4,1,iim,jjp1,2,llm,apppm,bpppm,0.01, &
                      6400000,fill,dum,220.)
                 !                         Positive Definite PPM
              ELSE IF (iadv(iq)==18) THEN
                 !                  Ss-prg PPM3d de Lin
                 CALL ppm3d(1,qppm(1,1,iq),psppm,psppm,unatppm,vnatppm, &
                      fluxwppm,dtbon,5,5,5,1,iim,jjp1,2,llm,apppm,bpppm,0.01, &
                      6400000,fill,dum,220.)
              END IF
           END DO
           !-----------------------------------------------------------------
           !               Ss-prg interface PPM3d-LMDZ.4
           !-----------------------------------------------------------------
           CALL interpost(q(1,1,iq),qppm(1,1,iq))
        END IF
        !----------------------------------------------------------------------

        !-----------------------------------------------------------------
        ! On impose une seule valeur du traceur au pôle Sud j=jjm+1=jjp1
        ! et Nord j=1
        !-----------------------------------------------------------------

        !                  call traceurpole(q(1,1,iq),massem)

        ! calcul du temps cpu pour un schema donne

        !                  call clock(t_final)
        !ym                  tps_cpu=t_final-t_initial
        !ym                  cpuadv(iq)=cpuadv(iq)+tps_cpu

     END DO


     !------------------------------------------------------------------
     !   on reinitialise a zero les flux de masse cumules
     !---------------------------------------------------
     iadvtr = 0


  END IF
  ! if iadvtr.EQ.iapp_tracvl                                 
  RETURN
END SUBROUTINE advtrac
