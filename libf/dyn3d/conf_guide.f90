module conf_guide_m

  IMPLICIT NONE

  REAL tau_min_u, tau_max_u
  REAL tau_min_v, tau_max_v
  REAL tau_min_t, tau_max_t
  REAL tau_min_q, tau_max_q
  REAL tau_min_p, tau_max_p
  LOGICAL ncep, ini_anal
  LOGICAL guide_u, guide_v, guide_t, guide_q, guide_p
  INTEGER online

contains

  SUBROUTINE conf_guide

    ! From LMDZ4/libf/dyn3d/conf_guide.F, version 1.1.1.1 2004/05/19 12:53:07
    !  Parametres de controle du run:

    use getparam, only: ini_getparam, getpar, fin_getparam
    use tau2alpha_m, only: lat_max_guide, lat_min_guide

    !-----------------------------------------------------------------------

    print *, "Call sequence information: conf_guide"
    call ini_getparam('guide.eff')

    call getpar('online',1,online,'Index de controle du guide')
    CALL getpar('ncep',.false.,ncep,'Coordonnee vert NCEP ou ECMWF')
    CALL getpar('ini_anal',.false.,ini_anal,'Initial = analyse')

    CALL getpar('guide_u',.true.,guide_u,'guidage de u')
    CALL getpar('guide_v',.true.,guide_v,'guidage de v')
    CALL getpar('guide_T',.true.,guide_T,'guidage de T')
    CALL getpar('guide_P',.true.,guide_P,'guidage de P')
    CALL getpar('guide_Q',.true.,guide_Q,'guidage de Q')

    !   Constantes de rappel. Unite : fraction de jour
    CALL getpar('tau_min_u',0.02,tau_min_u,'Cste de rappel min, u')
    CALL getpar('tau_max_u', 10.,tau_max_u,'Cste de rappel max, u')
    CALL getpar('tau_min_v',0.02,tau_min_v,'Cste de rappel min, v')
    CALL getpar('tau_max_v', 10.,tau_max_v,'Cste de rappel max, v')
    CALL getpar('tau_min_T',0.02,tau_min_T,'Cste de rappel min, T')
    CALL getpar('tau_max_T', 10.,tau_max_T,'Cste de rappel max, T')
    CALL getpar('tau_min_Q',0.02,tau_min_Q,'Cste de rappel min, Q')
    CALL getpar('tau_max_Q', 10.,tau_max_Q,'Cste de rappel max, Q')
    CALL getpar('tau_min_P',0.02,tau_min_P,'Cste de rappel min, P')
    CALL getpar('tau_max_P', 10.,tau_max_P,'Cste de rappel max, P')

    !   Latitude min et max pour le rappel.
    !   dans le cas ou on 'a les analyses que sur une bande de latitudes.
    CALL getpar('lat_min_guide',-90.,lat_min_guide &
         ,'Latitude minimum pour le guidage ')
    CALL getpar('lat_max_guide', 90.,lat_max_guide &
         ,'Latitude maximum pour le guidage ')

    CALL fin_getparam

  end SUBROUTINE conf_guide

end module conf_guide_m
