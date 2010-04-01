!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/conf_guide.F,v 1.1.1.1 2004/05/19 12:53:07 lmdzadmin Exp $
!
c
c
      SUBROUTINE conf_guide
c
      use getparam
      use guide_m
      IMPLICIT NONE

c-----------------------------------------------------------------------
c  Parametres de controle du run:
c-----------------------------------------------------------------------

      call getpar('guide.eff')

      call getpar('online',1,online,'Index de controle du guide')
      CALL getpar('ncep',.false.,ncep,'Coordonnee vert NCEP ou ECMWF')
      CALL getpar('ini_anal',.false.,ini_anal,'Initial = analyse')

      CALL getpar('guide_u',.true.,guide_u,'guidage de u')
      CALL getpar('guide_v',.true.,guide_v,'guidage de v')
      CALL getpar('guide_T',.true.,guide_T,'guidage de T')
      CALL getpar('guide_P',.true.,guide_P,'guidage de P')
      CALL getpar('guide_Q',.true.,guide_Q,'guidage de Q')

c   Constantes de rappel. Unite : fraction de jour
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

c   Latitude min et max pour le rappel.
c   dans le cas ou on 'a les analyses que sur une bande de latitudes.
      CALL getpar('lat_min_guide',-90.,lat_min_guide
     s     ,'Latitude minimum pour le guidage ')
      CALL getpar('lat_max_guide', 90.,lat_max_guide
     s     ,'Latitude maximum pour le guidage ')


      CALL getpar

      end
