module guide_m

  ! From dyn3d/guide.F,v 1.3 2005/05/25 13:10:09
  ! and dyn3d/guide.h,v 1.1.1.1 2004/05/19 12:53:06

  real tau_min_u,tau_max_u
  real tau_min_v,tau_max_v
  real tau_min_T,tau_max_T
  real tau_min_q,tau_max_q
  real tau_min_P,tau_max_P
  real aire_min,aire_max


  logical guide_u,guide_v,guide_T,guide_Q,guide_P
  real lat_min_guide,lat_max_guide

  LOGICAL ncep,ini_anal
  integer online

contains

  subroutine guide(itau,ucov,vcov,teta,q,masse,ps)

    use dimens_m
    use paramet_m
    use comconst
    use comdissnew
    use comvert
    use conf_gcm_m
    use logic
    use comgeom
    use serre
    use temps
    use tracstoke
    use ener
    use q_sat_m, only: q_sat
    use exner_hyb_m, only: exner_hyb
    use pression_m, only: pression
    use inigrads_m, only: inigrads

    IMPLICIT NONE

    !      ......   Version  du 10/01/98    ..........

    !             avec  coordonnees  verticales hybrides 
    !   avec nouveaux operat. dissipation * ( gradiv2,divgrad2,nxgraro2 )

    !=======================================================================
    !
    !   Auteur:  F.Hourdin
    !   -------
    !
    !   Objet:
    !   ------
    !
    !   GCM LMD nouvelle grille
    !
    !=======================================================================

    !   ...  Dans inigeom , nouveaux calculs pour les elongations  cu , cv 
    !        et possibilite d'appeler une fonction f(y)  a derivee tangente 
    !        hyperbolique a la  place de la fonction a derivee sinusoidale.         

    !   ...  Possibilite de choisir le shema de Van-leer pour l'advection de
    !         q  , en faisant iadv = 10  dans   traceur  (29/04/97) .
    !
    !-----------------------------------------------------------------------
    !   Declarations:
    !   -------------

    include "netcdf.inc"

    !   variables dynamiques
    REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
    REAL teta(ip1jmp1,llm)                 ! temperature potentielle 
    REAL q(ip1jmp1,llm)                 ! temperature potentielle 
    REAL ps(ip1jmp1)                       ! pression  au sol
    REAL masse(ip1jmp1,llm)                ! masse d'air

    !   common passe pour des sorties
    real dxdys(iip1,jjp1),dxdyu(iip1,jjp1),dxdyv(iip1,jjm)
    common/comdxdy/dxdys,dxdyu,dxdyv

    !   variables dynamiques pour les reanalyses.
    REAL ucovrea1(ip1jmp1,llm),vcovrea1(ip1jm,llm) !vts cov reas
    REAL tetarea1(ip1jmp1,llm)             ! temp pot  reales
    REAL qrea1(ip1jmp1,llm)             ! temp pot  reales
    REAL psrea1(ip1jmp1)             ! ps
    REAL ucovrea2(ip1jmp1,llm),vcovrea2(ip1jm,llm) !vts cov reas
    REAL tetarea2(ip1jmp1,llm)             ! temp pot  reales
    REAL qrea2(ip1jmp1,llm)             ! temp pot  reales
    REAL masserea2(ip1jmp1,llm)             ! masse
    REAL psrea2(ip1jmp1)             ! ps

    real alpha_q(ip1jmp1)
    real alpha_T(ip1jmp1),alpha_P(ip1jmp1)
    real alpha_u(ip1jmp1),alpha_v(ip1jm)
    real dday_step,toto,reste,itau_test
    INTEGER step_rea,count_no_rea

    !IM 180305   real aire_min,aire_max
    integer ilon,ilat
    real factt,ztau(ip1jmp1)

    INTEGER, intent(in):: itau
    integer ij, l
    integer ncidpl,varidpl,nlev,status
    integer rcod,rid 
    real ditau,tau,a
    save nlev

    !  TEST SUR QSAT
    real p(ip1jmp1,llmp1),pk(ip1jmp1,llm),pks(ip1jmp1)
    real pkf(ip1jmp1,llm)
    real pres(ip1jmp1,llm)

    real qsat(ip1jmp1,llm)
    real unskap
    real tnat(ip1jmp1,llm)
    !cccccccccccccccc


    LOGICAL first
    save first
    data first/.true./

    save ucovrea1,vcovrea1,tetarea1,psrea1,qrea1
    save ucovrea2,vcovrea2,tetarea2,masserea2,psrea2,qrea2

    save alpha_T,alpha_q,alpha_u,alpha_v,alpha_P,itau_test
    save step_rea,count_no_rea

    character*10 file
    integer igrads
    real dtgrads
    save igrads,dtgrads
    data igrads,dtgrads/2,100./

    print *,'Call sequence information: guide'

    !-----------------------------------------------------------------------
    ! calcul de l'humidite saturante
    !-----------------------------------------------------------------------
    CALL pression( ip1jmp1, ap, bp, ps, p )
    call massdair(p,masse)
    print*,'OK1'
    CALL exner_hyb(ps,p,pks,pk,pkf)
    print*,'OK2'
    tnat(:,:)=pk(:,:)*teta(:,:)/cpp
    print*,'OK3'
    unskap   = 1./ kappa
    pres(:,:)=preff*(pk(:,:)/cpp)**unskap
    print*,'OK4'
    qsat = q_sat(tnat, pres)

    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !   initialisations pour la lecture des reanalyses.
    !    alpha determine la part des injections de donnees a chaque etape
    !    alpha=1 signifie pas d'injection
    !    alpha=0 signifie injection totale
    !-----------------------------------------------------------------------

    print*,'ONLINE=',online
    if(online.eq.-1) then
       return
    endif

    if (first) then

       print*,'initialisation du guide '
       call conf_guide
       print*,'apres conf_guide'

       file='guide'
       call inigrads(igrads &
            ,rlonv,180./pi,-180.,180.,rlatu,-90.,90.,180./pi &
            ,presnivs,1. &
            ,dtgrads,file,'dyn_zon ')

       print* &
            ,'1: en-ligne, 0: hors-ligne (x=x_rea), -1: climat (x=x_gcm)'

       if(online.eq.-1) return
       if (online.eq.1) then

          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !  Constantes de temps de rappel en jour
          !  0.1 c'est en gros 2h30. 
          !  1e10  est une constante infinie donc en gros pas de guidage
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !   coordonnees du centre du zoom
          call coordij(clon,clat,ilon,ilat)
          !   aire de la maille au centre du zoom
          aire_min=aire(ilon+(ilat-1)*iip1)
          !   aire maximale de la maille
          aire_max=0.
          do ij=1,ip1jmp1
             aire_max=max(aire_max,aire(ij))
          enddo
          !  factt = pas de temps en fraction de jour
          factt=dtvr*iperiod/daysec

          !     subroutine tau2alpha(type,im,jm,factt,taumin,taumax,alpha)
          call tau2alpha(3,iip1,jjm ,factt,tau_min_v,tau_max_v,alpha_v)
          call tau2alpha(2,iip1,jjp1,factt,tau_min_u,tau_max_u,alpha_u)
          call tau2alpha(1,iip1,jjp1,factt,tau_min_T,tau_max_T,alpha_T)
          call tau2alpha(1,iip1,jjp1,factt,tau_min_P,tau_max_P,alpha_P)
          call tau2alpha(1,iip1,jjp1,factt,tau_min_q,tau_max_q,alpha_q)

          call dump2d(iip1,jjp1,aire,'AIRE MAILLe ')
          call dump2d(iip1,jjp1,alpha_u,'COEFF U   ')
          call dump2d(iip1,jjp1,alpha_T,'COEFF T   ')

          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !   Cas ou on force exactement par les variables analysees
       else
          alpha_T=0.
          alpha_u=0.
          alpha_v=0.
          alpha_P=0.
          !           physic=.false.
       endif

       itau_test=1001
       step_rea=1
       count_no_rea=0
       ncidpl=-99

       !    itau_test    montre si l'importation a deja ete faite au rang itau
       ! lecture d'un fichier netcdf pour determiner le nombre de niveaux
       if (guide_u) then
          if (ncidpl.eq.-99) ncidpl=NCOPN('u.nc',NCNOWRIT,rcod)
       endif
       !
       if (guide_v) then
          if (ncidpl.eq.-99) ncidpl=NCOPN('v.nc',NCNOWRIT,rcod)
       endif
       !
       if (guide_T) then
          if (ncidpl.eq.-99) ncidpl=NCOPN('T.nc',NCNOWRIT,rcod)
       endif
       !
       if (guide_Q) then
          if (ncidpl.eq.-99) ncidpl=NCOPN('hur.nc',NCNOWRIT,rcod)
       endif
       !
       if (ncep) then
          status=NF_INQ_DIMID(ncidpl,'LEVEL',rid)
       else
          status=NF_INQ_DIMID(ncidpl,'PRESSURE',rid)
       endif
       status=NF_INQ_DIMLEN(ncidpl,rid,nlev)
       print *,'nlev', nlev 
       call ncclos(ncidpl,rcod)
       !   Lecture du premier etat des reanalyses.
       call read_reanalyse(1,ps &
            ,ucovrea2,vcovrea2,tetarea2,qrea2,masserea2,psrea2,1,nlev)
       qrea2(:,:)=max(qrea2(:,:),0.1)


       !-----------------------------------------------------------------------
       !   Debut de l'integration temporelle:
       !   ----------------------------------

    endif ! first
    !
    !-----------------------------------------------------------------------
    !----- IMPORTATION DES VENTS,PRESSION ET TEMPERATURE REELS:
    !-----------------------------------------------------------------------

    ditau=real(itau)
    DDAY_step=real(day_step)
    write(*,*)'ditau,dday_step'
    write(*,*)ditau,dday_step
    toto=4*ditau/dday_step
    reste=toto-aint(toto)
    !     write(*,*)'toto,reste',toto,reste

    if (reste.eq.0.) then
       if (itau_test.eq.itau) then
          write(*,*)'deuxieme passage de advreel a itau=',itau
          stop
       else
          vcovrea1(:,:)=vcovrea2(:,:)
          ucovrea1(:,:)=ucovrea2(:,:)
          tetarea1(:,:)=tetarea2(:,:)
          qrea1(:,:)=qrea2(:,:)

          print*,'LECTURE REANALYSES, pas ',step_rea &
               ,'apres ',count_no_rea,' non lectures'
          step_rea=step_rea+1
          itau_test=itau
          call read_reanalyse(step_rea,ps &
               ,ucovrea2,vcovrea2,tetarea2,qrea2,masserea2,psrea2,1,nlev)
          qrea2(:,:)=max(qrea2(:,:),0.1)
          factt=dtvr*iperiod/daysec
          ztau(:)=factt/max(alpha_T(:),1.e-10)
          call wrgrads(igrads,1,aire   ,'aire      ','aire      ' )
          call wrgrads(igrads,1,dxdys  ,'dxdy      ','dxdy      ' )
          call wrgrads(igrads,1,alpha_u,'au        ','au        ' )
          call wrgrads(igrads,1,alpha_T,'at        ','at        ' )
          call wrgrads(igrads,1,ztau,'taut      ','taut      ' )
          call wrgrads(igrads,llm,ucov,'u         ','u         ' )
          call wrgrads(igrads,llm,ucovrea2,'ua        ','ua        ' )
          call wrgrads(igrads,llm,teta,'T         ','T         ' )
          call wrgrads(igrads,llm,tetarea2,'Ta        ','Ta        ' )
          call wrgrads(igrads,llm,qrea2,'Qa        ','Qa        ' )
          call wrgrads(igrads,llm,q,'Q         ','Q         ' )

          call wrgrads(igrads,llm,qsat,'QSAT      ','QSAT      ' )

       endif
    else
       count_no_rea=count_no_rea+1
    endif

    !-----------------------------------------------------------------------
    !   Guidage
    !    x_gcm = a * x_gcm + (1-a) * x_reanalyses
    !-----------------------------------------------------------------------

    if(ini_anal) print*,'ATTENTION !!! ON PART DU GUIDAGE'

    ditau=real(itau)
    dday_step=real(day_step)


    tau=4*ditau/dday_step
    tau=tau-aint(tau)

    !  ucov
    if (guide_u) then
       do l=1,llm
          do ij=1,ip1jmp1
             a=(1.-tau)*ucovrea1(ij,l)+tau*ucovrea2(ij,l)
             ucov(ij,l)=(1.-alpha_u(ij))*ucov(ij,l)+alpha_u(ij)*a
             if (first.and.ini_anal) ucov(ij,l)=a
          enddo
       enddo
    endif

    !  teta
    if (guide_T) then
       do l=1,llm
          do ij=1,ip1jmp1
             a=(1.-tau)*tetarea1(ij,l)+tau*tetarea2(ij,l)
             teta(ij,l)=(1.-alpha_T(ij))*teta(ij,l)+alpha_T(ij)*a
             if (first.and.ini_anal) teta(ij,l)=a
          enddo
       enddo
    endif

    !  P
    if (guide_P) then
       do ij=1,ip1jmp1
          a=(1.-tau)*psrea1(ij)+tau*psrea2(ij)
          ps(ij)=(1.-alpha_P(ij))*ps(ij)+alpha_P(ij)*a
          if (first.and.ini_anal) ps(ij)=a
       enddo
       CALL pression(ip1jmp1,ap,bp,ps,p)
       CALL massdair(p,masse)
    endif


    !  q
    if (guide_Q) then
       do l=1,llm
          do ij=1,ip1jmp1
             a=(1.-tau)*qrea1(ij,l)+tau*qrea2(ij,l)
             !   hum relative en % -> hum specif
             a=qsat(ij,l)*a*0.01
             q(ij,l)=(1.-alpha_Q(ij))*q(ij,l)+alpha_Q(ij)*a
             if (first.and.ini_anal) q(ij,l)=a
          enddo
       enddo
    endif

    ! vcov
    if (guide_v) then
       do l=1,llm
          do ij=1,ip1jm
             a=(1.-tau)*vcovrea1(ij,l)+tau*vcovrea2(ij,l)
             vcov(ij,l)=(1.-alpha_v(ij))*vcov(ij,l)+alpha_v(ij)*a
             if (first.and.ini_anal) vcov(ij,l)=a
          enddo
          if (first.and.ini_anal) vcov(ij,l)=a
       enddo
    endif

    !     call dump2d(iip1,jjp1,tetarea1,'TETA REA 1     ')
    !     call dump2d(iip1,jjp1,tetarea2,'TETA REA 2     ')
    !     call dump2d(iip1,jjp1,teta,'TETA           ')

    first=.false.

    return
  end subroutine guide

  !=======================================================================
  subroutine tau2alpha(type,pim,pjm,factt,taumin,taumax,alpha)
    !=======================================================================

    use dimens_m
    use paramet_m
    use comconst, only: pi
    use comgeom
    use serre
    implicit none

    !   arguments :
    integer type
    integer pim,pjm
    real factt,taumin,taumax
    real dxdy_,alpha(pim,pjm)
    real dxdy_min,dxdy_max

    !  local :
    real alphamin,alphamax,gamma,xi
    save gamma
    integer i,j,ilon,ilat

    logical first
    save first
    data first/.true./

    real zdx(iip1,jjp1),zdy(iip1,jjp1)

    real zlat
    real dxdys(iip1,jjp1),dxdyu(iip1,jjp1),dxdyv(iip1,jjm)
    common/comdxdy/dxdys,dxdyu,dxdyv

    if (first) then
       do j=2,jjm
          do i=2,iip1
             zdx(i,j)=0.5*(cu_2d(i-1,j)+cu_2d(i,j))/cos(rlatu(j))
          enddo
          zdx(1,j)=zdx(iip1,j)
       enddo
       do j=2,jjm
          do i=1,iip1
             zdy(i,j)=0.5*(cv_2d(i,j-1)+cv_2d(i,j))
          enddo
       enddo
       do i=1,iip1
          zdx(i,1)=zdx(i,2)
          zdx(i,jjp1)=zdx(i,jjm)
          zdy(i,1)=zdy(i,2)
          zdy(i,jjp1)=zdy(i,jjm)
       enddo
       do j=1,jjp1
          do i=1,iip1
             dxdys(i,j)=sqrt(zdx(i,j)*zdx(i,j)+zdy(i,j)*zdy(i,j))
          enddo
       enddo
       do j=1,jjp1
          do i=1,iim
             dxdyu(i,j)=0.5*(dxdys(i,j)+dxdys(i+1,j))
          enddo
          dxdyu(iip1,j)=dxdyu(1,j)
       enddo
       do j=1,jjm
          do i=1,iip1
             dxdyv(i,j)=0.5*(dxdys(i,j)+dxdys(i+1,j))
          enddo
       enddo

       call dump2d(iip1,jjp1,dxdys,'DX2DY2 SCAL  ')
       call dump2d(iip1,jjp1,dxdyu,'DX2DY2 U     ')
       call dump2d(iip1,jjp1,dxdyv,'DX2DY2 v     ')

       !   coordonnees du centre du zoom
       call coordij(clon,clat,ilon,ilat)
       !   aire de la maille au centre du zoom
       dxdy_min=dxdys(ilon,ilat)
       !   dxdy maximale de la maille
       dxdy_max=0.
       do j=1,jjp1
          do i=1,iip1
             dxdy_max=max(dxdy_max,dxdys(i,j))
          enddo
       enddo

       if (abs(grossismx-1.).lt.0.1.or.abs(grossismy-1.).lt.0.1) then
          print*,'ATTENTION modele peu zoome'
          print*,'ATTENTION on prend une constante de guidage cste'
          gamma=0.
       else
          gamma=(dxdy_max-2.*dxdy_min)/(dxdy_max-dxdy_min)
          print*,'gamma=',gamma
          if (gamma.lt.1.e-5) then
             print*,'gamma =',gamma,'<1e-5'
             stop
          endif
          print*,'gamma=',gamma
          gamma=log(0.5)/log(gamma)
       endif
    endif

    alphamin=factt/taumax
    alphamax=factt/taumin

    do j=1,pjm
       do i=1,pim
          if (type.eq.1) then
             dxdy_=dxdys(i,j)
             zlat=rlatu(j)*180./pi
          elseif (type.eq.2) then
             dxdy_=dxdyu(i,j)
             zlat=rlatu(j)*180./pi
          elseif (type.eq.3) then
             dxdy_=dxdyv(i,j)
             zlat=rlatv(j)*180./pi
          endif
          if (abs(grossismx-1.).lt.0.1.or.abs(grossismy-1.).lt.0.1) then
             !  pour une grille reguliere, xi=xxx**0=1 -> alpha=alphamin
             alpha(i,j)=alphamin
          else
             xi=((dxdy_max-dxdy_)/(dxdy_max-dxdy_min))**gamma
             xi=min(xi,1.)
             if(lat_min_guide.le.zlat .and. zlat.le.lat_max_guide) then
                alpha(i,j)=xi*alphamin+(1.-xi)*alphamax
             else
                alpha(i,j)=0.
             endif
          endif
       enddo
    enddo


    return
  end subroutine tau2alpha

end module guide_m
