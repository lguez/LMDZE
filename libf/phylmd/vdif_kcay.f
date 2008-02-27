!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/vdif_kcay.F,v 1.1 2004/06/22 11:45:36 lmdzadmin Exp $
!
      SUBROUTINE vdif_kcay(ngrid,dt,g,rconst,plev,temp
     s   ,zlev,zlay,u,v,teta,cd,q2,q2diag,km,kn,ustar
     s   ,l_mix)
      use dimens_m
      use dimphy
      IMPLICIT NONE
c.......................................................................
c.......................................................................
c
c dt : pas de temps
c g  : g
c zlev : altitude a chaque niveau (interface inferieure de la couche
c        de meme indice)
c zlay : altitude au centre de chaque couche
c u,v : vitesse au centre de chaque couche
c       (en entree : la valeur au debut du pas de temps)
c teta : temperature potentielle au centre de chaque couche
c        (en entree : la valeur au debut du pas de temps)
c cd : cdrag
c      (en entree : la valeur au debut du pas de temps)
c q2 : $q^2$ au bas de chaque couche
c      (en entree : la valeur au debut du pas de temps)
c      (en sortie : la valeur a la fin du pas de temps)
c km : diffusivite turbulente de quantite de mouvement (au bas de chaque
c      couche)
c      (en sortie : la valeur a la fin du pas de temps)
c kn : diffusivite turbulente des scalaires (au bas de chaque couche)
c      (en sortie : la valeur a la fin du pas de temps)
c 
c.......................................................................
      REAL dt,g,rconst
      real plev(klon,klev+1),temp(klon,klev)
      real ustar(klon),snstable
      REAL zlev(klon,klev+1)
      REAL zlay(klon,klev)
      REAL u(klon,klev)
      REAL v(klon,klev)
      REAL teta(klon,klev)
      REAL cd(klon)
      REAL q2(klon,klev+1),q2s(klon,klev+1)
      REAL q2diag(klon,klev+1)
      REAL km(klon,klev+1)
      REAL kn(klon,klev+1)
      real sq(klon),sqz(klon),zz(klon,klev+1),zq,long0(klon)

      integer l_mix,iii
c.......................................................................
c
c nlay : nombre de couches        
c nlev : nombre de niveaux
c ngrid : nombre de points de grille       
c unsdz : 1 sur l'epaisseur de couche
c unsdzdec : 1 sur la distance entre le centre de la couche et le
c            centre de la couche inferieure
c q : echelle de vitesse au bas de chaque couche
c     (valeur a la fin du pas de temps)
c
c.......................................................................
      INTEGER nlay,nlev,ngrid
      REAL unsdz(klon,klev)
      REAL unsdzdec(klon,klev+1)
      REAL q(klon,klev+1)

c.......................................................................
c
c kmpre : km au debut du pas de temps
c qcstat : q : solution stationnaire du probleme couple
c          (valeur a la fin du pas de temps)
c q2cstat : q2 : solution stationnaire du probleme couple
c           (valeur a la fin du pas de temps)
c
c.......................................................................
      REAL kmpre(klon,klev+1)
      REAL qcstat
      REAL q2cstat
      real sss,sssq
c.......................................................................
c
c long : longueur de melange calculee selon Blackadar
c
c.......................................................................
      REAL long(klon,klev+1)
c.......................................................................
c
c kmq3 : terme en q^3 dans le developpement de km
c        (valeur au debut du pas de temps)
c kmcstat : valeur de km solution stationnaire du systeme {q2 ; du/dz}
c           (valeur a la fin du pas de temps)
c knq3 : terme en q^3 dans le developpement de kn
c mcstat : valeur de m solution stationnaire du systeme {q2 ; du/dz}
c          (valeur a la fin du pas de temps)
c m2cstat : valeur de m2 solution stationnaire du systeme {q2 ; du/dz}
c           (valeur a la fin du pas de temps)
c m : valeur a la fin du pas de temps
c mpre : valeur au debut du pas de temps
c m2 : valeur a la fin du pas de temps
c n2 : valeur a la fin du pas de temps
c 
c.......................................................................
      REAL kmq3
      REAL kmcstat
      REAL knq3
      REAL mcstat
      REAL m2cstat
      REAL m(klon,klev+1)
      REAL mpre(klon,klev+1)
      REAL m2(klon,klev+1)
      REAL n2(klon,klev+1)
c.......................................................................
c
c gn : intermediaire pour les coefficients de stabilite
c gnmin : borne inferieure de gn (-0.23 ou -0.28)
c gnmax : borne superieure de gn (0.0233)
c gninf : vrai si gn est en dessous de sa borne inferieure
c gnsup : vrai si gn est en dessus de sa borne superieure
c gm : drole d'objet bien utile
c ri : nombre de Richardson
c sn : coefficient de stabilite pour n
c snq2 : premier terme du developement limite de sn en q2
c sm : coefficient de stabilite pour m
c smq2 : premier terme du developement limite de sm en q2
c
c.......................................................................
      REAL gn
      REAL gnmin
      REAL gnmax
      LOGICAL gninf
      LOGICAL gnsup
      REAL gm
c      REAL ri(klon,klev+1)
      REAL sn(klon,klev+1)
      REAL snq2(klon,klev+1)
      REAL sm(klon,klev+1)
      REAL smq2(klon,klev+1)
c.......................................................................
c
c kappa : consatnte de Von Karman (0.4)
c long00 : longueur de reference pour le calcul de long (160)
c a1,a2,b1,b2,c1 : constantes d'origine pour les  coefficients
c                  de stabilite (0.92/0.74/16.6/10.1/0.08)
c cn1,cn2 : constantes pour sn
c cm1,cm2,cm3,cm4 : constantes pour sm
c
c.......................................................................
      REAL kappa
      REAL long00
      REAL a1,a2,b1,b2,c1
      REAL cn1,cn2
      REAL cm1,cm2,cm3,cm4
c.......................................................................
c
c termq : termes en $q$ dans l'equation de q2
c termq3 : termes en $q^3$ dans l'equation de q2
c termqm2 : termes en $q*m^2$ dans l'equation de q2
c termq3m2 : termes en $q^3*m^2$ dans l'equation de q2
c
c.......................................................................
      REAL termq
      REAL termq3
      REAL termqm2
      REAL termq3m2
c.......................................................................
c
c q2min : borne inferieure de q2
c q2max : borne superieure de q2
c
c.......................................................................
      REAL q2min
      REAL q2max
c.......................................................................
c knmin : borne inferieure de kn
c kmmin : borne inferieure de km
c.......................................................................
      REAL knmin
      REAL kmmin
c.......................................................................
      INTEGER ilay,ilev,igrid
      REAL tmp1,tmp2
c.......................................................................
      PARAMETER (kappa=0.4E+0)
      PARAMETER (long00=160.E+0)
c     PARAMETER (gnmin=-10.E+0)
      PARAMETER (gnmin=-0.28)
      PARAMETER (gnmax=0.0233E+0)
      PARAMETER (a1=0.92E+0)
      PARAMETER (a2=0.74E+0)
      PARAMETER (b1=16.6E+0)
      PARAMETER (b2=10.1E+0)
      PARAMETER (c1=0.08E+0)
      PARAMETER (knmin=1.E-5)
      PARAMETER (kmmin=1.E-5)
      PARAMETER (q2min=1.e-5)
      PARAMETER (q2max=1.E+2)
      PARAMETER (nlay=klev)
      PARAMETER (nlev=klev+1)
c
      PARAMETER (
     &  cn1=a2*(1.E+0 -6.E+0 *a1/b1)
     &          )
      PARAMETER (
     &  cn2=-3.E+0 *a2*(6.E+0 *a1+b2)
     &          )
      PARAMETER (
     &  cm1=a1*(1.E+0 -3.E+0 *c1-6.E+0 *a1/b1)
     &          )
      PARAMETER (
     &  cm2=a1*(-3.E+0 *a2*((b2-3.E+0 *a2)*(1.E+0 -6.E+0 *a1/b1)
     &          -3.E+0 *c1*(b2+6.E+0 *a1)))
     &          )
      PARAMETER (
     &  cm3=-3.E+0 *a2*(6.E+0 *a1+b2)
     &          )
      PARAMETER (
     &  cm4=-9.E+0 *a1*a2
     &          )

      logical first
      save first
      data first/.true./
c.......................................................................
c  traitment des valeur de q2 en entree
c.......................................................................
c
c   Initialisation de q2

      call yamada(ngrid,dt,g,rconst,plev,temp
     s   ,zlev,zlay,u,v,teta,cd,q2diag,km,kn,ustar
     s   ,l_mix)
      if (first.and.1.eq.1) then
      first=.false.
      q2=q2diag
      endif

      DO ilev=1,nlev
                                                      DO igrid=1,ngrid 
        q2(igrid,ilev)=amax1(q2(igrid,ilev),q2min)
        q(igrid,ilev)=sqrt(q2(igrid,ilev))
                                                      ENDDO
      ENDDO
c
                                                      DO igrid=1,ngrid 
      tmp1=cd(igrid)*(u(igrid,1)**2+v(igrid,1)**2)
      q2(igrid,1)=b1**(2.E+0/3.E+0)*tmp1
      q2(igrid,1)=amax1(q2(igrid,1),q2min)
      q(igrid,1)=sqrt(q2(igrid,1))
                                                      ENDDO
c
c.......................................................................
c  les increments verticaux
c.......................................................................
c
c!!!!! allerte !!!!!c
c!!!!! zlev n'est pas declare a nlev !!!!!c
c!!!!! ---->
                                                      DO igrid=1,ngrid 
            zlev(igrid,nlev)=zlay(igrid,nlay)
     &             +( zlay(igrid,nlay) - zlev(igrid,nlev-1) )
                                                      ENDDO            
c!!!!! <----
c!!!!! allerte !!!!!c
c
      DO ilay=1,nlay
                                                      DO igrid=1,ngrid 
        unsdz(igrid,ilay)=1.E+0/(zlev(igrid,ilay+1)-zlev(igrid,ilay))
                                                      ENDDO
      ENDDO
                                                      DO igrid=1,ngrid 
      unsdzdec(igrid,1)=1.E+0/(zlay(igrid,1)-zlev(igrid,1))
                                                      ENDDO
      DO ilay=2,nlay
                                                      DO igrid=1,ngrid 
        unsdzdec(igrid,ilay)=1.E+0/(zlay(igrid,ilay)-zlay(igrid,ilay-1))
                                                      ENDDO
      ENDDO
                                                      DO igrid=1,ngrid 
      unsdzdec(igrid,nlay+1)=1.E+0/(zlev(igrid,nlay+1)-zlay(igrid,nlay))
                                                      ENDDO
c
c.......................................................................
c  le cisaillement et le gradient de temperature
c.......................................................................
c
                                                      DO igrid=1,ngrid 
      m2(igrid,1)=(unsdzdec(igrid,1)
     &                   *u(igrid,1))**2
     &                 +(unsdzdec(igrid,1)
     &                   *v(igrid,1))**2
      m(igrid,1)=sqrt(m2(igrid,1))
      mpre(igrid,1)=m(igrid,1)
                                                      ENDDO
c
c-----------------------------------------------------------------------
      DO ilev=2,nlev-1
                                                      DO igrid=1,ngrid 
c-----------------------------------------------------------------------
c
        n2(igrid,ilev)=g*unsdzdec(igrid,ilev)
     &                   *(teta(igrid,ilev)-teta(igrid,ilev-1))
     &                   /(teta(igrid,ilev)+teta(igrid,ilev-1)) *2.E+0
c       n2(igrid,ilev)=0.
c
c --->
c       on ne sais traiter que les cas stratifies. et l'ajustement
c       convectif est cense faire en sorte que seul des configurations
c       stratifiees soient rencontrees en entree de cette routine.
c       mais, bon ... on sait jamais (meme on sait que n2 prends
c       quelques valeurs negatives ... parfois) alors : 
c<---
c
        IF (n2(igrid,ilev).lt.0.E+0) THEN
          n2(igrid,ilev)=0.E+0
        ENDIF
c
        m2(igrid,ilev)=(unsdzdec(igrid,ilev)
     &                     *(u(igrid,ilev)-u(igrid,ilev-1)))**2
     &                   +(unsdzdec(igrid,ilev)
     &                     *(v(igrid,ilev)-v(igrid,ilev-1)))**2
        m(igrid,ilev)=sqrt(m2(igrid,ilev))
        mpre(igrid,ilev)=m(igrid,ilev)
c
c-----------------------------------------------------------------------
                                                      ENDDO
      ENDDO
c-----------------------------------------------------------------------
c
                                                      DO igrid=1,ngrid 
      m2(igrid,nlev)=m2(igrid,nlev-1)
      m(igrid,nlev)=m(igrid,nlev-1)
      mpre(igrid,nlev)=m(igrid,nlev)
                                                      ENDDO
c
c.......................................................................
c  calcul des fonctions de stabilite
c.......................................................................
c
      if (l_mix.eq.4) then
                                                      DO igrid=1,ngrid 
         sqz(igrid)=1.e-10
         sq(igrid)=1.e-10
                                                      ENDDO
         do ilev=2,nlev-1
                                                      DO igrid=1,ngrid 
           zq=sqrt(q2(igrid,ilev))
           sqz(igrid)
     .     =sqz(igrid)+zq*zlev(igrid,ilev)
     .     *(zlay(igrid,ilev)-zlay(igrid,ilev-1))
           sq(igrid)=sq(igrid)+zq*(zlay(igrid,ilev)-zlay(igrid,ilev-1))
                                                      ENDDO
         enddo
                                                      DO igrid=1,ngrid 
         long0(igrid)=0.2*sqz(igrid)/sq(igrid)
                                                      ENDDO
      else if (l_mix.eq.3) then
         long0(igrid)=long00
      endif

c (abd 5 2)      print*,'LONG0=',long0

c-----------------------------------------------------------------------
      DO ilev=2,nlev-1
                                                      DO igrid=1,ngrid 
c-----------------------------------------------------------------------
c
        tmp1=kappa*(zlev(igrid,ilev)-zlev(igrid,1))
        if (l_mix.ge.10) then
            long(igrid,ilev)=l_mix
        else
           long(igrid,ilev)=tmp1/(1.E+0 + tmp1/long0(igrid))
        endif
        long(igrid,ilev)=max(min(long(igrid,ilev)
     s    ,0.5*sqrt(q2(igrid,ilev))/sqrt(max(n2(igrid,ilev),1.e-10)))
     s    ,5.)

        gn=-long(igrid,ilev)**2 / q2(igrid,ilev)
     &                                           * n2(igrid,ilev)
        gm=long(igrid,ilev)**2 / q2(igrid,ilev)
     &                                           * m2(igrid,ilev)
c
        gninf=.false.
        gnsup=.false.
        long(igrid,ilev)=long(igrid,ilev)
        long(igrid,ilev)=long(igrid,ilev)
c
        IF (gn.lt.gnmin) THEN
          gninf=.true.
          gn=gnmin
        ENDIF
c
        IF (gn.gt.gnmax) THEN
          gnsup=.true.
          gn=gnmax
        ENDIF
c
        sn(igrid,ilev)=cn1/(1.E+0 +cn2*gn)
        sm(igrid,ilev)=
     &    (cm1+cm2*gn)
     &   /( (1.E+0 +cm3*gn)
     &     *(1.E+0 +cm4*gn) )
c
        IF ((gninf).or.(gnsup)) THEN
          snq2(igrid,ilev)=0.E+0
          smq2(igrid,ilev)=0.E+0
        ELSE
          snq2(igrid,ilev)=
     &     -gn
     &     *(-cn1*cn2/(1.E+0 +cn2*gn)**2 )
          smq2(igrid,ilev)=
     &     -gn
     &     *( cm2*(1.E+0 +cm3*gn)
     &           *(1.E+0 +cm4*gn)
     &       -( cm3*(1.E+0 +cm4*gn)
     &         +cm4*(1.E+0 +cm3*gn) )
     &       *(cm1+cm2*gn)            )
     &     /( (1.E+0 +cm3*gn)
     &       *(1.E+0 +cm4*gn) )**2
        ENDIF
c
c abd
c        if(ilev.le.57.and.ilev.ge.37) then
c            print*,'L=',ilev,'   GN=',gn,'  SM=',sm(igrid,ilev)
c        endif
c --->
c       la decomposition de Taylor en q2 n'a de sens que
c       dans les cas stratifies ou sn et sm sont quasi
c       proportionnels a q2. ailleurs on laisse le meme
c       algorithme car l'ajustement convectif fait le travail.
c       mais c'est delirant quand sn et snq2 n'ont pas le meme
c       signe : dans ces cas, on ne fait pas la decomposition.
c<---
c
        IF (snq2(igrid,ilev)*sn(igrid,ilev).le.0.E+0)
     &      snq2(igrid,ilev)=0.E+0
        IF (smq2(igrid,ilev)*sm(igrid,ilev).le.0.E+0)
     &      smq2(igrid,ilev)=0.E+0
c
C   Correction pour les couches stables.
C   Schema repris de JHoltzlag Boville, lui meme venant de...

        if (1.eq.1) then
        snstable=1.-zlev(igrid,ilev)
     s     /(700.*max(ustar(igrid),0.0001))
        snstable=1.-zlev(igrid,ilev)/400.
        snstable=max(snstable,0.)
        snstable=snstable*snstable

c abde       print*,'SN ',ilev,sn(1,ilev),snstable
        if (sn(igrid,ilev).lt.snstable) then
           sn(igrid,ilev)=snstable
           snq2(igrid,ilev)=0.
        endif

        if (sm(igrid,ilev).lt.snstable) then
           sm(igrid,ilev)=snstable
           smq2(igrid,ilev)=0.
        endif

        endif

c sn : coefficient de stabilite pour n
c snq2 : premier terme du developement limite de sn en q2
c-----------------------------------------------------------------------
                                                      ENDDO
      ENDDO
c-----------------------------------------------------------------------
c
c.......................................................................
c  calcul de km et kn au debut du pas de temps
c.......................................................................
c
                                                      DO igrid=1,ngrid 
      kn(igrid,1)=knmin
      km(igrid,1)=kmmin
      kmpre(igrid,1)=km(igrid,1)
                                                      ENDDO
c
c-----------------------------------------------------------------------
      DO ilev=2,nlev-1
                                                      DO igrid=1,ngrid 
c-----------------------------------------------------------------------
c
        kn(igrid,ilev)=long(igrid,ilev)*q(igrid,ilev)
     &                                         *sn(igrid,ilev)
        km(igrid,ilev)=long(igrid,ilev)*q(igrid,ilev)
     &                                         *sm(igrid,ilev)
        kmpre(igrid,ilev)=km(igrid,ilev)
c
c-----------------------------------------------------------------------
                                                      ENDDO
      ENDDO
c-----------------------------------------------------------------------
c
                                                      DO igrid=1,ngrid 
      kn(igrid,nlev)=kn(igrid,nlev-1)
      km(igrid,nlev)=km(igrid,nlev-1)
      kmpre(igrid,nlev)=km(igrid,nlev)
                                                      ENDDO
c
c.......................................................................
c  boucle sur les niveaux 2 a nlev-1
c.......................................................................
c
c---->
      DO 10001 ilev=2,nlev-1
c---->
      DO 10002 igrid=1,ngrid 
c
c.......................................................................
c
c  calcul des termes sources et puits de l'equation de q2
c  ------------------------------------------------------
c
        knq3=kn(igrid,ilev)*snq2(igrid,ilev)
     &                                    /sn(igrid,ilev)
        kmq3=km(igrid,ilev)*smq2(igrid,ilev)
     &                                    /sm(igrid,ilev)
c
        termq=0.E+0
        termq3=0.E+0
        termqm2=0.E+0
        termq3m2=0.E+0
c
        tmp1=dt*2.E+0 *km(igrid,ilev)*m2(igrid,ilev)
        tmp2=dt*2.E+0 *kmq3*m2(igrid,ilev)
        termqm2=termqm2
     &    +dt*2.E+0 *km(igrid,ilev)*m2(igrid,ilev)
     &    -dt*2.E+0 *kmq3*m2(igrid,ilev)
        termq3m2=termq3m2
     &    +dt*2.E+0 *kmq3*m2(igrid,ilev)
c 
        termq=termq
     &    -dt*2.E+0 *kn(igrid,ilev)*n2(igrid,ilev)
     &    +dt*2.E+0 *knq3*n2(igrid,ilev)
        termq3=termq3
     &    -dt*2.E+0 *knq3*n2(igrid,ilev)
c
        termq3=termq3
     &    -dt*2.E+0 *q(igrid,ilev)**3 / (b1*long(igrid,ilev))
c
c.......................................................................
c
c  resolution stationnaire couplee avec le gradient de vitesse local
c  -----------------------------------------------------------------
c
c  -----{on cherche le cisaillement qui annule l'equation de q^2
c        supposee en q3}
c
        tmp1=termq+termq3
        tmp2=termqm2+termq3m2
        m2cstat=m2(igrid,ilev)
     &      -(tmp1+tmp2)/(dt*2.E+0*km(igrid,ilev))
        mcstat=sqrt(m2cstat)

c  abde      print*,'M2 L=',ilev,mpre(igrid,ilev),mcstat
c
c  -----{puis on ecrit la valeur de q qui annule l'equation de m
c        supposee en q3}
c
        IF (ilev.eq.2) THEN
          kmcstat=1.E+0 / mcstat
     &    *( unsdz(igrid,ilev)*kmpre(igrid,ilev+1)
     &                        *mpre(igrid,ilev+1)
     &      +unsdz(igrid,ilev-1)
     &              *cd(igrid)
     &              *( sqrt(u(igrid,3)**2+v(igrid,3)**2)
     &                -mcstat/unsdzdec(igrid,ilev)
     &                -mpre(igrid,ilev+1)/unsdzdec(igrid,ilev+1) )**2)
     &      /( unsdz(igrid,ilev)+unsdz(igrid,ilev-1) )
        ELSE
          kmcstat=1.E+0 / mcstat
     &    *( unsdz(igrid,ilev)*kmpre(igrid,ilev+1)
     &                        *mpre(igrid,ilev+1)
     &      +unsdz(igrid,ilev-1)*kmpre(igrid,ilev-1)
     &                          *mpre(igrid,ilev-1) )
     &      /( unsdz(igrid,ilev)+unsdz(igrid,ilev-1) )
        ENDIF
        tmp2=kmcstat
     &      /( sm(igrid,ilev)/q2(igrid,ilev) )
     &      /long(igrid,ilev)
        qcstat=tmp2**(1.E+0/3.E+0)
        q2cstat=qcstat**2
c
c.......................................................................
c
c  choix de la solution finale
c  ---------------------------
c
          q(igrid,ilev)=qcstat
          q2(igrid,ilev)=q2cstat
          m(igrid,ilev)=mcstat
c abd       if(ilev.le.57.and.ilev.ge.37) then
c           print*,'L=',ilev,'   M2=',m2(igrid,ilev),m2cstat,
c     s     'N2=',n2(igrid,ilev)
c abd       endif
          m2(igrid,ilev)=m2cstat
c
c --->
c       pour des raisons simples q2 est minore 
c<---
c
        IF (q2(igrid,ilev).lt.q2min) THEN
          q2(igrid,ilev)=q2min
          q(igrid,ilev)=sqrt(q2min)
        ENDIF
c
c.......................................................................
c
c  calcul final de kn et km
c  ------------------------
c
        gn=-long(igrid,ilev)**2 / q2(igrid,ilev)
     &                                           * n2(igrid,ilev)
        IF (gn.lt.gnmin) gn=gnmin
        IF (gn.gt.gnmax) gn=gnmax
        sn(igrid,ilev)=cn1/(1.E+0 +cn2*gn)
        sm(igrid,ilev)=
     &    (cm1+cm2*gn)
     &   /( (1.E+0 +cm3*gn)*(1.E+0 +cm4*gn) )
        kn(igrid,ilev)=long(igrid,ilev)*q(igrid,ilev)
     &                 *sn(igrid,ilev)
        km(igrid,ilev)=long(igrid,ilev)*q(igrid,ilev)
     &                 *sm(igrid,ilev)
c abd
c        if(ilev.le.57.and.ilev.ge.37) then
c            print*,'L=',ilev,'   GN=',gn,'  SM=',sm(igrid,ilev)
c        endif
c
c.......................................................................
c
10002 CONTINUE
c
10001 CONTINUE
c
c.......................................................................
c
c
                                                      DO igrid=1,ngrid 
      kn(igrid,1)=knmin
      km(igrid,1)=kmmin
c     kn(igrid,1)=cd(igrid)
c     km(igrid,1)=cd(igrid)
      q2(igrid,nlev)=q2(igrid,nlev-1)
      q(igrid,nlev)=q(igrid,nlev-1)
      kn(igrid,nlev)=kn(igrid,nlev-1)
      km(igrid,nlev)=km(igrid,nlev-1)
                                                      ENDDO
c
c  CALCUL DE LA DIFFUSION VERTICALE DE Q2
      if (1.eq.1) then

        do ilev=2,klev-1
           sss=sss+plev(1,ilev-1)-plev(1,ilev+1)
           sssq=sssq+(plev(1,ilev-1)-plev(1,ilev+1))*q2(1,ilev)
        enddo
c        print*,'Q2moy avant',sssq/sss
c       print*,'Q2q20 ',(q2(1,ilev),ilev=1,10)
c       print*,'Q2km0 ',(km(1,ilev),ilev=1,10)
c   ! C'est quoi ca qu'etait dans l'original???
c       do igrid=1,ngrid
c          q2(igrid,1)=10.
c       enddo
c        q2s=q2
c       do iii=1,10
c       call vdif_q2(dt,g,rconst,plev,temp,km,q2)
c       do ilev=1,klev+1
c          write(iii+49,*) q2(1,ilev),zlev(1,ilev)
c       enddo
c       enddo
c       stop
c       do ilev=1,klev
c          print*,zlev(1,ilev),q2s(1,ilev),q2(1,ilev)
c       enddo
c        q2s=q2-q2s
c       do ilev=1,klev
c          print*,q2s(1,ilev),zlev(1,ilev)
c       enddo
        do ilev=2,klev-1
           sss=sss+plev(1,ilev-1)-plev(1,ilev+1)
           sssq=sssq+(plev(1,ilev-1)-plev(1,ilev+1))*q2(1,ilev)
        enddo
        print*,'Q2moy apres',sssq/sss
c
c
        do ilev=1,nlev
           do igrid=1,ngrid
              q2(igrid,ilev)=max(q2(igrid,ilev),q2min)
              q(igrid,ilev)=sqrt(q2(igrid,ilev))

c.......................................................................
c
c  calcul final de kn et km
c  ------------------------
c
        gn=-long(igrid,ilev)**2 / q2(igrid,ilev)
     &                                           * n2(igrid,ilev)
        IF (gn.lt.gnmin) gn=gnmin
        IF (gn.gt.gnmax) gn=gnmax
        sn(igrid,ilev)=cn1/(1.E+0 +cn2*gn)
        sm(igrid,ilev)=
     &    (cm1+cm2*gn)
     &   /( (1.E+0 +cm3*gn)*(1.E+0 +cm4*gn) )
C   Correction pour les couches stables.
C   Schema repris de JHoltzlag Boville, lui meme venant de...

        if (1.eq.1) then
        snstable=1.-zlev(igrid,ilev)
     s     /(700.*max(ustar(igrid),0.0001))
        snstable=1.-zlev(igrid,ilev)/400.
        snstable=max(snstable,0.)
        snstable=snstable*snstable

c abde      print*,'SN ',ilev,sn(1,ilev),snstable
        if (sn(igrid,ilev).lt.snstable) then
           sn(igrid,ilev)=snstable
           snq2(igrid,ilev)=0.
        endif

        if (sm(igrid,ilev).lt.snstable) then
           sm(igrid,ilev)=snstable
           smq2(igrid,ilev)=0.
        endif

        endif

c sn : coefficient de stabilite pour n
        kn(igrid,ilev)=long(igrid,ilev)*q(igrid,ilev)
     &                 *sn(igrid,ilev)
        km(igrid,ilev)=long(igrid,ilev)*q(igrid,ilev)
c
           enddo
        enddo
c       print*,'Q2km1 ',(km(1,ilev),ilev=1,10)

      endif

      RETURN
      END
