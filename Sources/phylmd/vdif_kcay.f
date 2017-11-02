module vdif_kcay_m

  IMPLICIT NONE

contains

  SUBROUTINE vdif_kcay(knon, dt, g, zlev, zlay, u, v, teta, cd, q2, q2diag, &
       km, kn, ustar, l_mix)

    ! From LMDZ4/libf/phylmd/vdif_kcay.F, version 1.1, 2004/06/22 11:45:36

    USE dimphy, ONLY: klev, klon
    use yamada_m, only: yamada

    INTEGER knon
    ! knon : nombre de points de grille 

    REAL, intent(in):: dt
    ! dt : pas de temps
    real, intent(in):: g
    ! g : g
    REAL zlev(klon, klev+1)
    ! zlev : altitude a chaque niveau (interface inferieure de la couche
    ! de meme indice)
    REAL zlay(klon, klev)
    ! zlay : altitude au centre de chaque couche
    REAL, intent(in):: u(klon, klev)
    REAL, intent(in):: v(klon, klev)
    ! u, v : vitesse au centre de chaque couche
    ! (en entree : la valeur au debut du pas de temps)
    REAL teta(klon, klev)
    ! teta : temperature potentielle au centre de chaque couche
    ! (en entree : la valeur au debut du pas de temps)
    REAL, intent(in):: cd (:) ! (knon) cdrag, valeur au debut du pas de temps
    REAL q2(klon, klev+1)
    ! q2 : $q^2$ au bas de chaque couche
    ! (en entree : la valeur au debut du pas de temps)
    ! (en sortie : la valeur a la fin du pas de temps)
    REAL q2diag(klon, klev+1)
    REAL km(klon, klev+1)
    ! km : diffusivite turbulente de quantite de mouvement (au bas de chaque
    ! couche)
    ! (en sortie : la valeur a la fin du pas de temps)
    REAL kn(klon, klev+1)
    ! kn : diffusivite turbulente des scalaires (au bas de chaque couche)
    ! (en sortie : la valeur a la fin du pas de temps)
    real, intent(in):: ustar(:) ! (knon)
    integer l_mix

    ! Local:

    real snstable
    real sq(klon), sqz(klon), zq, long0(klon)

    INTEGER nlay, nlev
    ! nlay : nombre de couches 
    ! nlev : nombre de niveaux
    REAL unsdz(klon, klev)
    ! unsdz : 1 sur l'epaisseur de couche
    REAL unsdzdec(klon, klev+1)
    ! unsdzdec : 1 sur la distance entre le centre de la couche et le
    ! centre de la couche inferieure
    REAL q(klon, klev+1)
    ! q : echelle de vitesse au bas de chaque couche
    ! (valeur a la fin du pas de temps)

    REAL kmpre(klon, klev+1)
    ! kmpre : km au debut du pas de temps
    REAL qcstat
    ! qcstat : q : solution stationnaire du probleme couple
    ! (valeur a la fin du pas de temps)
    REAL q2cstat
    ! q2cstat : q2 : solution stationnaire du probleme couple
    ! (valeur a la fin du pas de temps)

    REAL long(klon, klev+1)
    ! long : longueur de melange calculee selon Blackadar

    ! kmq3 : terme en q^3 dans le developpement de km
    ! (valeur au debut du pas de temps)
    ! kmcstat : valeur de km solution stationnaire du systeme {q2 ; du/dz}
    ! (valeur a la fin du pas de temps)
    ! knq3 : terme en q^3 dans le developpement de kn
    ! mcstat : valeur de m solution stationnaire du systeme {q2 ; du/dz}
    ! (valeur a la fin du pas de temps)
    ! m2cstat : valeur de m2 solution stationnaire du systeme {q2 ; du/dz}
    ! (valeur a la fin du pas de temps)
    ! m : valeur a la fin du pas de temps
    ! mpre : valeur au debut du pas de temps
    ! m2 : valeur a la fin du pas de temps
    ! n2 : valeur a la fin du pas de temps

    REAL kmq3
    REAL kmcstat
    REAL knq3
    REAL mcstat
    REAL m2cstat
    REAL m(klon, klev+1)
    REAL mpre(klon, klev+1)
    REAL m2(klon, klev+1)
    REAL n2(klon, klev+1)

    ! gn : intermediaire pour les coefficients de stabilite
    ! gnmin : borne inferieure de gn (-0.23 ou -0.28)
    ! gnmax : borne superieure de gn (0.0233)
    ! gninf : vrai si gn est en dessous de sa borne inferieure
    ! gnsup : vrai si gn est en dessus de sa borne superieure
    ! ri : nombre de Richardson
    ! sn : coefficient de stabilite pour n
    ! snq2 : premier terme du developement limite de sn en q2
    ! sm : coefficient de stabilite pour m
    ! smq2 : premier terme du developement limite de sm en q2

    REAL gn
    REAL gnmin
    REAL gnmax
    LOGICAL gninf
    LOGICAL gnsup
    REAL sn(klon, klev+1)
    REAL snq2(klon, klev+1)
    REAL sm(klon, klev+1)
    REAL smq2(klon, klev+1)

    ! kappa : consatnte de Von Karman (0.4)
    ! long00 : longueur de reference pour le calcul de long (160)
    ! a1, a2, b1, b2, c1 : constantes d'origine pour les coefficients
    ! de stabilite (0.92/0.74/16.6/10.1/0.08)
    ! cn1, cn2 : constantes pour sn
    ! cm1, cm2, cm3, cm4 : constantes pour sm

    REAL kappa
    REAL long00
    REAL a1, a2, b1, b2, c1
    REAL cn1, cn2
    REAL cm1, cm2, cm3, cm4

    ! termq : termes en $q$ dans l'equation de q2
    ! termq3 : termes en $q^3$ dans l'equation de q2
    ! termqm2 : termes en $q*m^2$ dans l'equation de q2
    ! termq3m2 : termes en $q^3*m^2$ dans l'equation de q2

    REAL termq
    REAL termq3
    REAL termqm2
    REAL termq3m2

    ! q2min : borne inferieure de q2
    REAL q2min

    ! knmin : borne inferieure de kn
    ! kmmin : borne inferieure de km

    REAL knmin
    REAL kmmin

    INTEGER ilay, ilev, igrid
    REAL tmp1, tmp2

    PARAMETER (kappa=0.4E+0)
    PARAMETER (long00=160.E+0)
    ! PARAMETER (gnmin=-10.E+0)
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
    PARAMETER (nlay=klev)
    PARAMETER (nlev=klev+1)

    PARAMETER (cn1=a2*(1.E+0 -6.E+0 *a1/b1))
    PARAMETER (cn2=-3.E+0 *a2*(6.E+0 *a1+b2))
    PARAMETER (cm1=a1*(1.E+0 -3.E+0 *c1-6.E+0 *a1/b1))
    PARAMETER (cm2=a1*(-3.E+0 *a2*((b2-3.E+0 *a2)*(1.E+0 -6.E+0 *a1/b1) &
         -3.E+0 *c1*(b2+6.E+0 *a1))))
    PARAMETER (cm3=-3.E+0 *a2*(6.E+0 *a1+b2))
    PARAMETER (cm4=-9.E+0 *a1*a2)

    logical:: first = .true.

    !------------------------------------------------------------

    ! traitment des valeur de q2 en entree

    ! Initialisation de q2

    call yamada(knon, g, zlev, zlay, u, v, teta, q2diag, km, kn)
    if (first.and.1.eq.1) then
       first=.false.
       q2=q2diag
    endif

    DO ilev=1, nlev
       DO igrid=1, knon 
          q2(igrid, ilev)=amax1(q2(igrid, ilev), q2min)
          q(igrid, ilev)=sqrt(q2(igrid, ilev))
       ENDDO
    ENDDO

    DO igrid=1, knon 
       tmp1=cd(igrid)*(u(igrid, 1)**2+v(igrid, 1)**2)
       q2(igrid, 1)=b1**(2.E+0/3.E+0)*tmp1
       q2(igrid, 1)=amax1(q2(igrid, 1), q2min)
       q(igrid, 1)=sqrt(q2(igrid, 1))
    ENDDO

    ! les increments verticaux

    ! allerte !c
    ! zlev n'est pas declare a nlev !c
    DO igrid=1, knon 
       zlev(igrid, nlev)=zlay(igrid, nlay) &
            +( zlay(igrid, nlay) - zlev(igrid, nlev-1) )
    ENDDO
    ! allerte !c

    DO ilay=1, nlay
       DO igrid=1, knon 
          unsdz(igrid, ilay)=1.E+0/(zlev(igrid, ilay+1)-zlev(igrid, ilay))
       ENDDO
    ENDDO
    DO igrid=1, knon 
       unsdzdec(igrid, 1)=1.E+0/(zlay(igrid, 1)-zlev(igrid, 1))
    ENDDO
    DO ilay=2, nlay
       DO igrid=1, knon 
          unsdzdec(igrid, ilay)=1.E+0/(zlay(igrid, ilay)-zlay(igrid, ilay-1))
       ENDDO
    ENDDO
    DO igrid=1, knon 
       unsdzdec(igrid, nlay+1)=1.E+0/(zlev(igrid, nlay+1)-zlay(igrid, nlay))
    ENDDO

    ! le cisaillement et le gradient de temperature

    DO igrid=1, knon 
       m2(igrid, 1)=(unsdzdec(igrid, 1) &
            *u(igrid, 1))**2 &
            +(unsdzdec(igrid, 1) &
            *v(igrid, 1))**2
       m(igrid, 1)=sqrt(m2(igrid, 1))
       mpre(igrid, 1)=m(igrid, 1)
    ENDDO

    DO ilev=2, nlev-1
       DO igrid=1, knon 

          n2(igrid, ilev)=g*unsdzdec(igrid, ilev) &
               *(teta(igrid, ilev)-teta(igrid, ilev-1)) &
               /(teta(igrid, ilev)+teta(igrid, ilev-1)) *2.E+0

          ! on ne sais traiter que les cas stratifies. et l'ajustement
          ! convectif est cense faire en sorte que seul des
          ! configurations stratifiees soient rencontrees en entree de
          ! cette routine. mais, bon ... on sait jamais (meme on sait
          ! que n2 prends quelques valeurs negatives ... parfois)
          ! alors :

          IF (n2(igrid, ilev).lt.0.E+0) THEN
             n2(igrid, ilev)=0.E+0
          ENDIF

          m2(igrid, ilev)=(unsdzdec(igrid, ilev) &
               *(u(igrid, ilev)-u(igrid, ilev-1)))**2 &
               +(unsdzdec(igrid, ilev) &
               *(v(igrid, ilev)-v(igrid, ilev-1)))**2
          m(igrid, ilev)=sqrt(m2(igrid, ilev))
          mpre(igrid, ilev)=m(igrid, ilev)

       ENDDO
    ENDDO

    DO igrid=1, knon 
       m2(igrid, nlev)=m2(igrid, nlev-1)
       m(igrid, nlev)=m(igrid, nlev-1)
       mpre(igrid, nlev)=m(igrid, nlev)
    ENDDO

    ! calcul des fonctions de stabilite

    if (l_mix.eq.4) then
       DO igrid=1, knon 
          sqz(igrid)=1.e-10
          sq(igrid)=1.e-10
       ENDDO
       do ilev=2, nlev-1
          DO igrid=1, knon 
             zq=sqrt(q2(igrid, ilev))
             sqz(igrid) &
                  =sqz(igrid)+zq*zlev(igrid, ilev) &
                  *(zlay(igrid, ilev)-zlay(igrid, ilev-1))
             sq(igrid)=sq(igrid)+zq*(zlay(igrid, ilev)-zlay(igrid, ilev-1))
          ENDDO
       enddo
       DO igrid=1, knon 
          long0(igrid)=0.2*sqz(igrid)/sq(igrid)
       ENDDO
    else if (l_mix.eq.3) then
       long0(igrid)=long00
    endif

    DO ilev=2, nlev-1
       DO igrid=1, knon 
          tmp1=kappa*(zlev(igrid, ilev)-zlev(igrid, 1))
          if (l_mix.ge.10) then
             long(igrid, ilev)=l_mix
          else
             long(igrid, ilev)=tmp1/(1.E+0 + tmp1/long0(igrid))
          endif
          long(igrid, ilev)=max(min(long(igrid, ilev) &
               , 0.5*sqrt(q2(igrid, ilev))/sqrt(max(n2(igrid, ilev), 1.e-10))) &
               , 5.)

          gn=-long(igrid, ilev)**2 / q2(igrid, ilev) &
               * n2(igrid, ilev)
          gninf=.false.
          gnsup=.false.
          long(igrid, ilev)=long(igrid, ilev)
          long(igrid, ilev)=long(igrid, ilev)

          IF (gn.lt.gnmin) THEN
             gninf=.true.
             gn=gnmin
          ENDIF

          IF (gn.gt.gnmax) THEN
             gnsup=.true.
             gn=gnmax
          ENDIF

          sn(igrid, ilev)=cn1/(1.E+0 +cn2*gn)
          sm(igrid, ilev)= &
               (cm1+cm2*gn) &
               /( (1.E+0 +cm3*gn) &
               *(1.E+0 +cm4*gn) )

          IF ((gninf).or.(gnsup)) THEN
             snq2(igrid, ilev)=0.E+0
             smq2(igrid, ilev)=0.E+0
          ELSE
             snq2(igrid, ilev)= &
                  -gn &
                  *(-cn1*cn2/(1.E+0 +cn2*gn)**2 )
             smq2(igrid, ilev)= &
                  -gn &
                  *( cm2*(1.E+0 +cm3*gn) &
                  *(1.E+0 +cm4*gn) &
                  -( cm3*(1.E+0 +cm4*gn) &
                  +cm4*(1.E+0 +cm3*gn) ) &
                  *(cm1+cm2*gn) ) &
                  /( (1.E+0 +cm3*gn) &
                  *(1.E+0 +cm4*gn) )**2
          ENDIF
          ! la decomposition de Taylor en q2 n'a de sens que dans les
          ! cas stratifies ou sn et sm sont quasi proportionnels a
          ! q2. ailleurs on laisse le meme algorithme car l'ajustement
          ! convectif fait le travail.  mais c'est delirant quand sn
          ! et snq2 n'ont pas le meme signe : dans ces cas, on ne fait
          ! pas la decomposition.

          IF (snq2(igrid, ilev)*sn(igrid, ilev).le.0.E+0) &
               snq2(igrid, ilev)=0.E+0
          IF (smq2(igrid, ilev)*sm(igrid, ilev).le.0.E+0) &
               smq2(igrid, ilev)=0.E+0

          ! Correction pour les couches stables.
          ! Schema repris de JHoltzlag Boville, lui meme venant de...

          snstable=1.-zlev(igrid, ilev) &
               /(700.*max(ustar(igrid), 0.0001))
          snstable=1.-zlev(igrid, ilev)/400.
          snstable=max(snstable, 0.)
          snstable=snstable*snstable

          if (sn(igrid, ilev).lt.snstable) then
             sn(igrid, ilev)=snstable
             snq2(igrid, ilev)=0.
          endif

          if (sm(igrid, ilev).lt.snstable) then
             sm(igrid, ilev)=snstable
             smq2(igrid, ilev)=0.
          endif

          ! sn : coefficient de stabilite pour n
          ! snq2 : premier terme du developement limite de sn en q2
       ENDDO
    ENDDO

    ! calcul de km et kn au debut du pas de temps

    DO igrid=1, knon 
       kn(igrid, 1)=knmin
       km(igrid, 1)=kmmin
       kmpre(igrid, 1)=km(igrid, 1)
    ENDDO

    DO ilev=2, nlev-1
       DO igrid=1, knon 
          kn(igrid, ilev)=long(igrid, ilev)*q(igrid, ilev) &
               *sn(igrid, ilev)
          km(igrid, ilev)=long(igrid, ilev)*q(igrid, ilev) &
               *sm(igrid, ilev)
          kmpre(igrid, ilev)=km(igrid, ilev)
       ENDDO
    ENDDO

    DO igrid=1, knon 
       kn(igrid, nlev)=kn(igrid, nlev-1)
       km(igrid, nlev)=km(igrid, nlev-1)
       kmpre(igrid, nlev)=km(igrid, nlev)
    ENDDO

    ! boucle sur les niveaux 2 a nlev-1

    DO ilev=2, nlev-1
       DO igrid=1, knon 
          ! calcul des termes sources et puits de l'equation de q2

          knq3=kn(igrid, ilev)*snq2(igrid, ilev) &
               /sn(igrid, ilev)
          kmq3=km(igrid, ilev)*smq2(igrid, ilev) &
               /sm(igrid, ilev)

          termq=0.E+0
          termq3=0.E+0
          termqm2=0.E+0
          termq3m2=0.E+0

          tmp1=dt*2.E+0 *km(igrid, ilev)*m2(igrid, ilev)
          tmp2=dt*2.E+0 *kmq3*m2(igrid, ilev)
          termqm2=termqm2 &
               +dt*2.E+0 *km(igrid, ilev)*m2(igrid, ilev) &
               -dt*2.E+0 *kmq3*m2(igrid, ilev)
          termq3m2=termq3m2 &
               +dt*2.E+0 *kmq3*m2(igrid, ilev)

          termq=termq &
               -dt*2.E+0 *kn(igrid, ilev)*n2(igrid, ilev) &
               +dt*2.E+0 *knq3*n2(igrid, ilev)
          termq3=termq3 &
               -dt*2.E+0 *knq3*n2(igrid, ilev)

          termq3=termq3 &
               -dt*2.E+0 *q(igrid, ilev)**3 / (b1*long(igrid, ilev))

          ! resolution stationnaire couplee avec le gradient de vitesse local

          ! on cherche le cisaillement qui annule l'equation de q^2
          ! supposee en q3

          tmp1=termq+termq3
          tmp2=termqm2+termq3m2
          m2cstat=m2(igrid, ilev) &
               -(tmp1+tmp2)/(dt*2.E+0*km(igrid, ilev))
          mcstat=sqrt(m2cstat)

          ! puis on ecrit la valeur de q qui annule l'equation de m
          ! supposee en q3

          IF (ilev.eq.2) THEN
             kmcstat=1.E+0 / mcstat &
                  *( unsdz(igrid, ilev)*kmpre(igrid, ilev+1) &
                  *mpre(igrid, ilev+1) &
                  +unsdz(igrid, ilev-1) &
                  *cd(igrid) &
                  *( sqrt(u(igrid, 3)**2+v(igrid, 3)**2) &
                  -mcstat/unsdzdec(igrid, ilev) &
                  -mpre(igrid, ilev+1)/unsdzdec(igrid, ilev+1) )**2) &
                  /( unsdz(igrid, ilev)+unsdz(igrid, ilev-1) )
          ELSE
             kmcstat=1.E+0 / mcstat &
                  *( unsdz(igrid, ilev)*kmpre(igrid, ilev+1) &
                  *mpre(igrid, ilev+1) &
                  +unsdz(igrid, ilev-1)*kmpre(igrid, ilev-1) &
                  *mpre(igrid, ilev-1) ) &
                  /( unsdz(igrid, ilev)+unsdz(igrid, ilev-1) )
          ENDIF
          tmp2=kmcstat &
               /( sm(igrid, ilev)/q2(igrid, ilev) ) &
               /long(igrid, ilev)
          qcstat=tmp2**(1.E+0/3.E+0)
          q2cstat=qcstat**2

          ! choix de la solution finale

          q(igrid, ilev)=qcstat
          q2(igrid, ilev)=q2cstat
          m(igrid, ilev)=mcstat
          m2(igrid, ilev)=m2cstat

          ! pour des raisons simples q2 est minore 

          IF (q2(igrid, ilev).lt.q2min) THEN
             q2(igrid, ilev)=q2min
             q(igrid, ilev)=sqrt(q2min)
          ENDIF

          ! calcul final de kn et km

          gn=-long(igrid, ilev)**2 / q2(igrid, ilev) &
               * n2(igrid, ilev)
          IF (gn.lt.gnmin) gn=gnmin
          IF (gn.gt.gnmax) gn=gnmax
          sn(igrid, ilev)=cn1/(1.E+0 +cn2*gn)
          sm(igrid, ilev)= &
               (cm1+cm2*gn) &
               /( (1.E+0 +cm3*gn)*(1.E+0 +cm4*gn) )
          kn(igrid, ilev)=long(igrid, ilev)*q(igrid, ilev) &
               *sn(igrid, ilev)
          km(igrid, ilev)=long(igrid, ilev)*q(igrid, ilev) &
               *sm(igrid, ilev)
       end DO
    end DO

    DO igrid=1, knon 
       kn(igrid, 1)=knmin
       km(igrid, 1)=kmmin
       q2(igrid, nlev)=q2(igrid, nlev-1)
       q(igrid, nlev)=q(igrid, nlev-1)
       kn(igrid, nlev)=kn(igrid, nlev-1)
       km(igrid, nlev)=km(igrid, nlev-1)
    ENDDO

    ! CALCUL DE LA DIFFUSION VERTICALE DE Q2
    do ilev=1, nlev
       do igrid=1, knon
          q2(igrid, ilev)=max(q2(igrid, ilev), q2min)
          q(igrid, ilev)=sqrt(q2(igrid, ilev))

          ! calcul final de kn et km

          gn=-long(igrid, ilev)**2 / q2(igrid, ilev) &
               * n2(igrid, ilev)
          IF (gn.lt.gnmin) gn=gnmin
          IF (gn.gt.gnmax) gn=gnmax
          sn(igrid, ilev)=cn1/(1.E+0 +cn2*gn)
          sm(igrid, ilev)= &
               (cm1+cm2*gn) &
               /( (1.E+0 +cm3*gn)*(1.E+0 +cm4*gn) )
          ! Correction pour les couches stables.
          ! Schema repris de JHoltzlag Boville, lui meme venant de...

          snstable=1.-zlev(igrid, ilev) &
               /(700.*max(ustar(igrid), 0.0001))
          snstable=1.-zlev(igrid, ilev)/400.
          snstable=max(snstable, 0.)
          snstable=snstable*snstable

          if (sn(igrid, ilev).lt.snstable) then
             sn(igrid, ilev)=snstable
             snq2(igrid, ilev)=0.
          endif

          if (sm(igrid, ilev).lt.snstable) then
             sm(igrid, ilev)=snstable
             smq2(igrid, ilev)=0.
          endif

          ! sn : coefficient de stabilite pour n
          kn(igrid, ilev)=long(igrid, ilev)*q(igrid, ilev) &
               *sn(igrid, ilev)
          km(igrid, ilev)=long(igrid, ilev)*q(igrid, ilev)
       enddo
    enddo

  END SUBROUTINE vdif_kcay

end module vdif_kcay_m
