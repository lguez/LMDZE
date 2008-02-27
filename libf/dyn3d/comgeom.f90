module comgeom

  use dimens_m, only: iim, jjm
  use paramet_m, only: ip1jmp1, ip1jm

  implicit none

  private iim, jjm, ip1jmp1, ip1jm

  real cu_2d(iim + 1, jjm + 1), cv_2d(iim + 1, jjm)
  real cu(ip1jmp1), cv(ip1jm)
  equivalence (cu, cu_2d), (cv, cv_2d)

  real unscu2_2d(iim + 1, jjm + 1)
  real unscu2(ip1jmp1)
  equivalence (unscu2, unscu2_2d)

  real unscv2_2d(iim + 1,jjm)
  real unscv2(ip1jm)
  equivalence (unscv2, unscv2_2d)

  real aire_2d(iim + 1,jjm + 1), airesurg_2d(iim + 1,jjm + 1)
  real aire(ip1jmp1), airesurg(ip1jmp1)
  equivalence (aire, aire_2d), (airesurg, airesurg_2d)

  real aireu_2d(iim + 1,jjm + 1)           
  real aireu(ip1jmp1)
  equivalence (aireu, aireu_2d)

  real airev_2d(iim + 1,jjm),unsaire_2d(iim + 1,jjm + 1)
  real airev(ip1jm),unsaire(ip1jmp1)
  equivalence (airev, airev_2d), (unsaire, unsaire_2d)

  real apoln,apols

  real unsairez_2d(iim + 1,jjm),airuscv2_2d(iim + 1,jjm)
  real unsairez(ip1jm),airuscv2(ip1jm)
  equivalence (unsairez, unsairez_2d), (airuscv2, airuscv2_2d)

  real airvscu2_2d(iim + 1,jjm)       
  real airvscu2(ip1jm)
  equivalence (airvscu2, airvscu2_2d)

  real aireij1_2d(iim + 1,jjm + 1),aireij2_2d(iim + 1,jjm + 1)
  real aireij1(ip1jmp1),aireij2(ip1jmp1)
  equivalence (aireij1, aireij1_2d), (aireij2, aireij2_2d)

  real aireij3(ip1jmp1)
  real aireij3_2d(iim + 1,jjm + 1)       
  equivalence (aireij3, aireij3_2d)

  real aireij4_2d(iim + 1,jjm + 1), alpha1_2d(iim + 1,jjm + 1)
  real aireij4(ip1jmp1), alpha1(ip1jmp1)
  equivalence (aireij4, aireij4_2d), (alpha1, alpha1_2d)

  real alpha2_2d(iim + 1,jjm + 1)         
  real alpha2(ip1jmp1)
  equivalence (alpha2, alpha2_2d)

  real alpha3_2d(iim + 1,jjm + 1), alpha4_2d(iim + 1,jjm + 1)
  real alpha3(ip1jmp1), alpha4(ip1jmp1)
  equivalence (alpha3, alpha3_2d), (alpha4, alpha4_2d)

  real alpha1p2_2d(iim + 1,jjm + 1)        
  real alpha1p2(ip1jmp1)
  equivalence (alpha1p2, alpha1p2_2d)

  real alpha1p4_2d(iim + 1,jjm + 1),alpha2p3_2d(iim + 1,jjm + 1)
  real alpha1p4(ip1jmp1),alpha2p3(ip1jmp1)
  equivalence (alpha1p4, alpha1p4_2d), (alpha2p3, alpha2p3_2d)

  real alpha3p4(ip1jmp1)
  real alpha3p4_2d(iim + 1,jjm + 1)    
  equivalence (alpha3p4, alpha3p4_2d)

  real fext_2d(iim + 1,jjm),constang_2d(iim + 1,jjm + 1)
  real fext(ip1jm),constang(ip1jmp1)
  equivalence (fext, fext_2d), (constang, constang_2d)

  real rlatu(jjm + 1)
  ! (latitudes of points of the "scalar" and "u" grid, in rad)

  real rlatv(jjm) 
  ! (latitudes of points of the "v" grid, in rad, in decreasing order)

  real rlonu(iim + 1) ! longitudes of points of the "u" grid, in rad

  real rlonv(iim + 1)
  ! (longitudes of points of the "scalar" and "v" grid, in rad)

  real cuvsurcv_2d(iim + 1,jjm),cvsurcuv_2d(iim + 1,jjm)  
  real cuvsurcv(ip1jm),cvsurcuv(ip1jm)
  equivalence (cuvsurcv, cuvsurcv_2d), (cvsurcuv, cvsurcuv_2d)

  real cvusurcu_2d(iim + 1,jjm + 1),cusurcvu_2d(iim + 1,jjm + 1)
  real cvusurcu(ip1jmp1),cusurcvu(ip1jmp1)
  equivalence (cvusurcu, cvusurcu_2d), (cusurcvu, cusurcvu_2d)

  real cuvscvgam1_2d(iim + 1,jjm)
  real cuvscvgam1(ip1jm)
  equivalence (cuvscvgam1, cuvscvgam1_2d)

  real cuvscvgam2_2d(iim + 1,jjm),cvuscugam1_2d(iim + 1,jjm + 1)
  real cuvscvgam2(ip1jm),cvuscugam1(ip1jmp1)
  equivalence (cuvscvgam2, cuvscvgam2_2d), (cvuscugam1, cvuscugam1_2d)

  real cvuscugam2_2d(iim + 1,jjm + 1),cvscuvgam_2d(iim + 1,jjm)
  real cvuscugam2(ip1jmp1),cvscuvgam(ip1jm)
  equivalence (cvuscugam2, cvuscugam2_2d), (cvscuvgam, cvscuvgam_2d)

  real cuscvugam(ip1jmp1)
  real cuscvugam_2d(iim + 1,jjm + 1) 
  equivalence (cuscvugam, cuscvugam_2d)

  real unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2                

  real unsair_gam1_2d(iim + 1,jjm + 1),unsair_gam2_2d(iim + 1,jjm + 1)
  real unsair_gam1(ip1jmp1),unsair_gam2(ip1jmp1)
  equivalence (unsair_gam1, unsair_gam1_2d), (unsair_gam2, unsair_gam2_2d)

  real unsairz_gam_2d(iim + 1,jjm)
  real unsairz_gam(ip1jm)
  equivalence (unsairz_gam, unsairz_gam_2d)

  real aivscu2gam_2d(iim + 1,jjm),aiuscv2gam_2d(iim + 1,jjm)
  real aivscu2gam(ip1jm),aiuscv2gam(ip1jm)
  equivalence (aivscu2gam, aivscu2gam_2d), (aiuscv2gam, aiuscv2gam_2d)

  real xprimu(iim + 1),xprimv(iim + 1)

  save

end module comgeom
