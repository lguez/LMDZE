!
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/comdissipn.h,v 1.1.1.1 2004/05/19 12:53:06 lmdzadmin Exp $
!
c-----------------------------------------------------------------------
c INCLUDE comdissipn.h

      REAL  tetaudiv, tetaurot, tetah, cdivu, crot, cdivh
c
      COMMON/comdissipn/ tetaudiv(llm),tetaurot(llm),tetah(llm)   ,
     1                        cdivu,      crot,         cdivh

c
c    Les parametres de ce common proviennent des calculs effectues dans 
c             Inidissip  .
c
c-----------------------------------------------------------------------
