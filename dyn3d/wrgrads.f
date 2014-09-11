
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/wrgrads.F,v 1.2 2004/06/22 11:45:30
! lmdzadmin Exp $

SUBROUTINE wrgrads(i_f, nl, field, name, titlevar)
  USE gradsdef
  IMPLICIT NONE

  ! Declarations
  ! i_f indice du fichier
  ! nl nombre de couches
  ! field   champ
  ! name    petit nom
  ! titlevar   Titre


  ! arguments
  INTEGER, INTENT(IN):: i_f
  integer nl
  REAL, INTENT(IN):: field(imx*jmx*lmx)
  CHARACTER(len=*) name, titlevar
  CHARACTER(len=10) file

  ! local

  INTEGER im, jm, lm, i, j, l, iv, iii, iji, iif, ijf

  LOGICAL writectl


  writectl = .FALSE.

  PRINT *, i_f, iid(i_f), jid(i_f), ifd(i_f), jfd(i_f)
  iii = iid(i_f)
  iji = jid(i_f)
  iif = ifd(i_f)
  ijf = jfd(i_f)
  im = iif - iii + 1
  jm = ijf - iji + 1
  lm = lmd(i_f)

  PRINT *, 'im,jm,lm,name,firsttime(i_f)'
  PRINT *, im, jm, lm, name, firsttime(i_f)

  IF (firsttime(i_f)) THEN
    IF (name==var(1,i_f)) THEN
      firsttime(i_f) = .FALSE.
      ivar(i_f) = 1
      PRINT *, 'fin de l initialiation de l ecriture du fichier'
      PRINT *, file
      PRINT *, 'fichier no: ', i_f
      PRINT *, 'unit ', unit(i_f)
      PRINT *, 'nvar  ', nvar(i_f)
      PRINT *, 'vars ', (var(iv,i_f), iv=1, nvar(i_f))
    ELSE
      ivar(i_f) = ivar(i_f) + 1
      nvar(i_f) = ivar(i_f)
      var(ivar(i_f), i_f) = name
      tvar(ivar(i_f), i_f) = trim(titlevar)
      nld(ivar(i_f), i_f) = nl
      PRINT *, 'initialisation ecriture de ', var(ivar(i_f), i_f)
      PRINT *, 'i_f ivar(i_f) nld ', i_f, ivar(i_f), nld(ivar(i_f), i_f)
    END IF
    writectl = .TRUE.
    itime(i_f) = 1
  ELSE
    ivar(i_f) = mod(ivar(i_f), nvar(i_f)) + 1
    IF (ivar(i_f)==nvar(i_f)) THEN
      writectl = .TRUE.
      itime(i_f) = itime(i_f) + 1
    END IF

    IF (var(ivar(i_f),i_f)/=name) THEN
      PRINT *, 'Il faut stoker la meme succession de champs a chaque'
      PRINT *, 'pas de temps'
      PRINT *, 'fichier no: ', i_f
      PRINT *, 'unit ', unit(i_f)
      PRINT *, 'nvar  ', nvar(i_f)
      PRINT *, 'vars ', (var(iv,i_f), iv=1, nvar(i_f))

      STOP
    END IF
  END IF

  PRINT *, 'ivar(i_f),nvar(i_f),var(ivar(i_f),i_f),writectl'
  PRINT *, ivar(i_f), nvar(i_f), var(ivar(i_f), i_f), writectl
  DO l = 1, nl
    irec(i_f) = irec(i_f) + 1
    ! print*,'Ecrit rec=',irec(i_f),iii,iif,iji,ijf,
    ! s (l-1)*imd(i_f)*jmd(i_f)+(iji-1)*imd(i_f)+iii
    ! s ,(l-1)*imd(i_f)*jmd(i_f)+(ijf-1)*imd(i_f)+iif
    WRITE (unit(i_f)+1, REC=irec(i_f))((field((l-1)*imd(i_f)*jmd(i_f)+ &
      (j-1)*imd(i_f)+i),i=iii,iif), j=iji, ijf)
  END DO
  IF (writectl) THEN

    file = fichier(i_f)
    ! WARNING! on reecrase le fichier .ctl a chaque ecriture
    OPEN (unit(i_f), FILE=trim(file)//'.ctl', FORM='formatted', &
      STATUS='unknown')
    WRITE (unit(i_f), '(a5,1x,a40)') 'DSET ', '^' // trim(file) // &
      '.dat'

    WRITE (unit(i_f), '(a12)') 'UNDEF 1.0E30'
    WRITE (unit(i_f), '(a5,1x,a40)') 'TITLE ', title(i_f)
    CALL formcoord(unit(i_f), im, xd(iii,i_f), 1., .FALSE., 'XDEF')
    CALL formcoord(unit(i_f), jm, yd(iji,i_f), 1., .TRUE., 'YDEF')
    CALL formcoord(unit(i_f), lm, zd(1,i_f), 1., .FALSE., 'ZDEF')
    WRITE (unit(i_f), '(a4,i10,a30)') 'TDEF ', itime(i_f), &
      ' LINEAR 02JAN1987 1MO '
    WRITE (unit(i_f), '(a4,2x,i5)') 'VARS', nvar(i_f)
    DO iv = 1, nvar(i_f)
      ! print*,'i_f,var(iv,i_f),nld(iv,i_f),nld(iv,i_f)-1/nld(iv,i_f)'
      ! print*,i_f,var(iv,i_f),nld(iv,i_f),nld(iv,i_f)-1/nld(iv,i_f)
      WRITE (unit(i_f), 1000) var(iv, i_f), nld(iv, i_f) - 1/nld(iv, i_f), 99, &
        tvar(iv, i_f)
    END DO
    WRITE (unit(i_f), '(a7)') 'ENDVARS'

1000 FORMAT (A5, 3X, I4, I3, 1X, A39)

    CLOSE (unit(i_f))

  END IF ! writectl

  RETURN

END SUBROUTINE wrgrads

