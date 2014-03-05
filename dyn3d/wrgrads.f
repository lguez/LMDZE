
! $Header: /home/cvsroot/LMDZ4/libf/dyn3d/wrgrads.F,v 1.2 2004/06/22 11:45:30
! lmdzadmin Exp $

SUBROUTINE wrgrads(if, nl, field, name, titlevar)
  USE gradsdef
  IMPLICIT NONE

  ! Declarations
  ! if indice du fichier
  ! nl nombre de couches
  ! field   champ
  ! name    petit nom
  ! titlevar   Titre


  ! arguments
  INTEGER if, nl
  REAL, INTENT (IN) :: field(imx*jmx*lmx)
  CHARACTER *10 name, file
  CHARACTER *10 titlevar

  ! local

  INTEGER im, jm, lm, i, j, l, lnblnk, iv, iii, iji, iif, ijf

  LOGICAL writectl


  writectl = .FALSE.

  PRINT *, if, iid(if), jid(if), ifd(if), jfd(if)
  iii = iid(if)
  iji = jid(if)
  iif = ifd(if)
  ijf = jfd(if)
  im = iif - iii + 1
  jm = ijf - iji + 1
  lm = lmd(if)

  PRINT *, 'im,jm,lm,name,firsttime(if)'
  PRINT *, im, jm, lm, name, firsttime(if)

  IF (firsttime(if)) THEN
    IF (name==var(1,if)) THEN
      firsttime(if) = .FALSE.
      ivar(if) = 1
      PRINT *, 'fin de l initialiation de l ecriture du fichier'
      PRINT *, file
      PRINT *, 'fichier no: ', if
      PRINT *, 'unit ', unit(if)
      PRINT *, 'nvar  ', nvar(if)
      PRINT *, 'vars ', (var(iv,if), iv=1, nvar(if))
    ELSE
      ivar(if) = ivar(if) + 1
      nvar(if) = ivar(if)
      var(ivar(if), if) = name
      tvar(ivar(if), if) = titlevar(1:lnblnk(titlevar))
      nld(ivar(if), if) = nl
      PRINT *, 'initialisation ecriture de ', var(ivar(if), if)
      PRINT *, 'if ivar(if) nld ', if, ivar(if), nld(ivar(if), if)
    END IF
    writectl = .TRUE.
    itime(if) = 1
  ELSE
    ivar(if) = mod(ivar(if), nvar(if)) + 1
    IF (ivar(if)==nvar(if)) THEN
      writectl = .TRUE.
      itime(if) = itime(if) + 1
    END IF

    IF (var(ivar(if),if)/=name) THEN
      PRINT *, 'Il faut stoker la meme succession de champs a chaque'
      PRINT *, 'pas de temps'
      PRINT *, 'fichier no: ', if
      PRINT *, 'unit ', unit(if)
      PRINT *, 'nvar  ', nvar(if)
      PRINT *, 'vars ', (var(iv,if), iv=1, nvar(if))

      STOP
    END IF
  END IF

  PRINT *, 'ivar(if),nvar(if),var(ivar(if),if),writectl'
  PRINT *, ivar(if), nvar(if), var(ivar(if), if), writectl
  DO l = 1, nl
    irec(if) = irec(if) + 1
    ! print*,'Ecrit rec=',irec(if),iii,iif,iji,ijf,
    ! s (l-1)*imd(if)*jmd(if)+(iji-1)*imd(if)+iii
    ! s ,(l-1)*imd(if)*jmd(if)+(ijf-1)*imd(if)+iif
    WRITE (unit(if)+1, REC=irec(if))((field((l-1)*imd(if)*jmd(if)+ &
      (j-1)*imd(if)+i),i=iii,iif), j=iji, ijf)
  END DO
  IF (writectl) THEN

    file = fichier(if)
    ! WARNING! on reecrase le fichier .ctl a chaque ecriture
    OPEN (unit(if), FILE=file(1:lnblnk(file))//'.ctl', FORM='formatted', &
      STATUS='unknown')
    WRITE (unit(if), '(a5,1x,a40)') 'DSET ', '^' // file(1:lnblnk(file)) // &
      '.dat'

    WRITE (unit(if), '(a12)') 'UNDEF 1.0E30'
    WRITE (unit(if), '(a5,1x,a40)') 'TITLE ', title(if)
    CALL formcoord(unit(if), im, xd(iii,if), 1., .FALSE., 'XDEF')
    CALL formcoord(unit(if), jm, yd(iji,if), 1., .TRUE., 'YDEF')
    CALL formcoord(unit(if), lm, zd(1,if), 1., .FALSE., 'ZDEF')
    WRITE (unit(if), '(a4,i10,a30)') 'TDEF ', itime(if), &
      ' LINEAR 02JAN1987 1MO '
    WRITE (unit(if), '(a4,2x,i5)') 'VARS', nvar(if)
    DO iv = 1, nvar(if)
      ! print*,'if,var(iv,if),nld(iv,if),nld(iv,if)-1/nld(iv,if)'
      ! print*,if,var(iv,if),nld(iv,if),nld(iv,if)-1/nld(iv,if)
      WRITE (unit(if), 1000) var(iv, if), nld(iv, if) - 1/nld(iv, if), 99, &
        tvar(iv, if)
    END DO
    WRITE (unit(if), '(a7)') 'ENDVARS'

1000 FORMAT (A5, 3X, I4, I3, 1X, A39)

    CLOSE (unit(if))

  END IF ! writectl

  RETURN

END SUBROUTINE wrgrads

