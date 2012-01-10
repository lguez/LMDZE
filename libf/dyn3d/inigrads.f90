module inigrads_m

  implicit none

contains

  subroutine inigrads(i_f, x, fx, xmin, xmax, y, ymin, ymax, fy, z, fz, dt, &
       file, titlel)

    ! From dyn3d/inigrads.F, version 1.1.1.1 2004/05/19 12:53:07

    use gradsdef, only: unit, title, ivar, fichier, firsttime, dtime, iid, &
         ifd, imd, xd, iid, jid, jfd, jmd, yd, lmd, zd, irec

    integer, intent(in):: i_f
    real, intent(in):: x(:), y(:), z(:), fx, fy, fz, dt
    real, intent(in):: xmin, xmax, ymin, ymax
    character(len=*), intent(in):: file, titlel

    ! Variables local to the procedure:
    integer im, jm, lm
    integer i, j, l
    integer:: nf = 0

    !-------------------------------------

    print *, 'Call sequence information: inigrads'

    im = size(x)
    jm = size(y)
    lm = size(z)

    unit(1)=66
    unit(2)=32
    unit(3)=34
    unit(4)=36
    unit(5)=38
    unit(6)=40
    unit(7)=42
    unit(8)=44
    unit(9)=46

    if (i_f.le.nf) stop'verifier les appels a inigrads'

    nf=i_f
    title(i_f)=titlel
    ivar(i_f)=0

    fichier(i_f)=trim(file)

    firsttime(i_f)=.true.
    dtime(i_f)=dt

    iid(i_f)=1
    ifd(i_f)=im
    imd(i_f)=im
    do i=1, im
       xd(i, i_f)=x(i)*fx
       if(xd(i, i_f).lt.xmin) iid(i_f)=i+1
       if(xd(i, i_f).le.xmax) ifd(i_f)=i
    enddo
    print *, 'On stoke du point ', iid(i_f), '  à ', ifd(i_f), ' en x'

    jid(i_f)=1
    jfd(i_f)=jm
    jmd(i_f)=jm
    do j=1, jm
       yd(j, i_f)=y(j)*fy
       if(yd(j, i_f).gt.ymax) jid(i_f)=j+1
       if(yd(j, i_f).ge.ymin) jfd(i_f)=j
    enddo
    print *, 'On stoke du point ', jid(i_f), '  à ', jfd(i_f), ' en y'

    print *, 'fichier(i_f)=', fichier(i_f)
    print *, 4 * (ifd(i_f) - iid(i_f)) * (jfd(i_f) - jid(i_f))
    print *, 'Opening "' // trim(file) // '.dat"...'

    OPEN(unit=unit(i_f)+1, FILE=trim(file)//'.dat', FORM='unformatted', &
         ACCESS='direct', RECL=4*(ifd(i_f)-iid(i_f)+1)*(jfd(i_f)-jid(i_f)+1))

    lmd(i_f)=lm
    do l=1, lm
       zd(l, i_f)=z(l)*fz
    enddo

    irec(i_f)=0
    print *, "i_f = ", i_f
    print *, "imd(i_f) = ", imd(i_f)
    print *, "jmd(i_f) = ", jmd(i_f)
    print *, "lmd(i_f) = ", lmd(i_f)

  end subroutine inigrads

end module inigrads_m
