module cv_uncompress_m

  implicit none

contains

  SUBROUTINE cv_uncompress(nloc,len,ncum,nd,idcum &
       ,iflag &
       ,precip,cbmf &
       ,ft,fq,fu,fv &
       ,Ma,qcondc             &
       ,iflag1 &
       ,precip1,cbmf1 &
       ,ft1,fq1,fu1,fv1 &
       ,Ma1,qcondc1             &
       )
    use cvparam

    ! inputs:
    integer, intent(in):: len, ncum, nd, nloc
    integer idcum(nloc)
    integer iflag(nloc)
    real precip(nloc), cbmf(nloc)
    real ft(nloc,nd), fq(nloc,nd), fu(nloc,nd), fv(nloc,nd)
    real Ma(nloc,nd)
    real qcondc(nloc,nd) !cld

    ! outputs:
    integer iflag1(len)
    real precip1(len), cbmf1(len)
    real ft1(len,nd), fq1(len,nd), fu1(len,nd), fv1(len,nd)
    real Ma1(len,nd)
    real qcondc1(len,nd) !cld

    ! local variables:
    integer i,k

    do  i=1,ncum
       precip1(idcum(i))=precip(i)
       cbmf1(idcum(i))=cbmf(i)
       iflag1(idcum(i))=iflag(i)
    end do

    do  k=1,nl
       do  i=1,ncum
          ft1(idcum(i),k)=ft(i,k)
          fq1(idcum(i),k)=fq(i,k)
          fu1(idcum(i),k)=fu(i,k)
          fv1(idcum(i),k)=fv(i,k)
          Ma1(idcum(i),k)=Ma(i,k)
          qcondc1(idcum(i),k)=qcondc(i,k)
       end do
    end do

  end SUBROUTINE cv_uncompress

end module cv_uncompress_m
