!
! $Header: /home/cvsroot/LMDZ4/libf/filtrez/acc.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
        subroutine acc(vec,d,im)
        dimension vec(im,im),d(im)
        do 10 j=1,im
        do 9 i=1,im
 9	d(i)=vec(i,j)*vec(i,j)
        sum=ssum(im,d,1)
        sum=sqrt(sum)
        do 10 i=1,im
 10	vec(i,j)=vec(i,j)/sum
        return
        end
