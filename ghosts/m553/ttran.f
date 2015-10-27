      subroutine ttran(xj,c,temp,nbf)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
cc
cmp   extended dummy xj,c,temp
cc
      dimension xj(nbf,2),c(nbf,2),temp(nbf,2)
      common / number / zero,pt5,one,two,four,eight
c
      do 30 i=1,nbf
         do 20 j=1,nbf
            xx=zero
            do 10 k=1,nbf
               xx=xx+xj(i,k)*c(k,j)
 10         continue
            temp(i,j)=xx
 20      continue
 30   continue
c
      do 60 i=1,nbf
         do 50 j=1,nbf
            xx=zero
            do 40 k=1,nbf
               xx=xx+temp(k,i)*c(k,j)
 40         continue
c old xj(i,j)=xx
            xj(j,i)=xx
 50      continue
 60   continue
c
      return
      end
