      subroutine rtran(xj,cr,cs,temp,nbfr,nbfs)
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
c     implicit real*8(a-h,o-p,r-z),integer*2(q)
cc
cmp   extended dummy xj,cr,cs,temp
cc
      dimension xj(nbfr,2),cr(nbfr,2),cs(nbfs,2),temp(nbfr,2)
      common / number / zero,pt5,one,two,four,eight
c
      do 30 i=1,nbfr
         do 20 j=1,nbfs
            xx=zero
            do 10 k=1,nbfs
               xx=xx+xj(i,k)*cs(k,j)
 10         continue
            temp(i,j)=xx
 20      continue
 30   continue
c
      do 60 i=1,nbfs
         do 50 j=1,nbfr
            xx=zero
            do 40 k=1,nbfr
               xx=xx+temp(k,i)*cr(k,j)
 40         continue
            xj(j,i)=xx
 50      continue
 60   continue
c
      return
      end
