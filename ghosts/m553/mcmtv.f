      subroutine mcmtv (hess,c,grad,nbfr,nbfs)
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
c     implicit real*8 (a-h,o-z)
cc
cmp   extended dummy hess,c,grad
cc
      dimension hess(2),c(2),grad(2)
      common / number / zero,pt5,one,two,four,eight
c
      ix=0
      do 20 i=1,nbfs
         xx=zero
         do 10 j=1,nbfr
            ix=ix+1
            xx=xx+hess(ix)*c(j)
 10      continue
         grad(i)=grad(i)+xx
 20   continue
c
      return
      end
