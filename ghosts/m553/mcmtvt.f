      subroutine mcmtvt(hess,c,grad,nbfr,nbfs)
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
         xx=c(i)
         do 15 j=1,nbfr
            ix=ix+1
            grad(j)=grad(j)+hess(ix)*xx
 15      continue
 20   continue
c
      return
      end
