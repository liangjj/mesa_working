*deck %W%  %G%
      subroutine mcmtvt(hess,c,grad,nbfr,nbfs)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
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
