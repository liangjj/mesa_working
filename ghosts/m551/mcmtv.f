*deck %W%  %G%
      subroutine mcmtv (hess,c,grad,nbfr,nbfs)
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
      implicit real*8 (a-h,o-z)
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
