*deck %W%  %G%
      subroutine ttran(xj,c,temp,nbf)
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
      implicit real*8(a-h,o-z)
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
