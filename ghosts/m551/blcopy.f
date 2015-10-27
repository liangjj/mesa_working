*deck %W%  %G%
      subroutine blcopy(txk,iread,kmkl,nbf,xk)
C
C***Begin prologue     blcopy
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
C***End prologue       blcopy
C
      implicit real*8(a-h,o-z)
c
      dimension txk(kmkl,nbf,nbf),xk(nbf,nbf)
c
      do 20 i=1,nbf
         do 10 j=1,nbf
            xk(j,i)=txk(iread,i,j)
 10      continue
 20   continue
c
      return
      end
