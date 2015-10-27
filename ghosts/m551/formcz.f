*deck %W%  %G%
      subroutine formcz(cz,nbf)
C
C***Begin prologue     formcz
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
C***End prologue       formcz
C
      implicit real*8 (a-h,o-z)
c
      dimension cz(nbf,2)
c
      common / number / zero,pt5,one,two,four,eight
c
      do 30 i=1,nbf
         do 20 j=1,nbf
            cz(i,j)=zero
 20      continue
         cz(i,i)=one
 30   continue
c
      return
      end
