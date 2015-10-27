*deck %W%  %G%
      subroutine avadd(a,b,n)
C
C***Begin prologue     avadd
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
C***End prologue       avadd
C
      implicit integer(a-z)
c
      real*8 a(n),b(n)
c
      do 1 i=1,n
         a(i)=a(i)+b(i)
 1    continue
c
      return
      end
