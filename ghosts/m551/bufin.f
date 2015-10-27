*deck %W%  %G%
      subroutine bufin(x,len,iu)
C
C***Begin prologue     bufin
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
C***End prologue       bufin
C
      implicit real*8(a-h,o-z)
c
      dimension x(len)
c
      call lnkerr(' BUFIN called ')
c
      return
      end
