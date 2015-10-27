*deck %W%  %G%
      subroutine getvec(c,eigval,num,guessc,itap10,i10)
C
C***Begin prologue     getvec
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
C***End prologue       getvec
C
      implicit integer (a-z)
c
      real*8 c(num,num),eigval(num)
      integer i10(200)
c
c     ----- fetch pointers, etc. -----
c
      call lnkerr(' in GETVEC ')
c
c
      return
      end
