      subroutine getvec(c,eigval,num,guessc,itap10,i10)
c
c***begin prologue     getvec
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
c***end prologue       getvec
c
      implicit integer (a-z)
c
      real*8 c(num,num),eigval(num)
      integer i10(200)
c
c     ----- fetch pointers, etc. -----
c
      call lnkerr(' in getvec ')
c
c
      return
      end
