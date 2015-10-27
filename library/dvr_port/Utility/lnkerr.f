*deck @(#)lnkerr.f	5.1  11/6/94
      subroutine lnkerr(string)
c
      implicit integer (a-z)
      common /io/ inp, iout
c
      character*(*) string
c
      write(iout,1) string
      stop
c
       return
 1     format(a80)
       end
