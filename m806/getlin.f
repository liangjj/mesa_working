*deck @(#)getlin.f	5.1  11/6/94
      subroutine getlin
c
      implicit integer (a-z)
c
      character*1 image
c
      common /tapes/  out,input,drttap
      common /buffer/ pos
      common /buffc/  image(80)
c
      pos=0
      read (input,1,end=2) image
    1 format (80a1)
      return
    2 write (out,3) input
    3 format (//,' encountered an end to tape',i2,' while scanning '
     #,          'for orbital codes')
      call lnkerr('could not find enough orbital codes')
      end
