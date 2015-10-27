*deck @(#)fixlag.f	1.1  11/30/90
      subroutine fixlag(xlagt,xlaga,nocnob)
      implicit integer(a-z)
c
      real*8 xlagt(nocnob),xlaga(nocnob)
c
      do 10 i=1,nocnob
         xlagt(i)=xlaga(i)+(2.0d+00)*xlagt(i)
   10 continue
c
      return
      end
