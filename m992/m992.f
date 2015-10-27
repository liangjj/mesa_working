*deck @(#)m992.f	5.1  11/6/94
      program m992
      implicit integer(a-z)
      common/io/inp, iout
c
c
      inp=5
      iout=6
      open(inp,file='inp',status='old')
      open(iout,file='out',status='unknown')
      call pm992
c
c
      stop
      end
