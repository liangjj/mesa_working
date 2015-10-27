*deck @(#)m994.f	5.1  11/6/94
      program m994
      implicit integer(a-z)
      common/io/inp,iout
      inp=5
      out=6
      open(inp,file='inp',status='old')
      open(iout,file='out',status='unknown')
c
c
      call pm994
c
c
      stop
      end
