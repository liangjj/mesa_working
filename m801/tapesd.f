*deck @(#)tapesd.f	5.1  11/6/94
      block data tapesd
c
      implicit integer (a-z)
c
      common /tapes/  out,errout,input,drttap
c
      data errout /6/
      data out    /6/
      data input  /5/
      data drttap /8/
      end
