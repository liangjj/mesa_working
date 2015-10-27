*deck @(#)bkspac.f	5.1  11/6/94
      subroutine bkspac
c
      implicit integer (a-z)
c
      character*1 image
c
      common /tapes/ out,errout,input,drttap
      common /buffer/ pos
      common /buffc / image(80)
c
      pos=pos-1
      if (pos.ge.0) return
      backspace input
      backspace input
      call getlin
      pos=79
c
c
      return
      end
