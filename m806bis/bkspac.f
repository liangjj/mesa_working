*deck @(#)bkspac.f	1.3  8/3/91
      subroutine bkspac
c
      implicit integer (a-z)
c
      character*1 image
c
      common /tapes/ out,input,drttap
      common /buffer/ pos
      common /buffc / image(80)
c
      pos=pos-1
      if (pos.ge.0) return
      backspace input
      backspace input
      call getlin
      pos=79
      return
      end
