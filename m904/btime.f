*deck @(#)btime.f	5.1  11/6/94
      subroutine btime(dt,bsec)
      real*8 bsec
      real*8 dum1,dum2
      character*24 dt
c
      call dattim(dt)
      call timing(bsec,dum1,dum2)
c

      return
      end
