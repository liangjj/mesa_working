*deck @(#)bloksd.f	1.2  7/30/91
      subroutine bloksd
c
      implicit integer (a-z)
c
      common /bloksz/ blksiz,absmax,maxsiz
c
c     data blksiz /30000/         ! default integral block size
c     data absmax /30000/         ! largest possible block size
c     data maxsiz /150000/        ! space for integrals and two vectors
      data blksiz /300000/
      data absmax /1000000/
      data maxsiz /1500000/
      return
      end
