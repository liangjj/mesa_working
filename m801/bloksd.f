*deck @(#)bloksd.f	5.1  11/6/94
      block data bloksd
c
      implicit integer (a-z)
c
      common /bloksz/ blksiz,absmax,maxsiz
c
c     data blksiz /30000/         ! default integral block size
c     data absmax /30000/         ! largest possible block size
c     data maxsiz /150000/        ! space for integrals and two vectors
c      data blksiz /100000/
c      data absmax /100000/
c      data maxsiz /4000000/
      data blksiz /1000000/
      data absmax /1000000/
      data maxsiz /10000000/
      end
