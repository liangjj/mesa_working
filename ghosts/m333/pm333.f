*deck @(#)pm333.f	1.1  11/30/90
      subroutine pm333(rcore,icore)
c
      implicit integer (a-z)
c
      character*4096 ops
      integer icore(*)
      real*8 rcore(*)
c
c..bhl..unicos      common // rcore(1)
c
c..bhl..unicos      equivalence (icore,rcore)
c
c     ----- open the read-write file -----
c
c..bhl..unicos      call drum
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- enable timing routines if requested -----
c
c
c     ----- call the routines to make the supermatrices -----
c
      call mn333(rcore,icore,maxcor)
c
c     ----- stop timing -----
c
c
c     ----- and exit with grace -----
c
      call chainx(0)
c
c
      stop
      end
