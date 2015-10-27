*deck @(#)pm891.f	5.1  11/6/94
      subroutine pm891(z,a)
c
c     write out dummy ci vectors.
      implicit integer (a-z)
c
      character*4096 ops
      character itoc*4
      integer a(1)
      real*8 z(1)
c
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- enable timing routines if requested -----
c
c
c     get the number of roots and the number of walks.
      call iosys('read integer nwks from rwf',1,nwks,0,' ')
      nroots=intkey(ops,'ci=nroots',1,' ')
      do 10 i=1,nroots
         call iosys('create real "ci root '//itoc(i)//'" on rwf',nwks,
     $               0,0,' ')
   10 continue
c
c     ----- stop timing -----
c
c
c     ----- and exit with grace -----
c
      call chainx(0)
c
c
      return
      end
