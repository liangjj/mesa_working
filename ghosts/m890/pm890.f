*deck @(#)pm890.f	1.1  11/30/90
      subroutine pm890(rcore,icore)
c
c
c***begin prologue     m820
c***date written       871027   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           guga integral ordering
c***author             saxe, paul (lanl)
c***source             @(#)pm890.f	1.1   11/30/90
c
c***purpose            to sort mo ordered mo integrals to guga order.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m820
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 namint
      logical logkey
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
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
      call iosys('open ints as old',0,0,0,namint)
c
      call mn890(rcore,icore,maxcor,'ints','ci','ci')
c
      call iosys('close ints',0,0,0,' ')
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
