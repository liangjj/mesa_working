*deck @(#)pm841.f	1.1  11/30/90
      subroutine pm841(rcore,icore)
c
c***begin prologue     m841
c***date written       871119   (yymmdd)
c***revision date      871228   (yymmdd)
c   28 december 1987   bhl at brl
c   revised to sort ao integrals into group order
c   revised version m841
c
c
c***keywords           gradient density matrix ordering
c                      of ao integrals
c***author             saxe, paul (lanl)
c***source             @(#)pm841.f	1.1   11/30/90
c
c***purpose            to sort the two-electron integrals into
c  the same order as the order of two-particle ao density matrix
c  needed by the gradient codes.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       m841
c
      implicit integer (a-z)
c
      character*4096 ops
      character*128 namint
      logical logkey
      integer icore(1)
      real*8 rcore(1)
c
c
      maxcor=1
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
      call mn841(rcore,icore,maxcor)
c
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
