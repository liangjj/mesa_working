*deck @(#)onel.f	5.1  11/6/94
      subroutine onel (h,tr,nnp)
c
c
      implicit integer(a-z)
c
      real*8 h(nnp), tr(nnp)
c
c     read in the t and v one electron integrals and add together
c     results stored as triangle in h
c
      call iosys ('read real "kinetic integrals" from rwf',
     $             nnp,tr,0,' ')
      call iosys ('read real "potential integrals" from rwf',
     $             nnp,h,0,' ')
      call vadd (h,h,tr,nnp)
c
c
      return
      end
