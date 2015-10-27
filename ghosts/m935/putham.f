*deck %W%  %G%
      subroutine putham(iu,x,len,iw)
      real*8 x(len)
c
      call iosys('write real "sorted hamiltonian" to kohn',
     $           len,x,iw,' ')
c
c
      return
      end
