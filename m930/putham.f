*deck @(#)putham.f	5.1  11/6/94
      subroutine putham(iu,x,len,iw,filtyp)
      implicit integer(a-z)
      integer x(len)
      character*(*) filtyp
c
      call iosys('write integer "sorted hamiltonian" to '//filtyp,
     $           len,x,iw,' ')
c
c
      return
      end
