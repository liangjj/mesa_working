*deck @(#)getham.f	5.1  11/6/94
      subroutine getham(iu,x,len,iw,filtyp)
      implicit integer(a-z)
      integer x(len)
      character*(*) filtyp
c
      call iosys('read integer "sorted hamiltonian" from '//filtyp,
     $            len,x,iw,' ')
c
c
      return
      end
