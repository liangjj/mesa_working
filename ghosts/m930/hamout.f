*deck %W%  %G%
      subroutine hamout(h,ntot,nspin,iu,enuc)
      real*8 h(nspin,*)
c
      do 1 i=1,nspin
         h(i,i)=h(i,i)+enuc
  1   continue
c
      call putham(iu,h,ntot,1)
c
      return
      end
