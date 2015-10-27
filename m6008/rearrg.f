*deck @(#)rearrg.f	1.1 9/8/91
c***begin prologue     rearrg
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rearrg, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
***purpose             re-arrange free-bound integrals
c***description        free-bound integrals re-ordered keeping
c***                   only variational orbitals.
c***                   note that in calling program thpb and thmb
c***                   were paseed as same varible. this necesitates
c***                   the two separate loops below.
c***references         schneider and rescigno, physical review
c***routines called    iosys, util and mdutil
c***end prologue       rearrg
      subroutine rearrg(hpb,hmb,thpb,thmb,ndim,n,m)
      implicit integer (a-z)
      real *8 hmb, thmb
      complex*16 hpb, thpb
      dimension hpb(m,ndim), thpb(m,n), hmb(m,ndim), thmb(m,n)
      do 10 i=1,m
         do 20 j=1,n
            thpb(i,j)=hpb(i,j)
   20    continue
   10 continue
      call cc2opy(thpb,hpb,n*m)
      do 30 i=1,m
         do 40 j=1,n
            thmb(i,j)=hmb(i,j)
   40    continue
   30 continue
      call copy(thmb,hmb,n*m)
      return
      end
