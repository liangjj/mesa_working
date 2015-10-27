*deck @(#)pvectr.f	5.1  11/6/94
      subroutine pvectr(civec,pvec,roots,nwks,ncsfs,nl2)
      implicit integer(a-z)
      real*8 civec(nwks,roots), pvec(ncsfs,roots*nl2)
      common/io/inp, iout
c
c -- build p-space ci vectors
c --       the p-space vectors reflect guga configuration order and
c --       guga orbital order. the resulting hpp will be in
c --       mo order
c
      ntot=ncsfs*nl2*roots
      call rzero(pvec,ntot)
c
      jx=0
      do 1 i=1,roots
         ix=ncsfs-nwks
         do 2 j=1,nl2
            jx=jx+1
            do 3 k=1,nwks
               pvec(ix+k,jx)=civec(k,i)
  3         continue
         ix=ix-nwks
  2      continue
  1   continue
c
      return
      end
