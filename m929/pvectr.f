*deck @(#)pvectr.f	5.1  11/6/94
      subroutine pvectr(civec,pvec,nroots,nwks,ncsfs,nl2)
      real*8 civec(nwks,*),pvec(ncsfs,*)
      character*80 title
      common/io/inp, iout
c
c -- build p-space ci vectors
c --       the p-space vectors reflect guga configuration order and
c --       guga orbital order. the resulting hpp will be in
c --       mo order
c
      npvec=nl2*nroots
      ntot=ncsfs*npvec
      call rzero(pvec,ntot)
c
      jx=0
      do 1 i=1,nroots
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
