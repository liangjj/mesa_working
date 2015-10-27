*deck @(#)svectr.f	5.1  11/6/94
      subroutine svectr(civec,pvec,roots,nwks,ncsfs,nl2,
     #                  nsplit,split)
      implicit integer(a-z)    
      real*8 civec(nwks,roots), pvec(ncsfs,roots*nl2)
      integer split(*),ioff(20)
      integer roots
      common /io/ inp,iout
c
c -- build p-space ci vectors when virtual space has been partitioned (split)
c --       the p-space vectors reflect guga configuration order and
c --       guga orbital order. the resulting hpp will be in
c --       mo order
c
      ix=1
      ic=0
      do 10 i=1,nsplit
         j=nsplit-i+1
         ic=ic+split(i)
         ioff(j)=ix
         ix=ix+split(j)
  10  continue
c
      if(ic.ne.nl2) then
          write(iout,*)' m929: svectr error  ic ne nl2'
          call lnkerr(' m929 svectr ')
      end if
c
      ntot=ncsfs*nl2*roots
      call rzero(pvec,ntot)
c
      jx=0
      do 1 i=1,roots
         do 2 j=1,nsplit
            ix=ncsfs-nwks*ioff(j)
            do 3 k=1,split(j)
               jx=jx+1
               do 4 l=1,nwks
                  pvec(ix+l,jx)=civec(l,i)
  4            continue
               ix=ix-nwks
  3         continue
  2      continue
  1   continue
c
      return
      end
