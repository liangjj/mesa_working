*deck @(#)fixup.f	1.1  11/30/90
      subroutine fixup(npass,dercen,cdint,nints,ndcen)
c
      implicit integer (a-z)
c
      real*8 cdint(nints,3,ndcen)
      integer dercen(4)
c
      if ((npass.eq.2.and.dercen(1).eq.dercen(2)).or.
     #    (npass.eq.3.and.dercen(1).eq.dercen(3))) then
         do 2 coord=1,3
            do 1 i=1,nints
               cdint(i,coord,1)=cdint(i,coord,1)+cdint(i,coord,2)
    1       continue
            call rzero(cdint(1,coord,2),nints)
    2    continue
      end if
c
c
      return
      end
