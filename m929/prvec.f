*deck @(#)prvec.f	1.2  7/30/91
      subroutine prvec(pvec,nroots,nl2,ncsfs)
      real*8 pvec(ncsfs,nl2,nroots)
      character*3 itoc
      character*80 title
      common /io/ inp,iout
c
      do 10 i=1,nroots
         title='          channel = '//itoc(i)
         call prntrm(title,pvec(1,1,i),ncsfs,nl2,ncsfs,nl2,iout)
 10   continue   
      return
      end
