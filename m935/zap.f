*deck @(#)zap.f	5.1  11/6/94
      subroutine zap(xvec,mdim,twalks,nrhs)
      implicit real*8 (a-h,o-z)
      real*8 xvec(mdim,nrhs)
      integer twalks
c
c
      npdim=mdim-twalks
      do 1 i=1,nrhs
         do 2 j=1,twalks
            xvec(npdim+j,i)=0.d0
  2      continue
  1   continue
c
c
      return
      end
