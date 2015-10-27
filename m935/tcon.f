*deck @(#)tcon.f	5.1  11/6/94
      subroutine tcon(c,iter,nrhs,icnvrg,stopj)
      implicit real*8(a-h,o-z)
      real*8 c(iter,nrhs)
      logical debug
c
      parameter (debug=.false.)
c
      common /io/ inp,iout
c
      if(debug) then
         write(iout,*)' small matrix soln '
         call matout(c,iter,nrhs,iter,nrhs,iout)
      endif
c
      icnvrg=1
      do 80 j=1,nrhs
         xx=sdot(iter,c(1,j),1,c(1,j),1)
         testc=abs(c(iter,j))/sqrt(xx)
         if(testc.gt.stopj) then
            icnvrg=0
            go to 81
         end if
  80  continue
  81  continue
c
      return
      end
