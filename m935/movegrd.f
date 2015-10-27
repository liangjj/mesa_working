*deck @(#)movegrd.f	5.1  11/6/94
      subroutine movegrd(gl,ogl,tgl,liter,mvec,iter,nrhs)
      implicit real*8(a-h,o-z)
      real*8 gl(iter,nrhs),ogl(liter,nrhs),tgl(mvec,nrhs)
      logical debug
c
      parameter (debug=.false.)
c
      common /io/ inp,iout
c
      if(debug) then
         write(iout,*)' movegrd ',liter,mvec,iter,nrhs
      endif
c
      if(liter.ne.0) then
         do 1 i=1,nrhs
            call scopy(liter,ogl(1,i),1,gl(1,i),1)
            call scopy(mvec,tgl(1,i),1,gl(liter+1,i),1)
  1      continue
      else
        ntot=nrhs*mvec
        call scopy(ntot,tgl,1,gl,1)
      end if

      return
      end
