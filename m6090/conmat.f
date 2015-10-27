      subroutine conmat(fns,ddfns,cfn,ddcfn,v,rhs,mtrx,mrhs,scr,
     1                  energy,wt,ipvt,nbfn,npts,prnt)
      implicit integer (a-z)
      real *8 fns, ddfns, v, wt, energy, phase
      complex *16 rhs, mtrx, mrhs, scr, cfn, ddcfn, eye, tandel
      character *80 title
      logical prnt
      dimension fns(npts,nbfn), ddfns(npts,nbfn), mtrx(nbfn+1,nbfn+1)
      dimension cfn(npts,2), ddcfn(npts,2), v(npts), ipvt(nbfn+1)
      dimension rhs(npts), mrhs(nbfn+1), wt(npts), scr(npts,nbfn+1)
      common /io/ inp, iout
      eye=cmplx(0.d0,1.d0)
      do 10 i=1,nbfn
         do 20 j=1,npts
            scr(j,i)=(-.5d0*ddfns(j,i)+(v(j)-energy)*fns(j,i))
     1                            *sqrt(wt(j))
   20    continue
   10 continue
      do 30 j=1,npts
         scr(j,nbfn+1)=(-.5d0*ddcfn(j,2)+(v(j)-energy)*cfn(j,2))
     1                             *sqrt(wt(j))
   30 continue
      if (prnt) then
         title='rectangular matrix'
         call prntcmn(title,scr,npts,nbfn+1,npts,nbfn+1,iout,'e')
      endif 
      call cebtc(mtrx,scr,scr,nbfn+1,npts,nbfn+1)
      do 40 i=1,nbfn+1
         do 50 j=1,npts
            scr(j,i)=scr(j,i)*sqrt(wt(j))
   50    continue
   40 continue
      call cebtc(mrhs,scr,rhs,nbfn+1,npts,1)
      do 60 i=1,nbfn+1
         do 70 j=1,npts
            scr(j,i)=scr(j,i)/wt(j)   
   70    continue
   60 continue   
      call cgefa(mtrx,nbfn+1,nbfn+1,ipvt,info)
      call cgesl(mtrx,nbfn+1,nbfn+1,ipvt,mrhs,0)
      if (prnt) then
          title='solution vector'
          call prntcmn(title,mrhs,nbfn+1,1,nbfn+1,1,iout,'e')
      endif
      tandel=mrhs(nbfn+1)/( 1.d0+ eye*mrhs(nbfn+1) )
      write(iout,80) tandel
c     phase=atan(tandel)
c     write(iout,90) phase
      return
   80 format(//,25x,'tangent of phase shift = ',e15.8,1x,e15.8)
   90 format(//,25x,'phase shift = ',e15.8)
      end






