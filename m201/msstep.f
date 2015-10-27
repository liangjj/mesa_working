*deck @(#)msstep.f	5.1 11/6/94
      subroutine msstep(nvar,mscycl,xnew,x,xlast,g,glast,s,alpha,
     $                  p,q,z,f,flast,goodstp,steppd,updd2e,srcd2e,
     $                  fncerr,iout)
c***begin prologue     msstep.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)msstep.f	5.1   11/6/94
c***purpose            
c***description
c     this routine computes the step to take for the function evaluation
c
c     input arguments:
c        nvar   ... the number of optimization variables.
c        mscycl ... the cycle counter for the optimization.
c                      if mscycl<2, initalize the optimization.
c        x      ... the current coordinates; real(nvar).
c        xlast  ... the last set of coordinates; real(nvar).
c        g      ... the gradient at x; real(nvar).
c        glast  ... the gradient at xlast; real(nvar).
c        s      ... the current approximation to the inverse hessian;
c                      real(nvar,nvar).
c        alpha  ... the scale factor for the previous step.
c        p      ... scratch; real(nvar).
c        q      ... scratch; real(nvar).
c        z      ... scratch; real(nvar).
c        f      ... the value of the function at x.
c        flast  ... the value of the function at xlast.
c        srcd2e ... source of the second derivatives.
c        fncerr ... estimate of error in value of function.
c        iout   ... unit number of device for posting informative messages.
c     output arguments:
c        xnew   ... the step to take from the current point; real(nvar).
c        s      ... the updated approximation to the inverse hessian;
c                      real(nvar,nvar).
c        alpha  ... the scale factor for the current step.
c        steppd ... a logical variable set to .false. if convergence
c                      appears unlikely.
c        updd2e ... a logical variable set to .false. if inverse
c                      hessian not updated.
c
c     the routine uses several parameters which should be set in
c     the calling routine.
c        delta1 ... 1.0d-08.  do all the others.
c
c***references
c     murtaugh and sargent, computer j., vol. 13, p. 185(1970).
c     vol.13,p.185, 1970.
c
c***routines called
c
c***end prologue       msstep.f
      implicit none
c     --- input variables -----
      integer nvar,mscycl,iout
      character*(*) srcd2e
      logical updd2e,goodstp
      real*8 f,flast,fncerr
c     --- input arrays (unmodified) ---
      real*8 x(nvar),xlast(nvar),g(nvar),glast(nvar)
c     --- input arrays (scratch) ---
      real*8 p(nvar),q(nvar),z(nvar)
c     --- output arrays ---
      real*8 xnew(nvar),s(nvar,nvar)
c     --- output variables ---
      logical steppd
      real*8 alpha
c     --- scratch arrays ---
c     --- local variables ---
      integer iguess,istype
      integer i
      real*8 delta1,delta2,epsiln,stpmin,flowb
      real*8 zero,half,one,two,three
      real*8 gkm1p,pnorm,test,diff,gkp,w,ztz,ztg,zee,c
      real*8 sdot,arrmax
c
      parameter (delta1=1.0d-08,delta2=1.0d-02,epsiln=1.0d-08)
      parameter (stpmin=1.0d-03,iguess=0,istype=1)
      parameter (flowb=0.0d+00)
      parameter (zero=0.0d+00,half=5.0d-01,one=1.0d+00,two=2.0d+00)
      parameter (three=3.0d+00)
c
 1000 format(5x,'msstep: previous step failed, redo from old point ',
     $          'with alpha=',f9.6,'.')
 1010 format(5x,'msstep: c is too small... resetting the inverse ',
     $          'hessian matrix.')
 1020 format(5x,'msstep: inverse hessian set to identity matrix.')
 1030 format(5x,'msstep: alpha(',f9.6,') is too small.',
     $       ' try a different starting point.')
 1040 format(5x,'msstep: hessian indicates a local maximum...',
     $       ' resetting the inverse hessian matrix.')
c
c     --- check the step
      goodstp=.false.
      steppd=.true.
      updd2e=.false.
      if(mscycl.gt.1) then
c        --- see if the last step did any good.
         diff=flast-f
         call embc(p,s,glast,nvar,nvar,1)
         gkm1p=sdot(nvar,glast,1,p,1)
         pnorm=arrmax(p,nvar)
         test=-(epsiln*alpha*gkm1p)
         if(test.le.zero.or.diff+fncerr.le.test) then
            gkp=sdot(nvar,g,1,p,1)
            if(abs(gkp).ge.delta2.and.alpha*pnorm.ge.delta2) then
               zee=three*diff/alpha +gkm1p +gkp
               w=sqrt(zee**2 -gkm1p*gkp)
               alpha=alpha*(one-(gkp+w-zee)/(gkp-gkm1p+two*w))
            else
               alpha=alpha*half
               if(alpha.le.stpmin) then
                  write(iout,1030) alpha
                  steppd=.false.
                  return
               endif
            endif
            call smul(xnew,p,alpha,nvar)
            write(iout,1000) alpha
            return
         endif
c
         goodstp=.true.
         if(srcd2e.ne.'analytic') then
c           --- update the inverse hessian.
            updd2e=.true.
            call vsub(q,g,glast,nvar)
            call ebc(z,s,q,nvar,nvar,1)
            do 20 i=1,nvar
               z(i)=alpha*p(i) -z(i)
   20       continue
            c=sdot(nvar,q,1,z,1)
            ztz=sdot(nvar,z,1,z,1)
            ztg=sdot(nvar,z,1,glast,1)
            if(abs(c).lt.delta2*ztz) then
               write(iout,1010)
               if(istype.eq.0) then
                  write(iout,1020)
                  call rzero(s,nvar*nvar)
                  do 30 i=1,nvar
                     s(i,i)=one
   30             continue
                  call rzero(z,nvar)
                  c=one
               else
                  c=ztz
               endif
            else if(ztg/c.gt.-delta1) then
               write(iout,1040)
               if(istype.eq.0) then
                  write(iout,1020)
                  call rzero(s,nvar*nvar)
                  do 40 i=1,nvar
                     s(i,i)=one
   40             continue
                  call rzero(z,nvar)
                  c=one
               else
                  c=ztz
               endif
            endif
            call vmove(p,z,nvar)
            call smul(p,p,one/c,nvar)
            call apbct(s,p,z,nvar,1,nvar)
         endif
      endif
c
c     --- generate the new direction and scale factor.
      call embc(p,s,g,nvar,nvar,1)
      if(iguess.eq.1) then
         alpha=min(one,two*abs((flowb-f)/gkm1p))
      else
         alpha=one
      endif
      call smul(xnew,p,alpha,nvar)
c
c
      return
      end
