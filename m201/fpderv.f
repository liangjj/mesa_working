*deck @(#)fpderv.f	5.1  11/6/94
      subroutine fpderv(nvar,newcyc,exit,fpcycl,curvar,nstep,mode,dump,
     $                  tstcrv,pool0,pool1,d1var,d2var,yold,d1vold,
     $                  delvar,xi,h,vname,fzero,f,f1)
c***begin prologue     fpderv.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fpderv.f	5.1   11/6/94
c***purpose            
c***description
c     newcyc is the cycle entry flag.  at the beginning of each pass
c     over the variables, it is .true..  exit is the cycle exit flag.
c     it is .false., until all variables have been stepped.
c
c     curvar is the pointer to the particular variable being stepped
c     either up or down.
c
c     nstep is used to determine when we have enough points to compute
c     derivatives.  it is incremented by the sum of mode and one.
c     if mode=0, nstep=1,2,3...2*nvar, where nvar is the number of
c     variables being optimized and both first and second derivatives
c     can be calculated.  if mode=1, nstep=2,4,6... and only first
c     derivatives are computed.  derivatives are only computed
c     when nstep is even.  thus, with mode=1, derivatives are computed
c     every cycle.  the appropriate formulae are used to compensate for the
c     missing values of the function in this case.
c
c***references
c
c***routines called
c
c***end prologue       fpderv.f
      implicit none
c     --- input variables -----
      integer nvar,fpcycl,nstep,curvar,mode
      logical newcyc,exit,dump,tstcrv
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
      real*8 yold(nvar),d1vold(nvar),delvar(nvar)
      real*8 xi(nvar),h(nvar,nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 pool0(nvar),d1var(nvar),d2var(nvar)
      real*8 pool1(nvar)
      real*8 fzero,f,f1(4)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,modulo
      character*16 tname
      real*8 rmsd,sdot,rmslo,rmshi,step2,stepsq,sign,dinc
      real*8 zero,one,two
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
      parameter (rmslo=1.0d-03, rmshi=1.6d-01)
c
      common/io/inp,iout
c
 2001 format(5x,8x,'variable',10x,'de/dx',8x,'de2/dx2')
 2002 format(5x,a16,5x,f10.6,5x,f10.6)
 2003 format(5x,'in deriv, variable ',i3,' incremented:  was ',e20.10,
     $       '  stepped by ',e20.10,'  is now ',e20.10,'.')
 2005 format(5x,'compute first and second derivatives:')
 2006 format(5x,'compute first derivatives.')
c
c
      exit=.false.
      if(newcyc) then
c        --- initial entry for a new pass.  set the pointers.
         nstep=0
         curvar=1
         sign=one
         fzero=f
         newcyc=.false.
c
c        --- calculate second derivatives if the rms displacement is not
c            between rmslo and rmshi.  if it is too low no progress is
c            being made, probably due to bad derivatives near the bottom.  if
c            it is too high, we've gone too far to assume the second
c            derivatives have remained constant.
         call vsub(yold,pool0,pool1,nvar)
         rmsd=sqrt(sdot(nvar,yold,1,yold,1)/float(nvar))
         if(rmsd.lt.rmslo.or.rmsd.gt.rmshi) then
            mode=0
            write(iout,2005)
         else
            mode=1
            write(iout,2006)
         endif
c
c        --- save pool1 in yold and transfer pool0 to pool1.
         call vmove(yold,pool1,nvar)
         call vmove(pool1,pool0,nvar)
c        save d1var in d1vold, if it is available.
         if(fpcycl.gt.0) call vmove(d1vold,d1var,nvar)
      else
c
c        --- a continuation.  can we compute derivatives yet?
         modulo=mod(nstep,2)
         if(modulo.ne.0) then
c           --- can't compute derivatives.  decrement the variable next time.
            f1(1)=f
            sign=-one
         else
c           --- which derivatives should we compute?
            if(mode.eq.0) then
c              --- do both first and second.
               f1(2)=f
               step2=delvar(curvar)+delvar(curvar)
               stepsq=delvar(curvar)*delvar(curvar)
               d1var(curvar)=(f1(1)-f1(2))/step2
               d2var(curvar)=(f1(1)+f1(2)-2*fzero)/stepsq
c              --- if the second derivative is negative, the function is
c                  tending toward a maximum rather than a minimum.
c                  unless the user wants to ignore the test, exit.
               if(d2var(curvar).le.zero) then
                  if(tstcrv) then
                     call fpdump(nvar,pool0,pool1,delvar,yold,d1var,
     $                           d2var,d1vold,xi,h,vname)
                     call lnkerr(' ')
                  endif
               endif
            else
c
c              --- just do first derivative.  correct it by the existing second
c                  derivative.
               f1(1)=f
               d1var(curvar)=(f1(1)-fzero)/delvar(curvar)
     $                      -d2var(curvar)*delvar(curvar)/two
            endif
c
c           --- prepare to do the next variable.
            sign=one
            pool0(curvar)=pool1(curvar)
            curvar=curvar+1
c
c           --- have we made a complete pass on the variables?
c               if so, test for convergence and return.
            if(curvar.gt.nvar) then
               exit=.true.
               write(iout,2001)
               do 10 i=1,nvar
                  call crjust(vname(i),tname)
                  write(iout,2002) tname,d1var(i),d2var(i)
   10          continue
            endif
         endif
      endif
c
c     --- if we must do more, set up the new variables and exit.
      if(.not.exit) then
         dinc=delvar(curvar)*sign
         pool0(curvar)=pool1(curvar)+dinc
         if(dump) 
     $      write(iout,2003) curvar,pool1(curvar),dinc,pool0(curvar)
         nstep=nstep+mode+1
      endif
c
c
      return
      end
