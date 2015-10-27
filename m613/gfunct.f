*deck @(#)gfunct.f	5.1 11/6/94
c***begin prologue     gfunct.f
c***date written       940724   (yymmdd)
c***revision date      11/6/94
c***keywords           gfunct, gaussian integrals
c***author             schneider, barry (lanl)
c***source             @(#)gfunct.f	5.1 11/6/94
c***purpose            basis gaussian integrals over finite intervals
c***description        integral zero to r of x**n exp(-alfa*x**2)
c***                   
c
c***references         bis notes
c
c***routines called
c***end prologue       &M&
      subroutine gfunct(fn,alfa,r,nmax,prnt)
      implicit integer (a-z)
      real *8 fn, alfa, r, gam0, arg, fac
      real *8 ex, fac1, pre
      real *8 sum, series
      logical prnt
      common /io/ inp, iout
      dimension fn(0:nmax)
      call rzero(fn,nmax+1)
c----------------------------------------------------------------------c
c           recursion relationship used for integral is                c
c                                                                      c
c   fn(n) = [ -t**(n-1)*exp(-alpha*t*t) +(n-1)*fn(n-2) ]/(2*alpha)     c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c      define some basic quantities and calculate fn(0) and fn(1)      c
c      if thats all that is required.  f(0) is computed from error     c
c      function called in gam0 and f(1) is analytic.                   c
c----------------------------------------------------------------------c
      arg=r*sqrt(alfa)
      ex=exp(-arg*arg)
      fac=1.d0/(2.d0*alfa)
      fn(0)=gam0(arg)
      if ( nmax.gt.0 ) then
           fn(1)=fac*(1.d+00-ex)
      endif
      if ( nmax.lt.2 ) then
           if (prnt) then
               write(iout,1) (fn(i),i=0,1)
           endif    
           return
      else 
c----------------------------------------------------------------------c
c     test value of arg.  1. if greater than one use forward           c 
c                            recursion.                                c
c                         2. if less than or equal to one use series   c
c                            for two large values of n and recur       c 
c                            downward                                  c
c----------------------------------------------------------------------c
           if (arg.gt.1.d+00) then
c----------------------------------------------------------------------c
c             forward to larger n                                      c
c----------------------------------------------------------------------c
               fac1=-r*fac*ex
               do 10 i=2,nmax
                  fn(i)=fac1+(i-1)*fac*fn(i-2)
                  fac1=fac1*r    
   10          continue
           else
c----------------------------------------------------------------------c
c             backward to smaller n                                    c
c----------------------------------------------------------------------c
c
c----------------------------------------------------------------------c
c            begin by computing last two values using series           c
c            and then going down to smaller values of n                c
c----------------------------------------------------------------------c
               if ( arg.ne.0.d+00 ) then
                    sum=series(arg,nmax)
                    pre=r**(nmax+1)
                    fn(nmax)=pre*sum*ex/(nmax+1)
                    if ( nmax.gt.0 ) then
                         sum=series(arg,nmax-1)
                         pre=pre/r
                         fn(nmax-1)=pre*sum*ex/nmax
                    endif
c        backward recursion
                    if (nmax.gt.1) then
                        pre=pre/r
                        fac1=pre*ex
                        do 20 i=nmax,2,-1
                           fn(i-2)=fac1+2.d+00*alfa*fn(i)
                           fn(i-2)=fn(i-2)/(i-1)
                           fac1=fac1/r
   20                   continue
                    endif
               endif
           endif
      endif
c
c----------------------------------------------------------------------c
c                      we are finished                                 c   
c----------------------------------------------------------------------c
      if (prnt) then
          write(iout,1) (fn(i),i=0,nmax)
      endif
      return
    1 format(/,1x,'fn integrals',(/,5e15.8))        
      end
