c----------------------------------------------------------------------c
c                    gaussian integral routines                        c
c----------------------------------------------------------------------c    
*deck bfunct
c***begin prologue     bfunct
c***date written       880528   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           bfunct, gaussian integrals
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose            basis gaussian integrals over finite intervals
c***description        integral zero to box of x**n exp(-alfa*(x-center)*
c***                   (x-center))
c***                   
c
c***references         bis notes
c
c***routines called
c***end prologue       bfunct
      subroutine bfunct(basint,fn,alfa,center,box,nmax)
      implicit integer (a-z)
      real *8 basint,fn,alfa,center,box,gam0,arg
      real *8 ex, fac1, fac2, pre
      real *8 sum, series
      common /io/ inp, iout
      dimension basint(nmax+1), fn(nmax+1,2), arg(2)
      maxv=nmax+1
      call rzero(basint,maxv)
      call rzero(fn,maxv*2)
c
c----------------------------------------------------------------------c
c           calculate integral from 0 to box of function               c
c           (x)**n *exp(-alfa*(x-center)**2                            c
c                        using                                         c
c           (x-center)**n* exp(alfa*(x-center)**2)                     c
c                        integrals                                     c       
c----------------------------------------------------------------------c
c
c----------------------------------------------------------------------c
c           recursion relationship used for integral is                c
c                                                                      c
c   fn(n) = [ -t**(n-1)*exp(-t*t) +(n-1)*fn(n-2) ]/2                   c
c           these are then summed to make final integrals              c
c----------------------------------------------------------------------c
      arg(2)=sqrt(alfa)
      arg(1)=arg(2)*(box-center)
      arg(2)=arg(2)*center      
c----------------------------------------------------------------------c
c           depending on the value of arg we use forward               c
c           recursion or series plus backward recursion                c
c----------------------------------------------------------------------c
      do 100 iarg=1,2
      ex=exp(-arg(iarg)*arg(iarg))
      if (arg(iarg).gt.1.d+00) then
c----------------------------------------------------------------------c
c           forward recursion on arg                                   c
c----------------------------------------------------------------------c
          fn(1,iarg)=gam0(arg(iarg))
          if (maxv.gt.1) then
              fn(2,iarg)=(1.d+00-ex)*.5d+00
          endif
          if (maxv.gt.2) then
              fac1=-arg(iarg)*ex*.5d+00
              do 10 i=3,maxv
              fn(i,iarg)=fac1+(i-2)*fn(i-2,iarg)*.5d+00
              fac1=fac1*arg(iarg)    
   10         continue
          endif
      else
c----------------------------------------------------------------------c
c            backward recursion on arg                                c
c----------------------------------------------------------------------c
c
c     calculate last two values using series
c
          if (arg(iarg).ne.0.d+00) then
              sum=series(arg(iarg),nmax)
              pre=arg(iarg)**maxv
              fn(maxv,iarg)=pre*sum*ex/maxv
              if (maxv.gt.1) then
                  sum=series(arg(iarg),nmax-1)
                  pre=pre/arg(iarg)
                  fn(maxv-1,iarg)=pre*sum*ex/(maxv-1)
              endif
c     recur
              if (maxv.gt.2) then
                  pre=pre/arg(iarg)
                  fac1=pre*ex
                  do 20 i=maxv,3,-1
                  fn(i-2,iarg)=fac1+2.d+00*fn(i,iarg)
                  fn(i-2,iarg)=fn(i-2,iarg)/(i-2)
                  fac1=fac1/arg(iarg)
   20             continue
              endif
          endif
      endif
  100 continue
c
c
c----------------------------------------------------------------------c
c                  make final integrals                                c
c----------------------------------------------------------------------c
      imul=1
      fac1=1.d+00/sqrt(alfa)
      fac2=fac1
      do 200 i=1,maxv
      fn(i,1)=fn(i,1)+imul*fn(i,2)
      fn(i,1)=fn(i,1)*fac1
      fac1=fac1*fac2
  200 imul=-imul
      if (center.ne.0.d+00) then
          basint(1)=fn(1,1)
          if (maxv.eq.1) return
          basint(2)=center*fn(1,1)+fn(2,1)
          if(maxv.eq.2) return
          fac1=center*center
          do 30 i=3,maxv
          fac2=fac1
          basint(i)=0.d+00
             do 40 j=1,i
             basint(i)=basint(i)+fac2*fn(j,1)
             fac2=fac2*(i-j)/(center*j)
   40        continue
             fac1=fac1*center
   30     continue
      else
      do 50 i=1,maxv
   50 basint(i)=fn(i,1)
      endif
c
c----------------------------------------------------------------------c
c                      we are finished                                 c   
c----------------------------------------------------------------------c
      return
      end
