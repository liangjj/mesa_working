*deck prep
c***begin prologue     prep
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, basis
c***author             schneider, barry (nsf)
c***source             
c***purpose            scale basis functions by integration weights
c***                   and fill derivative arrays
c***references       
c
c***routines called
c***end prologue       prep
      subroutine prep(f,df,ddf,flst,dflst,scr,wt,npts,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 f, df, ddf, flst, dflst, scr, wt
      dimension f(npts,n), df(npts,n), ddf(npts,n), flst(n), dflst(n)
      dimension scr(*), wt(npts)
      do 10 i=1,n
         flst(i)=f(npts,i)
         dflst(i)=df(npts,i)
   10 continue
      call vsqrt(scr,wt,npts)
      call vmmul(scr,f,f,npts,n)
      call vmmul(scr,ddf,ddf,npts,n)    
      return
      end
