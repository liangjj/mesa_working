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
      subroutine prep(fi,dfi,ddfi,fj,dfj,ddfj,flsti,flstj,dflsti,
     1                dflstj,scr,wt,npts,ni,nj)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fi, dfi, ddfi, fj, dfj, ddfj, flsti, dflsti, flstj, dflstj
      real*8 scr, wt
      dimension fi(npts,ni), dfi(npts,ni), ddfi(npts,ni)
      dimension fj(npts,nj), dfj(npts,nj), ddfj(npts,nj)
      dimension flsti(ni), dflsti(ni), flstj(nj), dflstj(nj)
      dimension scr(*), wt(npts)
      do 10 i=1,ni
         flsti(i)=fi(npts,i)
         dflsti(i)=dfi(npts,i)
   10 continue
      do 20 j=1,nj
         flstj(j)=fj(npts,j)
         dflstj(j)=dfj(npts,j)
   20 continue
      call vsqrt(scr,wt,npts)
      call vmmul(scr,fi,fi,npts,ni)
      call vmmul(scr,fj,fj,npts,nj)
      call vmmul(scr,ddfi,ddfi,npts,ni)    
      call vmmul(scr,ddfj,ddfj,npts,nj)                   
      return
      end
