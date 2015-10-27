*deck bser
c***begin prologue     bser
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate regular bessel functions for small argument
c***description        using series.
c***references       
c
c***routines called
c***end prologue       bser
      subroutine bser(x,j,dj,fact,m,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 x, j, dj, fact
      real*8 xfac, xfac2, pre, prex, term, sum, dsum
      real*8 xfacm, zero, one, half, cnvrg
      logical prnt
      character*80 title
      dimension fact(0:*)
      data zero, one, half, cnvrg / 0.d0, 1.d0, .5d0, 1.d-20 /
      mn=m-n 
      if(mn.lt.0) then
         write(iout,1) 
         return
      endif
      xfac=x*half
      xfacm=one/xfac
      xfac2=xfac*xfac
      pre=one/2**n
      sum=zero
      dsum=sum
      if(mn.ne.0) then
         prex=xfac**mn
      else
         prex=one
      endif
      imul=1
      do 10 i=0,100
	 if(prex.le.cnvrg) then
            go to 20
         endif
         term=imul*prex/(fact(i)*fact(m+i))
         sum=sum+term
         dsum=dsum+(m-n+i+i)*term*xfacm
         imul=-imul
         prex=prex*xfac2
 10   continue
      write(iout,2)
      call lnkerr('quit')
 20   j=pre*sum
      dj=pre*dsum*half
      return
 1    format(/,1x,'power of x must be ge 0')
 2    format(/,1x,'no series convergence')
      end



