*deck bkwrd
      subroutine bkwrd(f,wt,intgl,n)
c***begin prologue     bkwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            backward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on (a,b) as (a,r1), (a,r2)
c                      (a,r3)....(a,b) where ri are the quadrature points in
c                      the interval. the integral is initialized as intgl(0)
c                      which is either its last value or zero depending on the
c                      interval.   
c               
c***references
c
c***routines called
c
c***end prologue       bkwrd
c
      implicit integer (a-z)
      dimension f(n), wt(n-1,n), intgl(1:n)
      common /io/ inp, iout
      real*8 f, wt, intgl
      do 10 i=n-1,1,-1
         intgl(i)=intgl(i+1)
         do 20 j=1,n
            intgl(i)=intgl(i)+wt(i,j)*f(j)
   20    continue
   10 continue
      return
      end















