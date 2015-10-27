*deck forwrd
      subroutine forwrd(flm,psilm,j,y,wt,lndex,nl,ltop,n)
c***begin prologue     forwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            forward integration of indefinite integral
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
c***end prologue       forwrd
c
      implicit integer (a-z)
      dimension flm(n), wt(n-1,n), intgl(0:n-1)
      common /io/ inp, iout
      real*8 f, wt, intgl
      do 10 i=1,n-1
         intgl(i)=intgl(i-1)
         do 20 j=1,n
            intgl(i)=intgl(i)+wt(i,j)*f(j)
   20    continue
   10 continue
      return
      end















