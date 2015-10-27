*deck forwrd
      subroutine forwrd(f,wt,int0,intgl,n)
c***begin prologue     forwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            forward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on [a,b] as [a,r(1)], [a,r(2)]
c                      [a,r(3)]....[a,b] where r(i) are the quadrature points in
c                      the interval. the integral is initialized as intgl(0)
c                      which is either its last value or zero depending on the
c                      interval. intgl(1)=[a,r1] and intgl(n-1)=[a,b]   
c               
c***references
c
c***routines called
c
c***end prologue       forwrd
c
      implicit integer (a-z)
      dimension f(n), wt(n,n-1), intgl(1:n-1)
      common /io/ inp, iout
      real*8 f, wt, intgl, int0, sumf
      sumf=int0
      do 10 i=1,n-1
         do 20 j=1,n
            sumf=sumf+wt(j,i)*f(j)
   20    continue
         intgl(i)=sumf
   10 continue
      int0=intgl(n-1)
      return
      end















