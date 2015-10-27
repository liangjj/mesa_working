*deck bkwrd
      subroutine bkwrd(f,wt,int0,intgl,n)
c***begin prologue     bkwrd
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            backward integration of indefinite integral
c***
c***description        performs the set of indefinite integrals which result
c***                   from integrating a function on [a,b] as [a,r(1)], [a,r(2)]
c                      [a,r(3)]....[a,b] where r(i) are the quadrature points in
c                      the interval. the integral is initialized as intgl(n)
c                      which is either its last value or zero depending on the
c                      interval. intgl(n-1)=[r(n-1),b]   
c                                             
c***references
c
c***routines called
c
c***end prologue       bkwrd
c
      implicit integer (a-z)
      dimension f(n),wt(n,n-1), intgl(1:n-1)
      common /io/ inp, iout
      real*8 f, wt, intgl, int0, sumb
      sumb=int0
      do 10 i=n-1,1,-1
         do 20 j=1,n
            sumb=sumb+wt(j,i)*f(j)
   20    continue
         intgl(i)=sumb
   10 continue
      int0=intgl(1)
      return
      end















