*deck @(#)series.f	5.1 11/6/94
c***begin prologue     series.f
c***date written       880528   (yymmdd)
c***revision date      @(#)series.f	5.1
c***keywords           basint, gaussian integrals, series
c***author             schneider, barry (lanl)
c***source             @(#)series.f	5.1 11/6/94
c***purpose            series needed for calculation of integrals in basint
c***description        series used for small arguments in basint
c***                   
c***                   
c
c***references         bis notes
c
c***routines called    
c***end prologue       series.f
      function series(t,n)
      implicit integer (a-z)
      real *8 series, fac1, fac2, sum, term
      real *8 t
      common /io/ inp, iout
c------------------------------------------------------------------------c
c                    series used                                         c
c          1. +2.*t**2/(k+3) + (2*t**2)**2/(k+3)*(k+5)                   c
c                      + ........                                        c
c                                                                        c
c------------------------------------------------------------------------c
      fac1=2.d+00*t*t
      fac2=fac1
      fac3=n+3
      sum=1.d+00
      do 10 i=2,100
         term=fac1/fac3
         if(term.lt.1.d-10) go to 20
         sum=sum+term
         fac1=fac1*fac2
         fac3=fac3*(n+1+i+i)
   10 continue
      write (iout,30) term
   30 format(/,5x,'no convergence:last term in sum',1x,d15.8)
      stop
   20 series=sum
      return
      end
