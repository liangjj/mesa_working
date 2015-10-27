*deck num3pt.f
c***begin prologue     num3pt
c***date written       950720
c***revision date               (yymmdd)
c***keywords           num3pt
c***author             schneider, barry(nsf)
c***source             @(#)util
c***purpose            to set up the three point numerov formula
c***                   for a matrix solution of the second order
c***                   one-dimensional schroedinger equation.
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis. it is well
c***                   described in Kopal's book on numerical analysis.
c***
c***                   the basic numerov formula is:
c***                   y(i-1) -2.*y(i) +y(i+1) =
c***                          stp*stp*[ y''(i-1) + 10.*y''(i) +y''(i+1) ] /12.
c***                   if the equation to be solved is,           
c***                            y''  = g - f(x) y
c***                   we substitute in for the second derivatives to get,
c***                   [ 12. +     stp*stp*f(i-1) ]*y(i-1)
c***                 - [ 24. - 10.*stp*stp*f(i)   ]*y(i) 
c***                 + [ 12. +     stp*stp*f(i+1) ]*y(i+1)
c***                                  =
c***                               stp*stp*[ g(i-1) +10.*g(i) +g(i+1) ]       
c***                                         
c***
c***routines called  
c***end prologue       num3pt
      subroutine num3pt(band,f,stp,n)
      implicit integer (a-z)
      real*8  band, f, stp
      dimension band(n,-1:1), f(n)
      common /io/ inp, iout
c**********************************************************************c
c            calculate the diagonal, subdiagonal and superdiagonal     c
c            elements of the tri-diagonal matrix describing the        c      
c            numerov propagation of a second-order differential        c
c                                equation                              c
c**********************************************************************c
      do 10 i=1,n
         band(i,0) = - ( 24.d+00 -10.d+00*stp*stp*f(i) )
   10 continue
      do 20 i=2,n
         band(i,-1) = ( 12.d+00 + stp*stp*f(i-1) )
   20 continue
      do 30 i=1,n-1
         band(i,1) = ( 12.d+00 + stp*stp*f(i+1) )
   30 continue
      return
      end







