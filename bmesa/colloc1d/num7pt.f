*deck num7pt.f
c***begin prologue     num7pt
c***date written       950720
c***revision date               (yymmdd)
c***keywords           num7pt
c***author             schneider, barry(nsf)
c***source             @(#)util
c***purpose            to set up the seven point numerov formula
c***                   for a matrix solution of the second order
c***                   one-dimensional schroedinger equation.
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis. it is well
c***                   described in Kopal's book on numerical analysis.
c***
c***                   if the equation to be solved is,           
c***                            y''  = g - f(x) y
c***
c***                   the basic numerov formula is:
c***                   {619.
c***                   {465.  + 23. *stp*stp*f(i+2)}*y(i+2)  +
c***                   {465.  + 23. *stp*stp*f(i-2)}*y(i-2)  +
c***                   {1920. + 688.*stp*stp*f(i+1)}*y(i+1)  +
c***                   {1920. + 688.*stp*stp*f(i-1)}*y(i-1)  -
c***                   {4770. -2358.*stp*stp*f(i)}  *y(i)
c***                                 =
c***                   23. *stp*stp*{ g(i-2) + g(i+2) }     +
c***                   688.*stp*stp*{ g(i+1) + g(i-1) }     +
c***                   2358.*stp*stp*g(i)  

c***                                         
c***
c***routines called  
c***end prologue       num7pt
      subroutine num7pt(band,f,stp,n)
      implicit integer (a-z)
      real *8  band, f, stp
      dimension band(n,-3:3), f(n)
      common /io/ inp, iout
c**********************************************************************c
c            calculate the diagonal, subdiagonal and superdiagonal     c
c            elements of the penta-diagonal matrix describing the      c      
c            numerov propagation of a second-order differential        c
c                                equation                              c
c**********************************************************************c
      do 10 i=1,n
         band(i,0) = - ( 4770. -2358.*stp*stp*f(i) )
   10 continue
      do 20 i=2,n
         band(i,-1) =  ( 1920. + 688.*stp*stp*f(i-1) )
   20 continue
      do 30 i=1,n-1
         band(i,1)  =  ( 1920. + 688.*stp*stp*f(i+1) )
   30 continue
      do 40 i=3,n
         band(i,-2) = (  465.  + 23. *stp*stp*f(i-2) )
   40 continue
      do 50 i=1,n-2
         band(i,2) =  (  465.  + 23. *stp*stp*f(i+2) )
   50 continue
      return
      end







