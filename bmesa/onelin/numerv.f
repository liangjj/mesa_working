c $Header: numerv.f,v 1.2 92/12/12 09:34:57 bis Exp $
*deck numerv.f
c***begin prologue     numerv
c***date written       910803
c***revision date               (yymmdd)
c***keywords           numerov
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            to set up the three band numerov formula
c***                   for a matrix solution of the second order
c***                   one-dimensional schroedinger equation.
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis. it is well
c***                   described in Kopal's book on numerical analysis.
c***                                         ''
c***                   the equation here is y  + f(x) y = g
c***
c***                   first derivative formula needed for 
c***                   boundary condition is:
c***                   ( -3*y(0) + 4*y(1) -y(2) )/(2*stp)
c***routines called  
c***end prologue
      subroutine numerv(diag,sudiag,spdiag,f,s4,stp,m,last,nfd)
      implicit integer (a-z)
      real *8  diag, sudiag, spdiag, f, stp, s4
      dimension diag(0:m), sudiag(0:m), spdiag(0:m), f(0:m), s4(4)
      common /io/ inp, iout
c**********************************************************************c
c            calculate the diagonal, subdiagonal and superdiagonal     c
c            elements of the tri-diagonal matrix describing the        c      
c            numerov propagation of a second-order differential        c
c                                equation                              c
c            the m that is passed is enough larger than the n          c
c            of the finite difference equation that there should be    c
c            no problem with the matrix elements.                      c  
c**********************************************************************c
      do 10 i=1,m-1
         diag(i) = - ( 24.d+00 -10.d+00*stp*stp*f(i) )
   10 continue
      do 20 i=1,m
         sudiag(i) = ( 12.d+00 + stp*stp*f(i-1) )
   20 continue
      do 30 i=1,m-1
         spdiag(i) = ( 12.d+00 + stp*stp*f(i+1) )
   30 continue
      nfd=last-1
      s4(1)=sudiag(nfd)
      s4(2)=diag(nfd)
      s4(3)=spdiag(nfd)
      return
      end







