*deck num5pt.f
c***begin prologue     num5pt
c***date written       950720
c***revision date               (yymmdd)
c***keywords           num5pt
c***author             schneider, barry(nsf)
c***source             @(#)util
c***purpose            to set up the generalized five point numerov 
c***                   formula for a matrix solution of the second order
c***                   one-dimensional schroedinger equation.
c***
c***                   the generalized numerov formula does not assume
c***                   uniform step size.  it is given as:
c***                   a(i-2)*y(i-2) + a(i-1)*y(i-1) + a(i)*y(i) +
c***                                   a(i+1)*y(i+1) + a(i+2)*y(i+2) =
c***                   b(i-2)*y''(i-2) + b(i-1)*y''(i-1) + b(i)*y''(i) +
c***                                   + b(i+1)*y''(i+1) + b(i+2)*y''(i+2)
c***
c***                   if the equation to be solved is,           
c***                            y''  = g - f(x) y
c***                   we substitute in for the second derivatives to get
c***                   the final working formula.
c***
c***                   the construction of the a's and b's depends on the
c***                   points.  we define,
c***                   alpha = r(i)-r(i-1)  beta  =  r(i+1)-r(i)
c***                   gamma =  r(i)-r(i-2) delta =  r(i+2)-r(i)
c***                   and then get,
c***                   a(i-2) = 
c***                   a(i-1) = 
c***                   a(i)   = 
c***                   a(i+1) =
c***                   a(i+2) =
c***                   b(i-2) = 
c***                   b(i-1) = 
c***                   b(i)   = 
c***                   b(i+1) =
c***                   b(i+2) =
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis. it is well
c***                   described in Kopal's book on numerical analysis.
c***
c***routines called  
c***end prologue       num5pt
      subroutine num5pt(bandy,bandg,r,f,n,fixed)
      implicit integer (a-z)
      real*8  bandy, bandg, r, f
      logical fixed
      dimension bandy(n,-2:2), bandg(n,-2:2), r(n+1), f(n+1)
      common /io/ inp, iout
      if (.not.fixed) then
          do 10 i=2,n-1
             alpha=r(i)-r(i-1)
             beta=r(i+1)-r(i)
             gamma=r(i)-r(i-2)
             delta=r(i+2)-r(i)            
c             bandy(i,-2) =
c             bandy(i,-1) =
c             bandy(i,0) = 
c             bandy(i,1) = 
c             bandy(i,2) =
c             bandg(i,-2) = 
c             bandg(i,-1) = 
c             bandg(i,0)  = 
c             bandg(i,1) =
c             bandy(i,2) =
   10     continue
      else
          stp=r(2)-r(1)
          do 20 i=1,n
             bandy(i,0) = - ( 4770. -2358.*stp*stp*f(i) )
             bandg(i,0) = 2358.d0*stp*stp
   20     continue
          do 30 i=2,n
             bandy(i,-1) =  ( 1920. + 688.*stp*stp*f(i-1) )
             bandg(i,-1) = 688.d0*stp*stp
   30     continue
          do 40 i=1,n-1
             bandy(i,1)  =  ( 1920. + 688.*stp*stp*f(i+1) )
             bandg(i,1) = 688.d0*stp*stp
   40     continue
          do 50 i=3,n
             bandy(i,-2) = (  465.  + 23. *stp*stp*f(i-2) )
             bandg(i,-2) = 23.d0*stp*stp
   50     continue
          do 60 i=1,n-2
             bandy(i,2) =  (  465.  + 23. *stp*stp*f(i+2) )
             bandg(i,2) = 23.d0*stp*stp
   60     continue
      endif
      return
      end







