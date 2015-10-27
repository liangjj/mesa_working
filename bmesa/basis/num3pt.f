*deck num3pt.f
c***begin prologue     num3pt
c***date written       950720
c***revision date               (yymmdd)
c***keywords           num3pt
c***author             schneider, barry(nsf)
c***source             @(#)util
c***purpose            to set up the three point generalized
c***                   numerov formula for a matrix solution of the 
c***                   second order one-dimensional schroedinger equation.
c***
c***references         this is described in m. seaton, methods in
c***                   computational physics.
c***
c***                   the generalized numerov formula does not assume
c***                   uniform step size.  it is given as:
c***                   a(i-1)*y(i-1) +a(i)*y(i) + a(i+1)*y(i+1) =
c***                   [ b(i-1)*y''(i-1) + b(i)*y''(i) + b(i+1)*y''(i+1) ]
c***                   if the equation to be solved is,           
c***                            y''  = g - f(x) y
c***                   we substitute in for the second derivatives to get,
c***                   the final working formula.
c***                                         
c***                   the construction of the a's and b's depends on the
c***                   points.  we define,
c***                   alpha = r(i)-r(i-1) and beta =  r(i+1)-r(i)
c***                   and then get,
c***                   a(i-1) = beta
c***                   a(i)   = -(alpha + beta )
c***                   a(i+1) = alpha
c***                   b(i-1) = beta*( alpha**2 + alpha*beta -beta**2 ) /12.
c***                   b(i)   = ( alpha +beta )* ( alpha**2 +3.*alpha*beta +
c***                                               beta**2 )/12.
c***                   b(i+1) = alpha* ( -alpha**2 + alpha*beta +beta**2 )/12.
c***routines called  
c***end prologue       num3pt
      subroutine num3pt(bandy,bandg,r,f,n,fixed)
      implicit integer (a-z)
      real*8  r, bandy, bandg, alpha, beta, f, stp
      logical fixed
      dimension r(*), bandy(n,-1:1), bandg(n,-1:1), f(*)
      common /io/ inp, iout
c     the coefficients at the first and last point should never be needed.
      if (.not.fixed) then
          do 10 i=1,n
             alpha=r(i)-r(i-1)
             beta=r(i+1)-r(i)
             bandy(i,-1) = beta
             bandy(i,0) = -(alpha + beta )
             bandy(i,1) = alpha
             bandg(i,-1) = beta*( alpha*alpha + alpha*beta -beta*beta )
     1                                        /12.d0
             bandg(i,0)  = ( alpha +beta )*( alpha*alpha +
     2                       3.d0*alpha*beta +beta*beta )/12.d0
             bandg(i,1) = alpha*( -alpha*alpha + alpha*beta +beta*beta )
     1                                        /12.d0
             bandy(i,-1) = bandy(i,-1) + bandg(i,-1) * f(i-1)
             bandy(i,0)  = bandy(i,0)  + bandg(i,0)  * f(i)
             bandy(i,1)  = bandy(i,1)  + bandg(i,1)  * f(i+1)
   10     continue
      else
          stp=r(2)-r(1)
          do 20 i=1,n
             bandy(i,0) = - ( 24.d+00 -10.d+00*stp*stp*f(i) )
             bandg(i,0) = 10.d0*stp*stp
   20     continue
          do 30 i=2,n
             bandy(i,-1) = ( 12.d+00 + stp*stp*f(i-1) )
             bandg(i,-1)= stp*stp
   30     continue
          do 40 i=1,n-1
             bandy(i,1) = ( 12.d+00 + stp*stp*f(i+1) )
             bandg(i,1)= stp*stp
   40     continue
      endif
      return
      end







