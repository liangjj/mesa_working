*deck gnum3p.f
c***begin prologue     gnum3p
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
c***                   a(i-1) = beta, a(i) = -(alpha + beta ), a(i+1) = alpha
c***                   b(i-1) = beta*( alpha**2 + alpha*beta -beta**2 ) /12.
c***                   b(i)   = ( alpha +beta )* ( alpha**2 +3.*alpha*beta +
c***                                               beta**2 )/12.
c***                   b(i-1) = alpha* ( -alpha**2 + alpha*beta +beta**2 )/12.
c***routines called  
c***end prologue       gnum3p
      subroutine gnum3p(r,a,b,f,n)
      implicit integer (a-z)
      real*8  r, a, b, alpha, beta, f
      dimension r(n), a(n,-1:1), b(n,-1:1), f(n)
      common /io/ inp, iout
c     the coefficients at the first and last point should never be needed.
      do 10 i=2,n-1
         alpha=r(i)-r(i-1)
         beta=r(i+1)-r(i)
         a(i-1) = beta
         a(i) = -(alpha + beta )
         a(i+1) = alpha
         b(i-1) = beta*( alpha*alpha + alpha*beta -beta*beta ) /12.d0
         b(i)   = ( alpha +beta )*( alpha*alpha +3.d0*alpha*beta +
     1                                              beta*beta )/12.d0
         b(i-1) = alpha*( -alpha*alpha + alpha*beta +beta*beta )/12.d0
   10 continue
      return
      end







