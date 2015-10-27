h05459
s 00027/00000/00000
d D 1.1 94/02/16 20:34:59 mesa 1 0
c date and time created 94/02/16 20:34:59 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
      function gbslve(shift, n, a, b)
c
c       this procedure performs elimination to solve for the
c       n-th component of the solution delta to the equation
c
c             (jn - shift*identity) * delta  = en,
c
c       where en is the vector of all zeroes except for 1 in
c       the n-th position.
c
c       the matrix jn is symmetric tridiagonal, with diagonal
c       elements a(i), off-diagonal elements b(i).  this equation
c       must be solved to obtain the appropriate changes in the lower
c       2 by 2 submatrix of coefficients for orthogonal polynomials.
c
c
      implicit real*8 (a-h,o-z)
      dimension  a(n),b(n)
c
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      gbslve = 1.0d0  /alpha
      return
      end
E 1
