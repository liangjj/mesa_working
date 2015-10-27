*deck dplint
      subroutine dplint (n, x, y, c)
c***begin prologue  dplint
c***purpose  produce the polynomial which interpolates a set of discrete
c            data points.
c***library   slatec
c***category  e1b
c***type      double precision (polint-s, dplint-d)
c***keywords  polynomial interpolation
c***author  huddleston, r. e., (snll)
c***description
c
c     abstract
c        subroutine dplint is designed to produce the polynomial which
c     interpolates the data  (x(i),y(i)), i=1,...,n.  dplint sets up
c     information in the array c which can be used by subroutine dpolvl
c     to evaluate the polynomial and its derivatives and by subroutine
c     dpolcf to produce the coefficients.
c
c     formal parameters
c     *** all type real variables are double precision ***
c     n  - the number of data points  (n .ge. 1)
c     x  - the array of abscissas (all of which must be distinct)
c     y  - the array of ordinates
c     c  - an array of information used by subroutines
c     *******  dimensioning information  *******
c     arrays x,y, and c must be dimensioned at least n in the calling
c     program.
c
c***references  l. f. shampine, s. m. davenport and r. e. huddleston,
c                 curve fitting by polynomials in one variable, report
c                 sla-74-0270, sandia laboratories, june 1974.
c***routines called  xermsg
c***revision history  (yymmdd)
c   740601  date written
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dplint
      integer i,k,km1,n
      double precision dif,c(*),x(*),y(*)
c***first executable statement  dplint
      if (n .le. 0) go to 91
      c(1)=y(1)
      if(n .eq. 1) return
      do 10010 k=2,n
      c(k)=y(k)
      km1=k-1
      do 10010 i=1,km1
c     check for distinct x values
      dif = x(i)-x(k)
      if (dif .eq. 0.0) go to 92
      c(k) = (c(i)-c(k))/dif
10010 continue
      return
   91 call xermsg ('slatec', 'dplint', 'n is zero or negative.', 2, 1)
      return
   92 call xermsg ('slatec', 'dplint',
     +   'the abscissas are not distinct.', 2, 1)
      return
      end
