*deck fdjac3
      subroutine fdjac3 (fcn, m, n, x, fvec, fjac, ldfjac, iflag,
     +   epsfcn, wa)
c***begin prologue  fdjac3
c***subsidiary
c***purpose  subsidiary to snls1 and snls1e
c***library   slatec
c***type      single precision (fdjac3-s, dfdjc3-d)
c***author  (unknown)
c***description
c
c     this subroutine computes a forward-difference approximation
c     to the m by n jacobian matrix associated with a specified
c     problem of m functions in n variables.
c
c     the subroutine statement is
c
c       subroutine fdjac3(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(iflag,m,n,x,fvec,fjac,ldfjac)
c         integer ldfjac,m,n,iflag
c         real x(n),fvec(m),fjac(ldfjac,n)
c         ----------
c         when iflag.eq.1 calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac3.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       iflag is an integer variable which can be used to terminate
c         the execution of fdjac3. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c***see also  snls1, snls1e
c***routines called  r1mach
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c***end prologue  fdjac3
      integer m,n,ldfjac,iflag
      real epsfcn
      real x(*),fvec(*),fjac(ldfjac,*),wa(*)
      integer i,j
      real eps,epsmch,h,temp,zero
      real r1mach
      save zero
      data zero /0.0e0/
c***first executable statement  fdjac3
      epsmch = r1mach(4)
c
      eps = sqrt(max(epsfcn,epsmch))
c      set iflag=1 to indicate that function values
c           are to be returned by fcn.
      iflag = 1
      do 20 j = 1, n
         temp = x(j)
         h = eps*abs(temp)
         if (h .eq. zero) h = eps
         x(j) = temp + h
         call fcn(iflag,m,n,x,wa,fjac,ldfjac)
         if (iflag .lt. 0) go to 30
         x(j) = temp
         do 10 i = 1, m
            fjac(i,j) = (wa(i) - fvec(i))/h
   10       continue
   20    continue
   30 continue
      return
c
c     last card of subroutine fdjac3.
c
      end
