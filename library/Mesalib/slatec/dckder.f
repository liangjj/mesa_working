*deck dckder
      subroutine dckder (m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode,
     +   err)
c***begin prologue  dckder
c***purpose  check the gradients of m nonlinear functions in n
c            variables, evaluated at a point x, for consistency
c            with the functions themselves.
c***library   slatec
c***category  f3, g4c
c***type      double precision (chkder-s, dckder-d)
c***keywords  gradients, jacobian, minpack, nonlinear
c***author  hiebert, k. l. (snla)
c***description
c
c   this subroutine is a companion routine to dnsq and dnsqe. it may
c   be used to check the coding of the jacobian calculation.
c
c     subroutine dckder
c
c     this subroutine checks the gradients of m nonlinear functions
c     in n variables, evaluated at a point x, for consistency with
c     the functions themselves. the user must call dckder twice,
c     first with mode = 1 and then with mode = 2.
c
c     mode = 1. on input, x must contain the point of evaluation.
c               on output, xp is set to a neighboring point.
c
c     mode = 2. on input, fvec must contain the functions and the
c                         rows of fjac must contain the gradients
c                         of the respective functions each evaluated
c                         at x, and fvecp must contain the functions
c                         evaluated at xp.
c               on output, err contains measures of correctness of
c                          the respective gradients.
c
c     the subroutine does not perform reliably if cancellation or
c     rounding errors cause a severe loss of significance in the
c     evaluation of a function. therefore, none of the components
c     of x should be unusually small (in particular, zero) or any
c     other value which may cause loss of significance.
c
c     the subroutine statement is
c
c       subroutine dckder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables.
c
c       x is an input array of length n.
c
c       fvec is an array of length m. on input when mode = 2,
c         fvec must contain the functions evaluated at x.
c
c       fjac is an m by n array. on input when mode = 2,
c         the rows of fjac must contain the gradients of
c         the respective functions evaluated at x.
c
c       ldfjac is a positive integer input parameter not less than m
c         which specifies the leading dimension of the array fjac.
c
c       xp is an array of length n. on output when mode = 1,
c         xp is set to a neighboring point of x.
c
c       fvecp is an array of length m. on input when mode = 2,
c         fvecp must contain the functions evaluated at xp.
c
c       mode is an integer input variable set to 1 on the first call
c         and 2 on the second. other values of mode are equivalent
c         to mode = 1.
c
c       err is an array of length m. on output when mode = 2,
c         err contains measures of correctness of the respective
c         gradients. if there is no severe loss of significance,
c         then if err(i) is 1.0 the i-th gradient is correct,
c         while if err(i) is 0.0 the i-th gradient is incorrect.
c         for values of err between 0.0 and 1.0, the categorization
c         is less certain. in general, a value of err(i) greater
c         than 0.5 indicates that the i-th gradient is probably
c         correct, while a value of err(i) less than 0.5 indicates
c         that the i-th gradient is probably incorrect.
c
c***references  m. j. d. powell, a hybrid method for nonlinear equa-
c                 tions. in numerical methods for nonlinear algebraic
c                 equations, p. rabinowitz, editor.  gordon and breach,
c                 1988.
c***routines called  d1mach
c***revision history  (yymmdd)
c   800301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dckder
      integer i, j, ldfjac, m, mode, n
      double precision d1mach, eps, epsf, epslog, epsmch, err(*),
     1     factor, fjac(ldfjac,*), fvec(*), fvecp(*), one, temp, x(*),
     2     xp(*), zero
      save factor, one, zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
c
c     epsmch is the machine precision.
c
c***first executable statement  dckder
      epsmch = d1mach(4)
c
      eps = sqrt(epsmch)
c
      if (mode .eq. 2) go to 20
c
c        mode = 1.
c
         do 10 j = 1, n
            temp = eps*abs(x(j))
            if (temp .eq. zero) temp = eps
            xp(j) = x(j) + temp
   10       continue
         go to 70
   20 continue
c
c        mode = 2.
c
         epsf = factor*epsmch
         epslog = log10(eps)
         do 30 i = 1, m
            err(i) = zero
   30       continue
         do 50 j = 1, n
            temp = abs(x(j))
            if (temp .eq. zero) temp = one
            do 40 i = 1, m
               err(i) = err(i) + temp*fjac(i,j)
   40          continue
   50       continue
         do 60 i = 1, m
            temp = one
            if (fvec(i) .ne. zero .and. fvecp(i) .ne. zero
     1          .and. abs(fvecp(i)-fvec(i)) .ge. epsf*abs(fvec(i)))
     2         temp = eps*abs((fvecp(i)-fvec(i))/eps-err(i))
     3                /(abs(fvec(i)) + abs(fvecp(i)))
            err(i) = one
            if (temp .gt. epsmch .and. temp .lt. eps)
     1         err(i) = (log10(temp) - epslog)/epslog
            if (temp .ge. eps) err(i) = zero
   60       continue
   70 continue
c
      return
c
c     last card of subroutine dckder.
c
      end
