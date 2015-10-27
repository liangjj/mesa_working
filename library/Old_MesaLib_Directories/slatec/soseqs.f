*deck soseqs
      subroutine soseqs (fnc, n, s, rtolx, atolx, tolf, iflag, mxit,
     +   ncjs, nsrrc, nsri, iprint, fmax, c, nc, b, p, temp, x, y, fac,
     +   is)
c***begin prologue  soseqs
c***subsidiary
c***purpose  subsidiary to sos
c***library   slatec
c***type      single precision (soseqs-s, dsoseq-d)
c***author  (unknown)
c***description
c
c     soseqs solves a system of n simultaneous nonlinear equations.
c     see the comments in the interfacing routine sos for a more
c     detailed description of some of the items in the calling list.
c
c ********************************************************************
c
c   -input-
c     fnc -function subprogram which evaluates the equations
c     n   -number of equations
c     s   -solution vector of initial guesses
c     rtolx-relative error tolerance on solution components
c     atolx-absolute error tolerance on solution components
c     tolf-residual error tolerance
c     mxit-maximum number of allowable iterations.
c     ncjs-maximum number of consecutive iterative steps to perform
c          using the same triangular jacobian matrix approximation.
c     nsrrc-number of consecutive iterative steps for which the
c          limiting precision accuracy test must be satisfied
c          before the routine exits with iflag=4.
c     nsri-number of consecutive iterative steps for which the
c          diverging condition test must be satisfied before
c          the routine exits with iflag=7.
c     iprint-internal printing parameter.  you must set iprint=-1 if you
c          want the intermediate solution iterates and a residual norm
c          to be printed.
c     c   -internal work array, dimensioned at least n*(n+1)/2.
c     nc  -dimension of c array. nc  .ge.  n*(n+1)/2.
c     b   -internal work array, dimensioned n.
c     p   -internal work array, dimensioned n.
c     temp-internal work array, dimensioned n.
c     x   -internal work array, dimensioned n.
c     y   -internal work array, dimensioned n.
c     fac -internal work array, dimensioned n.
c     is  -internal work array, dimensioned n.
c
c   -output-
c     s   -solution vector
c     iflag-status indicator flag
c     mxit-the actual number of iterations performed
c     fmax-residual norm
c     c   -upper unit triangular matrix which approximates the
c          forward triangularization of the full jacobian matrix.
c          stored in a vector with dimension at least n*(n+1)/2.
c     b   -contains the residuals (function values) divided
c          by the corresponding components of the p vector
c     p   -array used to store the partial derivatives. after
c          each iteration p(k) contains the maximal derivative
c          occurring in the k-th reduced equation.
c     temp-array used to store the previous solution iterate.
c     x   -solution vector. contains the values achieved on the
c          last iteration loop upon exit from sos.
c     y   -array containing the solution increments.
c     fac -array containing factors used in computing numerical
c          derivatives.
c     is  -records the pivotal information (column interchanges)
c
c **********************************************************************
c *** three machine dependent parameters appear in this subroutine.
c
c *** the smallest positive magnitude, zero, is defined by the function
c *** routine r1mach(1).
c
c *** uro, the computer unit roundoff value, is defined by r1mach(3) for
c *** machines that round or r1mach(4) for machines that truncate.
c *** uro is the smallest positive number such that 1.+uro  .gt.  1.
c
c *** the output tape unit number, loun, is defined by the function
c *** i1mach(2).
c **********************************************************************
c
c***see also  sos
c***routines called  i1mach, r1mach, sossol
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  soseqs
c
c
      dimension s(*), c(nc), b(*), is(*), p(*), temp(*), x(*), y(*),
     1          fac(*)
c
c***first executable statement  soseqs
      uro = r1mach(4)
      loun = i1mach(2)
      zero = r1mach(1)
      re = max(rtolx,uro)
      sruro = sqrt(uro)
c
      iflag = 0
      np1 = n + 1
      icr = 0
      ic = 0
      itry = ncjs
      yn1 = 0.
      yn2 = 0.
      yn3 = 0.
      yns = 0.
      mit = 0
      fn1 = 0.
      fn2 = 0.
      fmxs = 0.
c
c     initialize the interchange (pivoting) vector and
c     save the current solution approximation for future use.
c
      do 10 k=1,n
        is(k) = k
        x(k) = s(k)
        temp(k) = x(k)
   10 continue
c
c
c    *****************************************
c    **** begin principal iteration loop  ****
c    *****************************************
c
      do 330 m=1,mxit
c
        do 20 k=1,n
          fac(k) = sruro
   20   continue
c
   30   kn = 1
        fmax = 0.
c
c
c    ******** begin subiteration loop defining the linearization of each
c    ******** equation which results in the construction of an upper
c    ******** triangular matrix approximating the forward
c    ******** triangularization of the full jacobian matrix
c
        do 170 k=1,n
          km1 = k - 1
c
c     back-solve a triangular linear system obtaining
c     improved solution values for k-1 of the variables
c     from the first k-1 equations. these variables are then
c     eliminated from the k-th equation.
c
          if (km1 .eq. 0) go to 50
          call sossol(k, n, km1, y, c, b, kn)
          do 40 j=1,km1
            js = is(j)
            x(js) = temp(js) + y(j)
   40     continue
c
c
c     evaluate the k-th equation and the intermediate computation
c     for the max norm of the residual vector.
c
   50     f = fnc(x,k)
          fmax = max(fmax,abs(f))
c
c     if we wish to perform several iterations using a fixed
c     factorization of an approximate jacobian,we need only
c     update the constant vector.
c
          if (itry .lt. ncjs) go to 160
c
c
          it = 0
c
c     compute partial derivatives that are required in the linearization
c     of the k-th reduced equation
c
          do 90 j=k,n
            item = is(j)
            hx = x(item)
            h = fac(item)*hx
            if (abs(h) .le. zero) h = fac(item)
            x(item) = hx + h
            if (km1 .eq. 0) go to 70
            y(j) = h
            call sossol(k, n, j, y, c, b, kn)
            do 60 l=1,km1
              ls = is(l)
              x(ls) = temp(ls) + y(l)
   60       continue
   70       fp = fnc(x,k)
            x(item) = hx
            fdif = fp - f
            if (abs(fdif) .gt. uro*abs(f)) go to 80
            fdif = 0.
            it = it + 1
   80       p(j) = fdif/h
   90     continue
c
          if (it .le. (n-k)) go to 110
c
c     all computed partial derivatives of the k-th equation
c     are effectively zero.try larger perturbations of the
c     independent variables.
c
          do 100 j=k,n
            isj = is(j)
            fact = 100.*fac(isj)
            if (fact .gt. 1.e+10) go to 340
            fac(isj) = fact
  100     continue
          go to 30
c
  110     if (k .eq. n) go to 160
c
c     achieve a pivoting effect by choosing the maximal derivative
c     element
c
          pmax = 0.
          do 120 j=k,n
            test = abs(p(j))
            if (test .le. pmax) go to 120
            pmax = test
            isv = j
  120     continue
          if (pmax .eq. 0.) go to 340
c
c     set up the coefficients for the k-th row of the triangular
c     linear system and save the partial derivative of
c     largest magnitude
c
          pmax = p(isv)
          kk = kn
          do 140 j=k,n
            if (j .eq. isv) go to 130
            c(kk) = -p(j)/pmax
  130       kk = kk + 1
  140     continue
          p(k) = pmax
c
c
          if (isv .eq. k) go to 160
c
c     interchange the two columns of c determined by the
c     pivotal strategy
c
          ksv = is(k)
          is(k) = is(isv)
          is(isv) = ksv
c
          kd = isv - k
          kj = k
          do 150 j=1,k
            csv = c(kj)
            jk = kj + kd
            c(kj) = c(jk)
            c(jk) = csv
            kj = kj + n - j
  150     continue
c
  160     kn = kn + np1 - k
c
c     store the components for the constant vector
c
          b(k) = -f/p(k)
c
  170   continue
c
c    ********
c    ******** end of loop creating the triangular linearization matrix
c    ********
c
c
c     solve the resulting triangular system for a new solution
c     approximation and obtain the solution increment norm.
c
        kn = kn - 1
        y(n) = b(n)
        if (n .gt. 1) call sossol(n, n, n, y, c, b, kn)
        xnorm = 0.
        ynorm = 0.
        do 180 j=1,n
          yj = y(j)
          ynorm = max(ynorm,abs(yj))
          js = is(j)
          x(js) = temp(js) + yj
          xnorm = max(xnorm,abs(x(js)))
  180   continue
c
c
c     print intermediate solution iterates and residual norm if desired
c
        if (iprint.ne.(-1)) go to 190
        mm = m - 1
        write (loun,1234) fmax, mm, (x(j),j=1,n)
 1234   format ('0residual norm =', e9.2, /1x, 'solution iterate',
     1   ' (', i3, ')', /(1x, 5e26.14))
  190   continue
c
c     test for convergence to a solution (relative and/or absolute error
c     comparison on successive approximations of each solution variable)
c
        do 200 j=1,n
          js = is(j)
          if (abs(y(j)) .gt. re*abs(x(js))+atolx) go to 210
  200   continue
        if (fmax .le. fmxs) iflag = 1
c
c     test for convergence to a solution based on residuals
c
  210   if (fmax .gt. tolf) go to 220
        iflag = iflag + 2
  220   if (iflag .gt. 0) go to 360
c
c
        if (m .gt. 1) go to 230
        fmin = fmax
        go to 280
c
c     save solution having minimum residual norm.
c
  230   if (fmax .ge. fmin) go to 250
        mit = m + 1
        yn1 = ynorm
        yn2 = yns
        fn1 = fmxs
        fmin = fmax
        do 240 j=1,n
          s(j) = x(j)
  240   continue
        ic = 0
c
c     test for limiting precision convergence.  very slowly convergent
c     problems may also be detected.
c
  250   if (ynorm .gt. sruro*xnorm) go to 260
        if ((fmax .lt. 0.2*fmxs) .or. (fmax .gt. 5.*fmxs)) go to 260
        if ((ynorm .lt. 0.2*yns) .or. (ynorm .gt. 5.*yns)) go to 260
        icr = icr + 1
        if (icr .lt. nsrrc) go to 270
        iflag = 4
        fmax = fmin
        go to 380
  260   icr = 0
c
c     test for divergence of the iterative scheme.
c
        if ((ynorm .le. 2.*yns) .and. (fmax .le. 2.*fmxs)) go to 270
        ic = ic + 1
        if (ic .lt. nsri) go to 280
        iflag = 7
        go to 360
  270   ic = 0
c
c     check to see if next iteration can use the old jacobian
c     factorization
c
  280   itry = itry - 1
        if (itry .eq. 0) go to 290
        if (20.*ynorm .gt. xnorm) go to 290
        if (ynorm .gt. 2.*yns) go to 290
        if (fmax .lt. 2.*fmxs) go to 300
  290   itry = ncjs
c
c     save the current solution approximation and the residual and
c     solution increment norms for use in the next iteration.
c
  300   do 310 j=1,n
          temp(j) = x(j)
  310   continue
        if (m.ne.mit) go to 320
        fn2 = fmax
        yn3 = ynorm
  320   fmxs = fmax
        yns = ynorm
c
c
  330 continue
c
c    *****************************************
c    **** end of principal iteration loop ****
c    *****************************************
c
c
c     too many iterations, convergence was not achieved.
      m = mxit
      iflag = 5
      if (yn1 .gt. 10.0*yn2 .or. yn3 .gt. 10.0*yn1) iflag = 6
      if (fn1 .gt. 5.0*fmin .or. fn2 .gt. 5.0*fmin) iflag = 6
      if (fmax .gt. 5.0*fmin) iflag = 6
      go to 360
c
c
c     a jacobian-related matrix is effectively singular.
  340 iflag = 8
      do 350 j=1,n
        s(j) = temp(j)
  350 continue
      go to 380
c
c
  360 do 370 j=1,n
        s(j) = x(j)
  370 continue
c
c
  380 mxit = m
      return
      end
