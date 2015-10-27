*deck dsoseq
      subroutine dsoseq (fnc, n, s, rtolx, atolx, tolf, iflag, mxit,
     +   ncjs, nsrrc, nsri, iprint, fmax, c, nc, b, p, temp, x, y, fac,
     +   is)
c***begin prologue  dsoseq
c***subsidiary
c***purpose  subsidiary to dsos
c***library   slatec
c***type      double precision (soseqs-s, dsoseq-d)
c***author  (unknown)
c***description
c
c     dsoseq solves a system of n simultaneous nonlinear equations.
c     see the comments in the interfacing routine dsos for a more
c     detailed description of some of the items in the calling list.
c
c **********************************************************************
c   -input-
c
c     fnc- function subprogram which evaluates the equations
c     n  -number of equations
c     s  -solution vector of initial guesses
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
c     iprint-internal printing parameter. you must set iprint=-1 if you
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
c     s    -solution vector
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
c          last iteration loop upon exit from dsos.
c     y   -array containing the solution increments.
c     fac -array containing factors used in computing numerical
c          derivatives.
c     is  -records the pivotal information (column interchanges)
c
c **********************************************************************
c *** three machine dependent parameters appear in this subroutine.
c
c *** the smallest positive magnitude, zero, is defined by the function
c *** routine d1mach(1).
c
c *** uro, the computer unit roundoff value, is defined by d1mach(3) for
c *** machines that round or d1mach(4) for machines that truncate.
c *** uro is the smallest positive number such that 1.+uro  .gt.  1.
c
c *** the output tape unit number, loun, is defined by the function
c *** i1mach(2).
c **********************************************************************
c
c***see also  dsos
c***routines called  d1mach, dsossl, i1mach
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dsoseq
c
c
      integer i1mach
      double precision d1mach
      integer ic, icr, iflag, iprint, is(*), isj, isv, it, item, itry,
     1     j, jk, js, k, kd, kj, kk, km1, kn, ksv, l, loun, ls, m, mit,
     2     mm, mxit, n, nc, ncjs, np1, nsri, nsrrc
      double precision atolx, b(*), c(*), csv, f, fac(*), fact, fdif,
     1     fmax, fmin, fmxs, fn1, fn2, fnc, fp, h, hx, p(*), pmax, re,
     2     rtolx, s(*), sruro, temp(*), test, tolf, uro, x(*), xnorm,
     3     y(*), yj, yn1, yn2, yn3, ynorm, yns, zero
c
c     begin block permitting ...exits to 430
c        begin block permitting ...exits to 410
c           begin block permitting ...exits to 390
c***first executable statement  dsoseq
               uro = d1mach(4)
               loun = i1mach(2)
               zero = d1mach(1)
               re = max(rtolx,uro)
               sruro = sqrt(uro)
c
               iflag = 0
               np1 = n + 1
               icr = 0
               ic = 0
               itry = ncjs
               yn1 = 0.0d0
               yn2 = 0.0d0
               yn3 = 0.0d0
               yns = 0.0d0
               mit = 0
               fn1 = 0.0d0
               fn2 = 0.0d0
               fmxs = 0.0d0
c
c              initialize the interchange (pivoting) vector and
c              save the current solution approximation for future use.
c
               do 10 k = 1, n
                  is(k) = k
                  x(k) = s(k)
                  temp(k) = x(k)
   10          continue
c
c
c              *********************************************************
c              **** begin principal iteration loop  ****
c              *********************************************************
c
               do 380 m = 1, mxit
c                 begin block permitting ...exits to 350
c                    begin block permitting ...exits to 240
c
                        do 20 k = 1, n
                           fac(k) = sruro
   20                   continue
c
   30                   continue
c                          begin block permitting ...exits to 180
                              kn = 1
                              fmax = 0.0d0
c
c
c                             ******** begin subiteration loop defining
c                             the linearization of each ********
c                             equation which results in the construction
c                             of an upper ******** triangular matrix
c                             approximating the forward ********
c                             triangularization of the full jacobian
c                             matrix
c
                              do 170 k = 1, n
c                                begin block permitting ...exits to 160
                                    km1 = k - 1
c
c                                   back-solve a triangular linear
c                                   system obtaining improved solution
c                                   values for k-1 of the variables from
c                                   the first k-1 equations. these
c                                   variables are then eliminated from
c                                   the k-th equation.
c
                                    if (km1 .eq. 0) go to 50
                                       call dsossl(k,n,km1,y,c,b,kn)
                                       do 40 j = 1, km1
                                          js = is(j)
                                          x(js) = temp(js) + y(j)
   40                                  continue
   50                               continue
c
c
c                                   evaluate the k-th equation and the
c                                   intermediate computation for the max
c                                   norm of the residual vector.
c
                                    f = fnc(x,k)
                                    fmax = max(fmax,abs(f))
c
c                                   if we wish to perform several
c                                   iterations using a fixed
c                                   factorization of an approximate
c                                   jacobian,we need only update the
c                                   constant vector.
c
c                                ...exit
                                    if (itry .lt. ncjs) go to 160
c
c
                                    it = 0
c
c                                   compute partial derivatives that are
c                                   required in the linearization of the
c                                   k-th reduced equation
c
                                    do 90 j = k, n
                                       item = is(j)
                                       hx = x(item)
                                       h = fac(item)*hx
                                       if (abs(h) .le. zero)
     1                                    h = fac(item)
                                       x(item) = hx + h
                                       if (km1 .eq. 0) go to 70
                                          y(j) = h
                                          call dsossl(k,n,j,y,c,b,kn)
                                          do 60 l = 1, km1
                                             ls = is(l)
                                             x(ls) = temp(ls) + y(l)
   60                                     continue
   70                                  continue
                                       fp = fnc(x,k)
                                       x(item) = hx
                                       fdif = fp - f
                                       if (abs(fdif) .gt. uro*abs(f))
     1                                    go to 80
                                          fdif = 0.0d0
                                          it = it + 1
   80                                  continue
                                       p(j) = fdif/h
   90                               continue
c
                                    if (it .le. (n - k)) go to 110
c
c                                      all computed partial derivatives
c                                      of the k-th equation are
c                                      effectively zero.try larger
c                                      perturbations of the independent
c                                      variables.
c
                                       do 100 j = k, n
                                          isj = is(j)
                                          fact = 100.0d0*fac(isj)
c           ..............................exit
                                          if (fact .gt. 1.0d10)
     1                                       go to 390
                                          fac(isj) = fact
  100                                  continue
c                          ............exit
                                       go to 180
  110                               continue
c
c                                ...exit
                                    if (k .eq. n) go to 160
c
c                                   achieve a pivoting effect by
c                                   choosing the maximal derivative
c                                   element
c
                                    pmax = 0.0d0
                                    do 130 j = k, n
                                       test = abs(p(j))
                                       if (test .le. pmax) go to 120
                                          pmax = test
                                          isv = j
  120                                  continue
  130                               continue
c           ........................exit
                                    if (pmax .eq. 0.0d0) go to 390
c
c                                   set up the coefficients for the k-th
c                                   row of the triangular linear system
c                                   and save the partial derivative of
c                                   largest magnitude
c
                                    pmax = p(isv)
                                    kk = kn
                                    do 140 j = k, n
                                       if (j .ne. isv)
     1                                    c(kk) = -p(j)/pmax
                                       kk = kk + 1
  140                               continue
                                    p(k) = pmax
c
c
c                                ...exit
                                    if (isv .eq. k) go to 160
c
c                                   interchange the two columns of c
c                                   determined by the pivotal strategy
c
                                    ksv = is(k)
                                    is(k) = is(isv)
                                    is(isv) = ksv
c
                                    kd = isv - k
                                    kj = k
                                    do 150 j = 1, k
                                       csv = c(kj)
                                       jk = kj + kd
                                       c(kj) = c(jk)
                                       c(jk) = csv
                                       kj = kj + n - j
  150                               continue
  160                            continue
c
                                 kn = kn + np1 - k
c
c                                store the components for the constant
c                                vector
c
                                 b(k) = -f/p(k)
c
  170                         continue
c                       ......exit
                              go to 190
  180                      continue
                        go to 30
  190                   continue
c
c                       ********
c                       ******** end of loop creating the triangular
c                       linearization matrix
c                       ********
c
c
c                        solve the resulting triangular system for a new
c                        solution approximation and obtain the solution
c                        increment norm.
c
                        kn = kn - 1
                        y(n) = b(n)
                        if (n .gt. 1) call dsossl(n,n,n,y,c,b,kn)
                        xnorm = 0.0d0
                        ynorm = 0.0d0
                        do 200 j = 1, n
                           yj = y(j)
                           ynorm = max(ynorm,abs(yj))
                           js = is(j)
                           x(js) = temp(js) + yj
                           xnorm = max(xnorm,abs(x(js)))
  200                   continue
c
c
c                       print intermediate solution iterates and
c                       residual norm if desired
c
                        if (iprint .ne. (-1)) go to 220
                           mm = m - 1
                           write (loun,210) fmax,mm,(x(j), j = 1, n)
  210                      format ('0residual norm =', d9.2, / 1x,
     1                             'solution iterate (', i3, ')', /
     2                             (1x, 5d26.14))
  220                   continue
c
c                       test for convergence to a solution (relative
c                       and/or absolute error comparison on successive
c                       approximations of each solution variable)
c
                        do 230 j = 1, n
                           js = is(j)
c                    ......exit
                           if (abs(y(j)) .gt. re*abs(x(js)) + atolx)
     1                        go to 240
  230                   continue
                        if (fmax .le. fmxs) iflag = 1
  240                continue
c
c                    test for convergence to a solution based on
c                    residuals
c
                     if (fmax .le. tolf) iflag = iflag + 2
c        ............exit
                     if (iflag .gt. 0) go to 410
c
c
                     if (m .gt. 1) go to 250
                        fmin = fmax
                     go to 330
  250                continue
c                       begin block permitting ...exits to 320
c
c                          save solution having minimum residual norm.
c
                           if (fmax .ge. fmin) go to 270
                              mit = m + 1
                              yn1 = ynorm
                              yn2 = yns
                              fn1 = fmxs
                              fmin = fmax
                              do 260 j = 1, n
                                 s(j) = x(j)
  260                         continue
                              ic = 0
  270                      continue
c
c                          test for limiting precision convergence. very
c                          slowly convergent problems may also be
c                          detected.
c
                           if (ynorm .gt. sruro*xnorm) go to 290
                           if (fmax .lt. 0.2d0*fmxs
     1                         .or. fmax .gt. 5.0d0*fmxs) go to 290
                           if (ynorm .lt. 0.2d0*yns
     1                         .or. ynorm .gt. 5.0d0*yns) go to 290
                              icr = icr + 1
                              if (icr .ge. nsrrc) go to 280
                                 ic = 0
c                       .........exit
                                 go to 320
  280                         continue
                              iflag = 4
                              fmax = fmin
c     ........................exit
                              go to 430
  290                      continue
                           icr = 0
c
c                          test for divergence of the iterative scheme.
c
                           if (ynorm .gt. 2.0d0*yns
     1                         .or. fmax .gt. 2.0d0*fmxs) go to 300
                              ic = 0
                           go to 310
  300                      continue
                              ic = ic + 1
c                       ......exit
                              if (ic .lt. nsri) go to 320
                              iflag = 7
c        .....................exit
                              go to 410
  310                      continue
  320                   continue
  330                continue
c
c                    check to see if next iteration can use the old
c                    jacobian factorization
c
                     itry = itry - 1
                     if (itry .eq. 0) go to 340
                     if (20.0d0*ynorm .gt. xnorm) go to 340
                     if (ynorm .gt. 2.0d0*yns) go to 340
c                 ......exit
                        if (fmax .lt. 2.0d0*fmxs) go to 350
  340                continue
                     itry = ncjs
  350             continue
c
c                 save the current solution approximation and the
c                 residual and solution increment norms for use in the
c                 next iteration.
c
                  do 360 j = 1, n
                     temp(j) = x(j)
  360             continue
                  if (m .ne. mit) go to 370
                     fn2 = fmax
                     yn3 = ynorm
  370             continue
                  fmxs = fmax
                  yns = ynorm
c
c
  380          continue
c
c              *********************************************************
c              **** end of principal iteration loop ****
c              *********************************************************
c
c
c               too many iterations, convergence was not achieved.
               m = mxit
               iflag = 5
               if (yn1 .gt. 10.0d0*yn2 .or. yn3 .gt. 10.0d0*yn1)
     1            iflag = 6
               if (fn1 .gt. 5.0d0*fmin .or. fn2 .gt. 5.0d0*fmin)
     1            iflag = 6
               if (fmax .gt. 5.0d0*fmin) iflag = 6
c        ......exit
               go to 410
  390       continue
c
c
c           a jacobian-related matrix is effectively singular.
            iflag = 8
            do 400 j = 1, n
               s(j) = temp(j)
  400       continue
c     ......exit
            go to 430
  410    continue
c
c
         do 420 j = 1, n
            s(j) = x(j)
  420    continue
  430 continue
c
c
      mxit = m
      return
      end
