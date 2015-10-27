*deck dwnlsm
      subroutine dwnlsm (w, mdw, mme, ma, n, l, prgopt, x, rnorm, mode,
     +   ipivot, itype, wd, h, scale, z, temp, d)
c***begin prologue  dwnlsm
c***subsidiary
c***purpose  subsidiary to dwnnls
c***library   slatec
c***type      double precision (wnlsm-s, dwnlsm-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     this is a companion subprogram to dwnnls.
c     the documentation for dwnnls has complete usage instructions.
c
c     in addition to the parameters discussed in the prologue to
c     subroutine dwnnls, the following work arrays are used in
c     subroutine dwnlsm  (they are passed through the calling
c     sequence from dwnnls for purposes of variable dimensioning).
c     their contents will in general be of no interest to the user.
c
c     variables of type real are double precision.
c
c         ipivot(*)
c            an array of length n.  upon completion it contains the
c         pivoting information for the cols of w(*,*).
c
c         itype(*)
c            an array of length m which is used to keep track
c         of the classification of the equations.  itype(i)=0
c         denotes equation i as an equality constraint.
c         itype(i)=1 denotes equation i as a least squares
c         equation.
c
c         wd(*)
c            an array of length n.  upon completion it contains the
c         dual solution vector.
c
c         h(*)
c            an array of length n.  upon completion it contains the
c         pivot scalars of the householder transformations performed
c         in the case krank.lt.l.
c
c         scale(*)
c            an array of length m which is used by the subroutine
c         to store the diagonal matrix of weights.
c         these are used to apply the modified givens
c         transformations.
c
c         z(*),temp(*)
c            working arrays of length n.
c
c         d(*)
c            an array of length n that contains the
c         column scaling for the matrix (e).
c                                       (a)
c
c***see also  dwnnls
c***routines called  d1mach, dasum, daxpy, dcopy, dh12, dnrm2, drotm,
c                    drotmg, dscal, dswap, dwnlit, idamax, xermsg
c***revision history  (yymmdd)
c   790701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890618  completely restructured and revised.  (wrb & rwc)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  fixed an error message.  (rwc)
c   900604  dp version created from sp version.  (rwc)
c   900911  restriction on value of alamda included.  (wrb)
c***end prologue  dwnlsm
      integer ipivot(*), itype(*), l, ma, mdw, mme, mode, n
      double precision d(*), h(*), prgopt(*), rnorm, scale(*), temp(*),
     *   w(mdw,*), wd(*), x(*), z(*)
c
      external d1mach, dasum, daxpy, dcopy, dh12, dnrm2, drotm, drotmg,
     *   dscal, dswap, dwnlit, idamax, xermsg
      double precision d1mach, dasum, dnrm2
      integer idamax
c
      double precision alamda, alpha, alsq, amax, blowup, bnorm,
     *   dope(3), drelpr, eanorm, fac, sm, sparam(5), t, tau, wmax, z2,
     *   zz
      integer i, idope(3), imax, isol, itemp, iter, itmax, iwmax, j,
     *   jcon, jp, key, krank, l1, last, link, m, me, next, niv, nlink,
     *   nopt, nsoln, ntimes
      logical done, feasbl, first, hitcon, pos
c
      save drelpr, first
      data first /.true./
c***first executable statement  dwnlsm
c
c     initialize variables.
c     drelpr is the precision for the particular machine
c     being used.  this logic avoids resetting it every entry.
c
      if (first) drelpr = d1mach(4)
      first = .false.
c
c     set the nominal tolerance used in the code.
c
      tau = sqrt(drelpr)
c
      m = ma + mme
      me = mme
      mode = 2
c
c     to process option vector
c
      fac = 1.d-4
c
c     set the nominal blow up factor used in the code.
c
      blowup = tau
c
c     the nominal column scaling used in the code is
c     the identity scaling.
c
      call dcopy (n, 1.d0, 0, d, 1)
c
c     define bound for number of options to change.
c
      nopt = 1000
c
c     define bound for positive value of link.
c
      nlink = 100000
      ntimes = 0
      last = 1
      link = prgopt(1)
      if (link.le.0 .or. link.gt.nlink) then
         call xermsg ('slatec', 'dwnlsm',
     +      'in dwnnls, the option vector is undefined', 3, 1)
         return
      endif
c
  100 if (link.gt.1) then
         ntimes = ntimes + 1
         if (ntimes.gt.nopt) then
         call xermsg ('slatec', 'dwnlsm',
     +      'in dwnnls, the links in the option vector are cycling.',
     +      3, 1)
            return
         endif
c
         key = prgopt(last+1)
         if (key.eq.6 .and. prgopt(last+2).ne.0.d0) then
            do 110 j = 1,n
               t = dnrm2(m,w(1,j),1)
               if (t.ne.0.d0) t = 1.d0/t
               d(j) = t
  110       continue
         endif
c
         if (key.eq.7) call dcopy (n, prgopt(last+2), 1, d, 1)
         if (key.eq.8) tau = max(drelpr,prgopt(last+2))
         if (key.eq.9) blowup = max(drelpr,prgopt(last+2))
c
         next = prgopt(link)
         if (next.le.0 .or. next.gt.nlink) then
            call xermsg ('slatec', 'dwnlsm',
     +         'in dwnnls, the option vector is undefined', 3, 1)
            return
         endif
c
         last = link
         link = next
         go to 100
      endif
c
      do 120 j = 1,n
         call dscal (m, d(j), w(1,j), 1)
  120 continue
c
c     process option vector
c
      done = .false.
      iter = 0
      itmax = 3*(n-l)
      mode = 0
      nsoln = l
      l1 = min(m,l)
c
c     compute scale factor to apply to equality constraint equations.
c
      do 130 j = 1,n
         wd(j) = dasum(m,w(1,j),1)
  130 continue
c
      imax = idamax(n,wd,1)
      eanorm = wd(imax)
      bnorm = dasum(m,w(1,n+1),1)
      alamda = eanorm/(drelpr*fac)
c
c     on machines, such as the vaxes using d floating, with a very
c     limited exponent range for double precision values, the previously
c     computed value of alamda may cause an overflow condition.
c     therefore, this code further limits the value of alamda.
c
      alamda = min(alamda,sqrt(d1mach(2)))
c
c     define scaling diagonal matrix for modified givens usage and
c     classify equation types.
c
      alsq = alamda**2
      do 140 i = 1,m
c
c        when equation i is heavily weighted itype(i)=0,
c        else itype(i)=1.
c
         if (i.le.me) then
            t = alsq
            itemp = 0
         else
            t = 1.d0
            itemp = 1
         endif
         scale(i) = t
         itype(i) = itemp
  140 continue
c
c     set the solution vector x(*) to zero and the column interchange
c     matrix to the identity.
c
      call dcopy (n, 0.d0, 0, x, 1)
      do 150 i = 1,n
         ipivot(i) = i
  150 continue
c
c     perform initial triangularization in the submatrix
c     corresponding to the unconstrained variables.
c     set first l components of dual vector to zero because
c     these correspond to the unconstrained variables.
c
      call dcopy (l, 0.d0, 0, wd, 1)
c
c     the arrays idope(*) and dope(*) are used to pass
c     information to dwnlit().  this was done to avoid
c     a long calling sequence or the use of common.
c
      idope(1) = me
      idope(2) = nsoln
      idope(3) = l1
c
      dope(1) = alsq
      dope(2) = eanorm
      dope(3) = tau
      call dwnlit (w, mdw, m, n, l, ipivot, itype, h, scale, rnorm,
     +            idope, dope, done)
      me    = idope(1)
      krank = idope(2)
      niv   = idope(3)
c
c     perform wnnls algorithm using the following steps.
c
c     until(done)
c        compute search direction and feasible point
c        when (hitcon) add constraints
c        else perform multiplier test and drop a constraint
c        fin
c     compute-final-solution
c
c     to compute search direction and feasible point,
c     solve the triangular system of currently non-active
c     variables and store the solution in z(*).
c
c     to solve system
c     copy right hand side into temp vector to use overwriting method.
c
  160 if (done) go to 330
      isol = l + 1
      if (nsoln.ge.isol) then
         call dcopy (niv, w(1,n+1), 1, temp, 1)
         do 170 j = nsoln,isol,-1
            if (j.gt.krank) then
               i = niv - nsoln + j
            else
               i = j
            endif
c
            if (j.gt.krank .and. j.le.l) then
               z(j) = 0.d0
            else
               z(j) = temp(i)/w(i,j)
               call daxpy (i-1, -z(j), w(1,j), 1, temp, 1)
            endif
  170    continue
      endif
c
c     increment iteration counter and check against maximum number
c     of iterations.
c
      iter = iter + 1
      if (iter.gt.itmax) then
         mode = 1
         done = .true.
      endif
c
c     check to see if any constraints have become active.
c     if so, calculate an interpolation factor so that all
c     active constraints are removed from the basis.
c
      alpha = 2.d0
      hitcon = .false.
      do 180 j = l+1,nsoln
         zz = z(j)
         if (zz.le.0.d0) then
            t = x(j)/(x(j)-zz)
            if (t.lt.alpha) then
               alpha = t
               jcon = j
            endif
            hitcon = .true.
         endif
  180 continue
c
c     compute search direction and feasible point
c
      if (hitcon) then
c
c        to add constraints, use computed alpha to interpolate between
c        last feasible solution x(*) and current unconstrained (and
c        infeasible) solution z(*).
c
         do 190 j = l+1,nsoln
            x(j) = x(j) + alpha*(z(j)-x(j))
  190    continue
         feasbl = .false.
c
c        remove column jcon and shift columns jcon+1 through n to the
c        left.  swap column jcon into the n th position.  this achieves
c        upper hessenberg form for the nonactive constraints and
c        leaves an upper hessenberg matrix to retriangularize.
c
  200    do 210 i = 1,m
            t = w(i,jcon)
            call dcopy (n-jcon, w(i, jcon+1), mdw, w(i, jcon), mdw)
            w(i,n) = t
  210    continue
c
c        update permuted index vector to reflect this shift and swap.
c
         itemp = ipivot(jcon)
         do 220 i = jcon,n - 1
            ipivot(i) = ipivot(i+1)
  220    continue
         ipivot(n) = itemp
c
c        similarly permute x(*) vector.
c
         call dcopy (n-jcon, x(jcon+1), 1, x(jcon), 1)
         x(n) = 0.d0
         nsoln = nsoln - 1
         niv = niv - 1
c
c        retriangularize upper hessenberg matrix after adding
c        constraints.
c
         i = krank + jcon - l
         do 230 j = jcon,nsoln
            if (itype(i).eq.0 .and. itype(i+1).eq.0) then
c
c              zero ip1 to i in column j
c
               if (w(i+1,j).ne.0.d0) then
                  call drotmg (scale(i), scale(i+1), w(i,j), w(i+1,j),
     +                         sparam)
                  w(i+1,j) = 0.d0
                  call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw,
     +                        sparam)
               endif
            elseif (itype(i).eq.1 .and. itype(i+1).eq.1) then
c
c              zero ip1 to i in column j
c
               if (w(i+1,j).ne.0.d0) then
                  call drotmg (scale(i), scale(i+1), w(i,j), w(i+1,j),
     +                         sparam)
                  w(i+1,j) = 0.d0
                  call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw,
     +                        sparam)
               endif
            elseif (itype(i).eq.1 .and. itype(i+1).eq.0) then
               call dswap (n+1, w(i,1), mdw, w(i+1,1), mdw)
               call dswap (1, scale(i), 1, scale(i+1), 1)
               itemp = itype(i+1)
               itype(i+1) = itype(i)
               itype(i) = itemp
c
c              swapped row was formerly a pivot element, so it will
c              be large enough to perform elimination.
c              zero ip1 to i in column j.
c
               if (w(i+1,j).ne.0.d0) then
                  call drotmg (scale(i), scale(i+1), w(i,j), w(i+1,j),
     +                         sparam)
                  w(i+1,j) = 0.d0
                  call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw,
     +                        sparam)
               endif
            elseif (itype(i).eq.0 .and. itype(i+1).eq.1) then
               if (scale(i)*w(i,j)**2/alsq.gt.(tau*eanorm)**2) then
c
c                 zero ip1 to i in column j
c
                  if (w(i+1,j).ne.0.d0) then
                     call drotmg (scale(i), scale(i+1), w(i,j),
     +                            w(i+1,j), sparam)
                     w(i+1,j) = 0.d0
                     call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw,
     +                           sparam)
                  endif
               else
                  call dswap (n+1, w(i,1), mdw, w(i+1,1), mdw)
                  call dswap (1, scale(i), 1, scale(i+1), 1)
                  itemp = itype(i+1)
                  itype(i+1) = itype(i)
                  itype(i) = itemp
                  w(i+1,j) = 0.d0
               endif
            endif
            i = i + 1
  230    continue
c
c        see if the remaining coefficients in the solution set are
c        feasible.  they should be because of the way alpha was
c        determined.  if any are infeasible, it is due to roundoff
c        error.  any that are non-positive will be set to zero and
c        removed from the solution set.
c
         do 240 jcon = l+1,nsoln
            if (x(jcon).le.0.d0) go to 250
  240    continue
         feasbl = .true.
  250    if (.not.feasbl) go to 200
      else
c
c        to perform multiplier test and drop a constraint.
c
         call dcopy (nsoln, z, 1, x, 1)
         if (nsoln.lt.n) call dcopy (n-nsoln, 0.d0, 0, x(nsoln+1), 1)
c
c        reclassify least squares equations as equalities as necessary.
c
         i = niv + 1
  260    if (i.le.me) then
            if (itype(i).eq.0) then
               i = i + 1
            else
               call dswap (n+1, w(i,1), mdw, w(me,1), mdw)
               call dswap (1, scale(i), 1, scale(me), 1)
               itemp = itype(i)
               itype(i) = itype(me)
               itype(me) = itemp
               me = me - 1
            endif
            go to 260
         endif
c
c        form inner product vector wd(*) of dual coefficients.
c
         do 280 j = nsoln+1,n
            sm = 0.d0
            do 270 i = nsoln+1,m
               sm = sm + scale(i)*w(i,j)*w(i,n+1)
  270       continue
            wd(j) = sm
  280    continue
c
c        find j such that wd(j)=wmax is maximum.  this determines
c        that the incoming column j will reduce the residual vector
c        and be positive.
c
  290    wmax = 0.d0
         iwmax = nsoln + 1
         do 300 j = nsoln+1,n
            if (wd(j).gt.wmax) then
               wmax = wd(j)
               iwmax = j
            endif
  300    continue
         if (wmax.le.0.d0) go to 330
c
c        set dual coefficients to zero for incoming column.
c
         wd(iwmax) = 0.d0
c
c        wmax .gt. 0.d0, so okay to move column iwmax to solution set.
c        perform transformation to retriangularize, and test for near
c        linear dependence.
c
c        swap column iwmax into nsoln-th position to maintain upper
c        hessenberg form of adjacent columns, and add new column to
c        triangular decomposition.
c
         nsoln = nsoln + 1
         niv = niv + 1
         if (nsoln.ne.iwmax) then
            call dswap (m, w(1,nsoln), 1, w(1,iwmax), 1)
            wd(iwmax) = wd(nsoln)
            wd(nsoln) = 0.d0
            itemp = ipivot(nsoln)
            ipivot(nsoln) = ipivot(iwmax)
            ipivot(iwmax) = itemp
         endif
c
c        reduce column nsoln so that the matrix of nonactive constraints
c        variables is triangular.
c
         do 320 j = m,niv+1,-1
            jp = j - 1
c
c           when operating near the me line, test to see if the pivot
c           element is near zero.  if so, use the largest element above
c           it as the pivot.  this is to maintain the sharp interface
c           between weighted and non-weighted rows in all cases.
c
            if (j.eq.me+1) then
               imax = me
               amax = scale(me)*w(me,nsoln)**2
               do 310 jp = j - 1,niv,-1
                  t = scale(jp)*w(jp,nsoln)**2
                  if (t.gt.amax) then
                     imax = jp
                     amax = t
                  endif
  310          continue
               jp = imax
            endif
c
            if (w(j,nsoln).ne.0.d0) then
               call drotmg (scale(jp), scale(j), w(jp,nsoln),
     +                      w(j,nsoln), sparam)
               w(j,nsoln) = 0.d0
               call drotm (n+1-nsoln, w(jp,nsoln+1), mdw, w(j,nsoln+1),
     +                     mdw, sparam)
            endif
  320    continue
c
c        solve for z(nsoln)=proposed new value for x(nsoln).  test if
c        this is nonpositive or too large.  if this was true or if the
c        pivot term was zero, reject the column as dependent.
c
         if (w(niv,nsoln).ne.0.d0) then
            isol = niv
            z2 = w(isol,n+1)/w(isol,nsoln)
            z(nsoln) = z2
            pos = z2 .gt. 0.d0
            if (z2*eanorm.ge.bnorm .and. pos) then
               pos = .not. (blowup*z2*eanorm.ge.bnorm)
            endif
c
c           try to add row me+1 as an additional equality constraint.
c           check size of proposed new solution component.
c           reject it if it is too large.
c
         elseif (niv.le.me .and. w(me+1,nsoln).ne.0.d0) then
            isol = me + 1
            if (pos) then
c
c              swap rows me+1 and niv, and scale factors for these rows.
c
               call dswap (n+1, w(me+1,1), mdw, w(niv,1), mdw)
               call dswap (1, scale(me+1), 1, scale(niv), 1)
               itemp = itype(me+1)
               itype(me+1) = itype(niv)
               itype(niv) = itemp
               me = me + 1
            endif
         else
            pos = .false.
         endif
c
         if (.not.pos) then
            nsoln = nsoln - 1
            niv = niv - 1
         endif
         if (.not.(pos.or.done)) go to 290
      endif
      go to 160
c
c     else perform multiplier test and drop a constraint.  to compute
c     final solution.  solve system, store results in x(*).
c
c     copy right hand side into temp vector to use overwriting method.
c
  330 isol = 1
      if (nsoln.ge.isol) then
         call dcopy (niv, w(1,n+1), 1, temp, 1)
         do 340 j = nsoln,isol,-1
            if (j.gt.krank) then
               i = niv - nsoln + j
            else
               i = j
            endif
c
            if (j.gt.krank .and. j.le.l) then
               z(j) = 0.d0
            else
               z(j) = temp(i)/w(i,j)
               call daxpy (i-1, -z(j), w(1,j), 1, temp, 1)
            endif
  340    continue
      endif
c
c     solve system.
c
      call dcopy (nsoln, z, 1, x, 1)
c
c     apply householder transformations to x(*) if krank.lt.l
c
      if (krank.lt.l) then
         do 350 i = 1,krank
            call dh12 (2, i, krank+1, l, w(i,1), mdw, h(i), x, 1, 1, 1)
  350    continue
      endif
c
c     fill in trailing zeroes for constrained variables not in solution.
c
      if (nsoln.lt.n) call dcopy (n-nsoln, 0.d0, 0, x(nsoln+1), 1)
c
c     permute solution vector to natural order.
c
      do 380 i = 1,n
         j = i
  360    if (ipivot(j).eq.i) go to 370
         j = j + 1
         go to 360
c
  370    ipivot(j) = ipivot(i)
         ipivot(i) = j
         call dswap (1, x(j), 1, x(i), 1)
  380 continue
c
c     rescale the solution using the column scaling.
c
      do 390 j = 1,n
         x(j) = x(j)*d(j)
  390 continue
c
      do 400 i = nsoln+1,m
         t = w(i,n+1)
         if (i.le.me) t = t/alamda
         t = (scale(i)*t)*t
         rnorm = rnorm + t
  400 continue
c
      rnorm = sqrt(rnorm)
      return
      end
