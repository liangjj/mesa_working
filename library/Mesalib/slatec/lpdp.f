*deck lpdp
      subroutine lpdp (a, mda, m, n1, n2, prgopt, x, wnorm, mode, ws,
     +   is)
c***begin prologue  lpdp
c***subsidiary
c***purpose  subsidiary to lsei
c***library   slatec
c***type      single precision (lpdp-s, dlpdp-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     dimension a(mda,n+1),prgopt(*),x(n),ws((m+2)*(n+7)),is(m+n+1),
c     where n=n1+n2.  this is a slight overestimate for ws(*).
c
c     determine an n1-vector w, and
c               an n2-vector z
c     which minimizes the euclidean length of w
c     subject to g*w+h*z .ge. y.
c     this is the least projected distance problem, lpdp.
c     the matrices g and h are of respective
c     dimensions m by n1 and m by n2.
c
c     called by subprogram lsi( ).
c
c     the matrix
c                (g h y)
c
c     occupies rows 1,...,m and cols 1,...,n1+n2+1 of a(*,*).
c
c     the solution (w) is returned in x(*).
c                  (z)
c
c     the value of mode indicates the status of
c     the computation after returning to the user.
c
c          mode=1  the solution was successfully obtained.
c
c          mode=2  the inequalities are inconsistent.
c
c***see also  lsei
c***routines called  scopy, sdot, snrm2, sscal, wnnls
c***revision history  (yymmdd)
c   790701  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated the author section.  (wrb)
c***end prologue  lpdp
c
c     subroutines called
c
c     wnnls         solves a nonnegatively constrained linear least
c                   squares problem with linear equality constraints.
c                   part of this package.
c
c++
c     sdot,         subroutines from the blas package.
c     sscal,snrm2,  see trans. math. soft., vol. 5, no. 3, p. 308.
c     scopy
c
      real             a(mda,*), prgopt(*), ws(*), wnorm, x(*)
      integer is(*)
      real             fac, one, rnorm, sc, ynorm, zero
      real             sdot, snrm2
      save zero, one, fac
      data zero, one /0.e0,1.e0/, fac /0.1e0/
c***first executable statement  lpdp
      n = n1 + n2
      mode = 1
      if (.not.(m.le.0)) go to 20
      if (.not.(n.gt.0)) go to 10
      x(1) = zero
      call scopy(n, x, 0, x, 1)
   10 wnorm = zero
      return
   20 np1 = n + 1
c
c     scale nonzero rows of inequality matrix to have length one.
      do 40 i=1,m
        sc = snrm2(n,a(i,1),mda)
        if (.not.(sc.ne.zero)) go to 30
        sc = one/sc
        call sscal(np1, sc, a(i,1), mda)
   30   continue
   40 continue
c
c     scale rt.-side vector to have length one (or zero).
      ynorm = snrm2(m,a(1,np1),1)
      if (.not.(ynorm.ne.zero)) go to 50
      sc = one/ynorm
      call sscal(m, sc, a(1,np1), 1)
c
c     scale cols of matrix h.
   50 j = n1 + 1
   60 if (.not.(j.le.n)) go to 70
      sc = snrm2(m,a(1,j),1)
      if (sc.ne.zero) sc = one/sc
      call sscal(m, sc, a(1,j), 1)
      x(j) = sc
      j = j + 1
      go to 60
   70 if (.not.(n1.gt.0)) go to 130
c
c     copy transpose of (h g y) to work array ws(*).
      iw = 0
      do 80 i=1,m
c
c     move col of transpose of h into work array.
        call scopy(n2, a(i,n1+1), mda, ws(iw+1), 1)
        iw = iw + n2
c
c     move col of transpose of g into work array.
        call scopy(n1, a(i,1), mda, ws(iw+1), 1)
        iw = iw + n1
c
c     move component of vector y into work array.
        ws(iw+1) = a(i,np1)
        iw = iw + 1
   80 continue
      ws(iw+1) = zero
      call scopy(n, ws(iw+1), 0, ws(iw+1), 1)
      iw = iw + n
      ws(iw+1) = one
      iw = iw + 1
c
c     solve eu=f subject to (transpose of h)u=0, u.ge.0.  the
c     matrix e = transpose of (g y), and the (n+1)-vector
c     f = transpose of (0,...,0,1).
      ix = iw + 1
      iw = iw + m
c
c     do not check lengths of work arrays in this usage of wnnls( ).
      is(1) = 0
      is(2) = 0
      call wnnls(ws, np1, n2, np1-n2, m, 0, prgopt, ws(ix), rnorm,
     1 modew, is, ws(iw+1))
c
c     compute the components of the soln denoted above by w.
      sc = one - sdot(m,a(1,np1),1,ws(ix),1)
      if (.not.(one+fac*abs(sc).ne.one .and. rnorm.gt.zero)) go to 110
      sc = one/sc
      do 90 j=1,n1
        x(j) = sc*sdot(m,a(1,j),1,ws(ix),1)
   90 continue
c
c     compute the vector q=y-gw.  overwrite y with this vector.
      do 100 i=1,m
        a(i,np1) = a(i,np1) - sdot(n1,a(i,1),mda,x,1)
  100 continue
      go to 120
  110 mode = 2
      return
  120 continue
  130 if (.not.(n2.gt.0)) go to 180
c
c     copy transpose of (h q) to work array ws(*).
      iw = 0
      do 140 i=1,m
        call scopy(n2, a(i,n1+1), mda, ws(iw+1), 1)
        iw = iw + n2
        ws(iw+1) = a(i,np1)
        iw = iw + 1
  140 continue
      ws(iw+1) = zero
      call scopy(n2, ws(iw+1), 0, ws(iw+1), 1)
      iw = iw + n2
      ws(iw+1) = one
      iw = iw + 1
      ix = iw + 1
      iw = iw + m
c
c     solve rv=s subject to v.ge.0.  the matrix r =(transpose
c     of (h q)), where q=y-gw.  the (n2+1)-vector s =(transpose
c     of (0,...,0,1)).
c
c     do not check lengths of work arrays in this usage of wnnls( ).
      is(1) = 0
      is(2) = 0
      call wnnls(ws, n2+1, 0, n2+1, m, 0, prgopt, ws(ix), rnorm, modew,
     1 is, ws(iw+1))
c
c     compute the components of the soln denoted above by z.
      sc = one - sdot(m,a(1,np1),1,ws(ix),1)
      if (.not.(one+fac*abs(sc).ne.one .and. rnorm.gt.zero)) go to 160
      sc = one/sc
      do 150 j=1,n2
        l = n1 + j
        x(l) = sc*sdot(m,a(1,l),1,ws(ix),1)*x(l)
  150 continue
      go to 170
  160 mode = 2
      return
  170 continue
c
c     account for scaling of rt.-side vector in solution.
  180 call sscal(n, ynorm, x, 1)
      wnorm = snrm2(n1,x,1)
      return
      end
