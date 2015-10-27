*deck bisect
      subroutine bisect (n, eps1, d, e, e2, lb, ub, mm, m, w, ind, ierr,
     +   rv4, rv5)
c***begin prologue  bisect
c***purpose  compute the eigenvalues of a symmetric tridiagonal matrix
c            in a given interval using sturm sequencing.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (bisect-s)
c***keywords  eigenvalues, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the bisection technique
c     in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix which lie in a specified interval,
c     using bisection.
c
c     on input
c
c        n is the order of the matrix.  n is an integer variable.
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix.
c          eps1 is a real variable.
c
c        d contains the diagonal elements of the input matrix.
c          d is a one-dimensional real array, dimensioned d(n).
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c          e is a one-dimensional real array, dimensioned e(n).
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.  e2 is a one-dimensional real array,
c          dimensioned e2(n).
c
c        lb and ub define the interval to be searched for eigenvalues.
c          if lb is not less than ub, no eigenvalues will be found.
c          lb and ub are real variables.
c
c        mm should be set to an upper bound for the number of
c          eigenvalues in the interval.  warning - if more than
c          mm eigenvalues are determined to lie in the interval,
c          an error return is made with no eigenvalues found.
c          mm is an integer variable.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        m is the number of eigenvalues determined to lie in (lb,ub).
c          m is an integer variable.
c
c        w contains the m eigenvalues in ascending order.
c          w is a one-dimensional real array, dimensioned w(mm).
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c          ind is an one-dimensional integer array, dimensioned ind(mm).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          3*n+1      if m exceeds mm.  in this case, m contains the
c                     number of eigenvalues determined to lie in
c                     (lb,ub).
c
c        rv4 and rv5 are one-dimensional real arrays used for temporary
c          storage, dimensioned rv4(n) and rv5(n).
c
c     the algol procedure sturmcnt contained in tristurm
c     appears in bisect in-line.
c
c     note that subroutine tql1 or imtql1 is generally faster than
c     bisect, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  r1mach
c***revision history  (yymmdd)
c   760101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bisect
c
      integer i,j,k,l,m,n,p,q,r,s,ii,mm,m1,m2,tag,ierr,isturm
      real d(*),e(*),e2(*),w(*),rv4(*),rv5(*)
      real u,v,lb,t1,t2,ub,xu,x0,x1,eps1,machep,s1,s2
      integer ind(*)
      logical first
c
      save first, machep
      data first /.true./
c***first executable statement  bisect
      if (first) then
         machep = r1mach(4)
      endif
      first = .false.
c
      ierr = 0
      tag = 0
      t1 = lb
      t2 = ub
c     .......... look for small sub-diagonal entries ..........
      do 40 i = 1, n
         if (i .eq. 1) go to 20
         s1 = abs(d(i)) + abs(d(i-1))
         s2 = s1 + abs(e(i))
         if (s2 .gt. s1) go to 40
   20    e2(i) = 0.0e0
   40 continue
c     .......... determine the number of eigenvalues
c                in the interval ..........
      p = 1
      q = n
      x1 = ub
      isturm = 1
      go to 320
   60 m = s
      x1 = lb
      isturm = 2
      go to 320
   80 m = m - s
      if (m .gt. mm) go to 980
      q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0e0
c
      do 120 q = p, n
         x1 = u
         u = 0.0e0
         v = 0.0e0
         if (q .eq. n) go to 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
         x0 = max(d(q)+(x1+u),x0)
         if (v .eq. 0.0e0) go to 140
  120 continue
c
  140 x1 = max(abs(xu),abs(x0)) * machep
      if (eps1 .le. 0.0e0) eps1 = -x1
      if (p .ne. q) go to 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q-p+1)
      lb = max(t1,xu-x1)
      ub = min(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5e0
         s1 = 2.0e0*(abs(xu) + abs(x0) + abs(eps1))
         s2 = s1 + abs(x0 - xu)
         if (s2 .eq. s1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0e0
c
         do 340 i = p, q
            if (u .ne. 0.0e0) go to 325
            v = abs(e(i)) / machep
            if (e2(i) .eq. 0.0e0) v = 0.0e0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0e0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     .......... order eigenvalues tagged with their
c                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- underestimate of number of
c                eigenvalues in interval ..........
  980 ierr = 3 * n + 1
 1001 lb = t1
      ub = t2
      return
      end
