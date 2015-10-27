*deck tridib
      subroutine tridib (n, eps1, d, e, e2, lb, ub, m11, m, w, ind,
     +   ierr, rv4, rv5)
c***begin prologue  tridib
c***purpose  compute the eigenvalues of a symmetric tridiagonal matrix
c            in a given interval using sturm sequencing.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (tridib-s)
c***keywords  eigenvalues of a real symmetric matrix, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure bisect,
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix between specified boundary indices,
c     using bisection.
c
c     on input
c
c        n is the order of the matrix.  n is an integer variable.
c
c        eps1 is an absolute error tolerance for the computed eigen-
c          values.  if the input eps1 is non-positive, it is reset for
c          each submatrix to a default value, namely, minus the product
c          of the relative machine precision and the 1-norm of the
c          submatrix.  eps1 is a real variable.
c
c        d contains the diagonal elements of the symmetric tridiagonal
c          matrix.  d is a one-dimensional real array, dimensioned d(n).
c
c        e contains the subdiagonal elements of the symmetric
c          tridiagonal matrix in its last n-1 positions.  e(1) is
c          arbitrary.  e is a one-dimensional real array, dimensioned
c          e(n).
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.  e2 is a one-dimensional real array,
c          dimensioned e2(n).
c
c        m11 specifies the lower boundary index for the set of desired
c          eigenvalues.  m11 is an integer variable.
c
c        m specifies the number of eigenvalues desired.  the upper
c          boundary index m22 is then obtained as m22=m11+m-1.
c          m is an integer variable.
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
c        lb and ub define an interval containing exactly the desired
c          eigenvalues.  lb and ub are real variables.
c
c        w contains, in its first m positions, the eigenvalues
c          between indices m11 and m22 in ascending order.
c          w is a one-dimensional real array, dimensioned w(m).
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c          ind is an one-dimensional integer array, dimensioned ind(m).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          3*n+1      if multiple eigenvalues at index m11 make
c                     unique selection of lb impossible,
c          3*n+2      if multiple eigenvalues at index m22 make
c                     unique selection of ub impossible.
c
c        rv4 and rv5 are one-dimensional real arrays used for temporary
c          storage of the lower and upper bounds for the eigenvalues in
c          the bisection process.  rv4 and rv5 are dimensioned rv4(n)
c          and rv5(n).
c
c     note that subroutine tql1, imtql1, or tqlrat is generally faster
c     than tridib, if more than n/4 eigenvalues are to be found.
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
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  tridib
c
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      real*8 d(*),e(*),e2(*),w(*),rv4(*),rv5(*)
      real*8 u,v,lb,t1,t2,ub,xu,x0,x1,eps1,machep,s1,s2
      integer ind(*)
      logical first
c
      save first, machep
      data first /.true./
c***first executable statement  tridib
      if (first) then
         machep = r1mach(4)
      endif
      first = .false.
c
      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0d0
c     .......... look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues ..........
      do 40 i = 1, n
         x1 = u
         u = 0.0d0
         if (i .ne. n) u = abs(e(i+1))
         xu = min(d(i)-(x1+u),xu)
         x0 = max(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         s1 = abs(d(i)) + abs(d(i-1))
         s2 = s1 + abs(e(i))
         if (s2 .gt. s1) go to 40
   20    e2(i) = 0.0d0
   40 continue
c
      x1 = max(abs(xu),abs(x0)) * machep * n
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     .......... determine an interval containing exactly
c                the desired eigenvalues ..........
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5d0
      if (x1 .eq. v) go to 980
      go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0d0
c
      do 120 q = p, n
         x1 = u
         u = 0.0d0
         v = 0.0d0
         if (q .eq. n) go to 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
         x0 = max(d(q)+(x1+u),x0)
         if (v .eq. 0.0d0) go to 140
  120 continue
c
  140 x1 = max(abs(xu),abs(x0)) * machep
      if (eps1 .le. 0.0d0) eps1 = -x1
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
  300    x1 = (xu + x0) * 0.5d0
         s1 = abs(xu) + abs(x0) + abs(eps1)
         s2 = s1 + abs(x0-xu)/2.0d0
         if (s2 .eq. s1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0d0
c
         do 340 i = p, q
            if (u .ne. 0.0d0) go to 325
            v = abs(e(i)) / machep
            if (e2(i) .eq. 0.0d0) v = 0.0d0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0d0) s = s + 1
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
c     .......... set error -- interval cannot be found containing
c                exactly the desired eigenvalues ..........
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
      end
