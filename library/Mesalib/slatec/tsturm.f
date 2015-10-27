*deck tsturm
      subroutine tsturm (nm, n, eps1, d, e, e2, lb, ub, mm, m, w, z,
     +   ierr, rv1, rv2, rv3, rv4, rv5, rv6)
c***begin prologue  tsturm
c***purpose  find those eigenvalues of a symmetric tridiagonal matrix
c            in a given interval and their associated eigenvectors by
c            sturm sequencing.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (tsturm-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix which lie in a specified interval and their
c     associated eigenvectors, using bisection and inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, z, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        eps1 is an absolute error tolerance for the computed eigen-
c          values.  it should be chosen so that the accuracy of these
c          eigenvalues is commensurate with relative perturbations of
c          the order of the relative machine precision in the matrix
c          elements.  if the input eps1 is non-positive, it is reset
c          for each submatrix to a default value, namely, minus the
c          product of the relative machine precision and the 1-norm of
c          the submatrix.  eps1 is a real variable.
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
c        lb and ub define the interval to be searched for eigenvalues.
c          if lb is not less than ub, no eigenvalues will be found.
c          lb and ub are real variables.
c
c        mm should be set to an upper bound for the number of
c          eigenvalues in the interval.  mm is an integer variable.
c          warning -  if more than mm eigenvalues are determined to lie
c          in the interval, an error return is made with no values or
c          vectors found.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded as
c          negligible, have been replaced by zero causing the matrix to
c          split into a direct sum of submatrices.  e2(1) is also set
c          to zero.
c
c        m is the number of eigenvalues determined to lie in (lb,ub).
c          m is an integer variable.
c
c        w contains the m eigenvalues in ascending order if the matrix
c          does not split.  if the matrix splits, the eigenvalues are
c          in ascending order for each submatrix.  if a vector error
c          exit is made, w contains those values already found.  w is a
c          one-dimensional real array, dimensioned w(mm).
c
c        z contains the associated set of orthonormal eigenvectors.
c          if an error exit is made, z contains those vectors already
c          found.  z is a one-dimensional real array, dimensioned
c          z(nm,mm).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          3*n+1      if m exceeds mm no eigenvalues or eigenvectors
c                     are computed,
c          4*n+j      if the eigenvector corresponding to the j-th
c                     eigenvalue fails to converge in 5 iterations, then
c                     the eigenvalues and eigenvectors in w and z should
c                     be correct for indices 1, 2, ..., j-1.
c
c        rv1, rv2, rv3, rv4, rv5, and rv6 are temporary storage arrays,
c          dimensioned rv1(n), rv2(n), rv3(n), rv4(n), rv5(n), and
c          rv6(n).
c
c     the algol procedure sturmcnt contained in tristurm
c     appears in tsturm in-line.
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
c***end prologue  tsturm
c
      integer i,j,k,m,n,p,q,r,s,ii,ip,jj,mm,m1,m2,nm,its
      integer ierr,group,isturm
      real d(*),e(*),e2(*),w(*),z(nm,*)
      real rv1(*),rv2(*),rv3(*),rv4(*),rv5(*),rv6(*)
      real u,v,lb,t1,t2,ub,uk,xu,x0,x1,eps1,eps2,eps3,eps4
      real norm,machep,s1,s2
      logical first
c
      save first, machep
      data first /.true./
c***first executable statement  tsturm
      if (first) then
         machep = r1mach(4)
      endif
      first = .false.
c
      ierr = 0
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
      r = r + 1
c
      do 160 i = 1, n
  160 z(i,r) = 0.0e0
c
      w(r) = d(p)
      z(p,r) = 1.0e0
      go to 940
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
c     .......... find vectors by inverse iteration ..........
      norm = abs(d(p))
      ip = p + 1
c
      do 500 i = ip, q
  500 norm = max(norm, abs(d(i)) + abs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
      eps2 = 1.0e-3 * norm
      uk = sqrt(real(q-p+5))
      eps3 = uk * machep * norm
      eps4 = uk * eps3
      uk = eps4 / sqrt(uk)
      group = 0
      s = p
c
      do 920 k = m1, m2
         r = r + 1
         its = 1
         w(r) = rv5(k)
         x1 = rv5(k)
c     .......... look for close or coincident roots ..........
         if (k .eq. m1) go to 520
         if (x1 - x0 .ge. eps2) group = -1
         group = group + 1
         if (x1 .le. x0) x1 = x0 + eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0e0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (abs(e(i)) .lt. abs(u)) go to 540
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0e0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0e0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0e0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0e0
         rv3(q) = 0.0e0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
c
         do 680 jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0e0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0e0
c
         do 720 i = p, q
  720    norm = norm + abs(rv6(i))
c
         if (norm .ge. 1.0e0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 960
         if (norm .ne. 0.0e0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0e0
c
         do 860 i = p, q
  860    u = u + rv6(i)**2
c
         xu = 1.0e0 / sqrt(u)
c
         do 880 i = 1, n
  880    z(i,r) = 0.0e0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- non-converged eigenvector ..........
  960 ierr = 4 * n + r
      go to 1001
c     .......... set error -- underestimate of number of
c                eigenvalues in interval ..........
  980 ierr = 3 * n + 1
 1001 lb = t1
      ub = t2
      return
      end
