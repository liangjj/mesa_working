*deck minfit
      subroutine minfit (nm, m, n, a, w, ip, b, ierr, rv1)
c***begin prologue  minfit
c***purpose  compute the singular value decomposition of a rectangular
c            matrix and solve the related linear least squares problem.
c***library   slatec (eispack)
c***category  d9
c***type      single precision (minfit-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure minfit,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines, towards the solution of the linear
c                                                        t
c     system ax=b, the singular value decomposition a=usv  of a real
c                                         t
c     m by n rectangular matrix, forming u b rather than u.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a and b, as declared in the calling
c          program dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.  nm is an integer
c          variable.
c
c        m is the number of rows of a and b.  m is an integer variable.
c
c        n is the number of columns of a and the order of v.  n is an
c          integer variable.
c
c        a contains the rectangular coefficient matrix of the system.
c          a is a two-dimensional real array, dimensioned a(nm,n).
c
c        ip is the number of columns of b.  ip can be zero.
c
c        b contains the constant column matrix of the system if ip is
c          not zero.  otherwise, b is not referenced.  b is a two-
c          dimensional real array, dimensioned b(nm,ip).
c
c     on output
c
c        a has been overwritten by the matrix v (orthogonal) of the
c          decomposition in its first n rows and columns.  if an
c          error exit is made, the columns of v corresponding to
c          indices of correct singular values should be correct.
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1, ierr+2, ..., n.  w is a one-dimensional
c          real array, dimensioned w(n).
c
c                                   t
c        b has been overwritten by u b.  if an error exit is made,
c                       t
c          the rows of u b corresponding to indices of correct singular
c          values should be correct.
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c                     the singular values should be correct for
c                     indices ierr+1, ierr+2, ..., n.
c
c        rv1 is a one-dimensional real array used for temporary storage,
c          dimensioned rv1(n).
c
c     calls pythag(a,b) for sqrt(a**2 + b**2).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  pythag
c***revision history  (yymmdd)
c   760101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  minfit
c
      integer i,j,k,l,m,n,ii,ip,i1,kk,k1,ll,l1,m1,nm,its,ierr
      real a(nm,*),w(*),b(nm,ip),rv1(*)
      real c,f,g,h,s,x,y,z,scale,s1
      real pythag
c
c***first executable statement  minfit
      ierr = 0
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0e0
      scale = 0.0e0
      s1 = 0.0e0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0e0
         s = 0.0e0
         scale = 0.0e0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + abs(a(k,i))
c
         if (scale .eq. 0.0e0) go to 210
c
         do 130 k = i, m
            a(k,i) = a(k,i) / scale
            s = s + a(k,i)**2
  130    continue
c
         f = a(i,i)
         g = -sign(sqrt(s),f)
         h = f * g - s
         a(i,i) = f - g
         if (i .eq. n) go to 160
c
         do 150 j = l, n
            s = 0.0e0
c
            do 140 k = i, m
  140       s = s + a(k,i) * a(k,j)
c
            f = s / h
c
            do 150 k = i, m
               a(k,j) = a(k,j) + f * a(k,i)
  150    continue
c
  160    if (ip .eq. 0) go to 190
c
         do 180 j = 1, ip
            s = 0.0e0
c
            do 170 k = i, m
  170       s = s + a(k,i) * b(k,j)
c
            f = s / h
c
            do 180 k = i, m
               b(k,j) = b(k,j) + f * a(k,i)
  180    continue
c
  190    do 200 k = i, m
  200    a(k,i) = scale * a(k,i)
c
  210    w(i) = scale * g
         g = 0.0e0
         s = 0.0e0
         scale = 0.0e0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + abs(a(i,k))
c
         if (scale .eq. 0.0e0) go to 290
c
         do 230 k = l, n
            a(i,k) = a(i,k) / scale
            s = s + a(i,k)**2
  230    continue
c
         f = a(i,l)
         g = -sign(sqrt(s),f)
         h = f * g - s
         a(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = a(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0e0
c
            do 250 k = l, n
  250       s = s + a(j,k) * a(i,k)
c
            do 260 k = l, n
               a(j,k) = a(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    a(i,k) = scale * a(i,k)
c
  290    s1 = max(s1,abs(w(i))+abs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations.
c                for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0e0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    a(j,i) = (a(i,j) / a(i,l)) / g
c
         do 350 j = l, n
            s = 0.0e0
c
            do 340 k = l, n
  340       s = s + a(i,k) * a(k,j)
c
            do 350 k = l, n
               a(k,j) = a(k,j) + s * a(k,i)
  350    continue
c
  360    do 380 j = l, n
            a(i,j) = 0.0e0
            a(j,i) = 0.0e0
  380    continue
c
  390    a(i,i) = 1.0e0
         g = rv1(i)
         l = i
  400 continue
c
      if (m .ge. n .or. ip .eq. 0) go to 510
      m1 = m + 1
c
      do 500 i = m1, n
c
         do 500 j = 1, ip
            b(i,j) = 0.0e0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
  510 continue
c     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            if (s1 + abs(rv1(l)) .eq. s1) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            if (s1 + abs(w(l1)) .eq. s1) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0e0
         s = 1.0e0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            if (s1 + abs(f) .eq. s1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (ip .eq. 0) go to 560
c
            do 550 j = 1, ip
               y = b(l1,j)
               z = b(i,j)
               b(l1,j) = y * c + z * s
               b(i,j) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5e0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0e0)
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
c     .......... next qr transformation ..........
         c = 1.0e0
         s = 1.0e0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
c
            do 570 j = 1, n
               x = a(j,i1)
               z = a(j,i)
               a(j,i1) = x * c + z * s
               a(j,i) = -x * s + z * c
  570       continue
c
            z = pythag(f,h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0e0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (ip .eq. 0) go to 600
c
            do 590 j = 1, ip
               y = b(i1,j)
               z = b(i,j)
               b(i1,j) = y * c + z * s
               b(i,j) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0e0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0e0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
c
         do 690 j = 1, n
  690    a(j,k) = -a(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
