*deck qzvec
      subroutine qzvec (nm, n, a, b, alfr, alfi, beta, z)
c***begin prologue  qzvec
c***purpose  the optional fourth step of the qz algorithm for
c            generalized eigenproblems.  accepts a matrix in
c            quasi-triangular form and another in upper triangular
c            and computes the eigenvectors of the triangular problem
c            and transforms them back to the original coordinates
c            usually preceded by qzhes, qzit, and qzval.
c***library   slatec (eispack)
c***category  d4c3
c***type      single precision (qzvec-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is the optional fourth step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them in
c     quasi-triangular form (in which each 2-by-2 block corresponds to
c     a pair of complex eigenvalues) and the other in upper triangular
c     form.  it computes the eigenvectors of the triangular problem and
c     transforms the results back to the original coordinate system.
c     it is usually preceded by  qzhes,  qzit, and  qzval.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a, b, and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrices a and b.  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        a contains a real upper quasi-triangular matrix.  a is a two-
c          dimensional real array, dimensioned a(nm,n).
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.  b is a two-dimensional real
c          array, dimensioned b(nm,n).
c
c        alfr, alfi, and beta are one-dimensional real arrays with
c          components whose ratios ((alfr+i*alfi)/beta) are the
c          generalized eigenvalues.  they are usually obtained from
c          qzval.  they are dimensioned alfr(n), alfi(n), and beta(n).
c
c        z contains the transformation matrix produced in the reductions
c          by  qzhes,  qzit, and  qzval,  if performed.  if the
c          eigenvectors of the triangular problem are desired, z must
c          contain the identity matrix.  z is a two-dimensional real
c          array, dimensioned z(nm,n).
c
c     on output
c
c        a is unaltered.  its subdiagonal elements provide information
c           about the storage of the complex eigenvectors.
c
c        b has been destroyed.
c
c        alfr, alfi, and beta are unaltered.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if alfi(j) .eq. 0.0, the j-th eigenvalue is real and
c            the j-th column of z contains its eigenvector.
c          if alfi(j) .ne. 0.0, the j-th eigenvalue is complex.
c            if alfi(j) .gt. 0.0, the eigenvalue is the first of
c              a complex pair and the j-th and (j+1)-th columns
c              of z contain its eigenvector.
c            if alfi(j) .lt. 0.0, the eigenvalue is the second of
c              a complex pair and the (j-1)-th and j-th columns
c              of z contain the conjugate of its eigenvector.
c          each eigenvector is normalized so that the modulus
c          of its largest component is 1.0 .
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  (none)
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  qzvec
c
      integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2
      real a(nm,*),b(nm,*),alfr(*),alfi(*),beta(*),z(nm,*)
      real d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2
      real w1,x1,zz,z1,alfm,almi,almr,betm,epsb
c
c***first executable statement  qzvec
      epsb = b(n,1)
      isw = 1
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 795
         if (alfi(en) .ne. 0.0e0) go to 710
c     .......... real vector ..........
         m = en
         b(en,en) = 1.0e0
         if (na .eq. 0) go to 800
         alfm = alfr(m)
         betm = beta(m)
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = betm * a(i,i) - alfm * b(i,i)
            r = 0.0e0
c
            do 610 j = m, en
  610       r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
c
            if (i .eq. 1 .or. isw .eq. 2) go to 630
            if (betm * a(i,i-1) .eq. 0.0e0) go to 630
            zz = w
            s = r
            go to 690
  630       m = i
            if (isw .eq. 2) go to 640
c     .......... real 1-by-1 block ..........
            t = w
            if (w .eq. 0.0e0) t = epsb
            b(i,en) = -r / t
            go to 700
c     .......... real 2-by-2 block ..........
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)
            y = betm * a(i+1,i)
            q = w * zz - x * y
            t = (x * s - zz * r) / q
            b(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            b(i+1,en) = (-r - w * t) / x
            go to 690
  650       b(i+1,en) = (-s - y * t) / zz
  690       isw = 3 - isw
  700    continue
c     .......... end real vector ..........
         go to 800
c     .......... complex vector ..........
  710    m = na
         almr = alfr(m)
         almi = alfi(m)
         betm = beta(m)
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         y = betm * a(en,na)
         b(na,na) = -almi * b(en,en) / y
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
         b(en,na) = 0.0e0
         b(en,en) = 1.0e0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 795
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 790 ii = 1, enm2
            i = na - ii
            w = betm * a(i,i) - almr * b(i,i)
            w1 = -almi * b(i,i)
            ra = 0.0e0
            sa = 0.0e0
c
            do 760 j = m, en
               x = betm * a(i,j) - almr * b(i,j)
               x1 = -almi * b(i,j)
               ra = ra + x * b(j,na) - x1 * b(j,en)
               sa = sa + x * b(j,en) + x1 * b(j,na)
  760       continue
c
            if (i .eq. 1 .or. isw .eq. 2) go to 770
            if (betm * a(i,i-1) .eq. 0.0e0) go to 770
            zz = w
            z1 = w1
            r = ra
            s = sa
            isw = 2
            go to 790
  770       m = i
            if (isw .eq. 2) go to 780
c     .......... complex 1-by-1 block ..........
            tr = -ra
            ti = -sa
  773       dr = w
            di = w1
c     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
  775       if (abs(di) .gt. abs(dr)) go to 777
            rr = di / dr
            d = dr + di * rr
            t1 = (tr + ti * rr) / d
            t2 = (ti - tr * rr) / d
            go to (787,782), isw
  777       rr = dr / di
            d = dr * rr + di
            t1 = (tr * rr + ti) / d
            t2 = (ti * rr - tr) / d
            go to (787,782), isw
c     .......... complex 2-by-2 block ..........
  780       x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = -almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if (dr .eq. 0.0e0 .and. di .eq. 0.0e0) dr = epsb
            go to 775
  782       b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1
            if (abs(y) .gt. abs(w) + abs(w1)) go to 785
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
            go to 773
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
  787       b(i,na) = t1
            b(i,en) = t2
  790    continue
c     .......... end complex vector ..........
  795    isw = 3 - isw
  800 continue
c     .......... end back substitution.
c                transform to original coordinate system.
c                for j=n step -1 until 1 do -- ..........
      do 880 jj = 1, n
         j = n + 1 - jj
c
         do 880 i = 1, n
            zz = 0.0e0
c
            do 860 k = 1, j
  860       zz = zz + z(i,k) * b(k,j)
c
            z(i,j) = zz
  880 continue
c     .......... normalize so that modulus of largest
c                component of each vector is 1.
c                (isw is 1 initially from before) ..........
      do 950 j = 1, n
         d = 0.0e0
         if (isw .eq. 2) go to 920
         if (alfi(j) .ne. 0.0e0) go to 945
c
         do 890 i = 1, n
            if (abs(z(i,j)) .gt. d) d = abs(z(i,j))
  890    continue
c
         do 900 i = 1, n
  900    z(i,j) = z(i,j) / d
c
         go to 950
c
  920    do 930 i = 1, n
            r = abs(z(i,j-1)) + abs(z(i,j))
            if (r .ne. 0.0e0) r = r * sqrt((z(i,j-1)/r)**2
     1                                     +(z(i,j)/r)**2)
            if (r .gt. d) d = r
  930    continue
c
         do 940 i = 1, n
            z(i,j-1) = z(i,j-1) / d
            z(i,j) = z(i,j) / d
  940    continue
c
  945    isw = 3 - isw
  950 continue
c
      return
      end
