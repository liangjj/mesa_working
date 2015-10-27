*deck dlsi
      subroutine dlsi (w, mdw, ma, mg, n, prgopt, x, rnorm, mode, ws,
     +   ip)
c***begin prologue  dlsi
c***subsidiary
c***purpose  subsidiary to dlsei
c***library   slatec
c***type      double precision (lsi-s, dlsi-d)
c***author  hanson, r. j., (snla)
c***description
c
c     this is a companion subprogram to dlsei.  the documentation for
c     dlsei has complete usage instructions.
c
c     solve..
c              ax = b,  a  ma by n  (least squares equations)
c     subject to..
c
c              gx.ge.h, g  mg by n  (inequality constraints)
c
c     input..
c
c      w(*,*) contains  (a b) in rows 1,...,ma+mg, cols 1,...,n+1.
c                       (g h)
c
c     mdw,ma,mg,n
c              contain (resp) var. dimension of w(*,*),
c              and matrix dimensions.
c
c     prgopt(*),
c              program option vector.
c
c     output..
c
c      x(*),rnorm
c
c              solution vector(unless mode=2), length of ax-b.
c
c      mode
c              =0   inequality constraints are compatible.
c              =2   inequality constraints contradictory.
c
c      ws(*),
c              working storage of dimension k+n+(mg+2)*(n+7),
c              where k=max(ma+mg,n).
c      ip(mg+2*n+1)
c              integer working storage
c
c***routines called  d1mach, dasum, daxpy, dcopy, ddot, dh12, dhfti,
c                    dlpdp, dscal, dswap
c***revision history  (yymmdd)
c   790701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890618  completely restructured and extensively revised (wrb & rwc)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900604  dp version created from sp version.  (rwc)
c   920422  changed call to dhfti to include variable ma.  (wrb)
c***end prologue  dlsi
      integer ip(*), ma, mdw, mg, mode, n
      double precision prgopt(*), rnorm, w(mdw,*), ws(*), x(*)
c
      external d1mach, dasum, daxpy, dcopy, ddot, dh12, dhfti, dlpdp,
     *   dscal, dswap
      double precision d1mach, dasum, ddot
c
      double precision anorm, drelpr, fac, gam, rb, tau, tol, xnorm
      integer i, j, k, key, krank, krm1, krp1, l, last, link, m, map1,
     *   mdlpdp, minman, n1, n2, n3, next, np1
      logical cov, first, sclcov
c
      save drelpr, first
      data first /.true./
c
c***first executable statement  dlsi
c
c     set the nominal tolerance used in the code.
c
      if (first) drelpr = d1mach(4)
      first = .false.
      tol = sqrt(drelpr)
c
      mode = 0
      rnorm = 0.d0
      m = ma + mg
      np1 = n + 1
      krank = 0
      if (n.le.0 .or. m.le.0) go to 370
c
c     to process option vector.
c
      cov = .false.
      sclcov = .true.
      last = 1
      link = prgopt(1)
c
  100 if (link.gt.1) then
         key = prgopt(last+1)
         if (key.eq.1) cov = prgopt(last+2) .ne. 0.d0
         if (key.eq.10) sclcov = prgopt(last+2) .eq. 0.d0
         if (key.eq.5) tol = max(drelpr,prgopt(last+2))
         next = prgopt(link)
         last = link
         link = next
         go to 100
      endif
c
c     compute matrix norm of least squares equations.
c
      anorm = 0.d0
      do 110 j = 1,n
         anorm = max(anorm,dasum(ma,w(1,j),1))
  110 continue
c
c     set tolerance for dhfti( ) rank test.
c
      tau = tol*anorm
c
c     compute householder orthogonal decomposition of matrix.
c
      call dcopy (n, 0.d0, 0, ws, 1)
      call dcopy (ma, w(1, np1), 1, ws, 1)
      k = max(m,n)
      minman = min(ma,n)
      n1 = k + 1
      n2 = n1 + n
      call dhfti (w, mdw, ma, n, ws, ma, 1, tau, krank, rnorm, ws(n2),
     +           ws(n1), ip)
      fac = 1.d0
      gam = ma - krank
      if (krank.lt.ma .and. sclcov) fac = rnorm**2/gam
c
c     reduce to dlpdp and solve.
c
      map1 = ma + 1
c
c     compute inequality rt-hand side for dlpdp.
c
      if (ma.lt.m) then
         if (minman.gt.0) then
            do 120 i = map1,m
               w(i,np1) = w(i,np1) - ddot(n,w(i,1),mdw,ws,1)
  120       continue
c
c           apply permutations to col. of inequality constraint matrix.
c
            do 130 i = 1,minman
               call dswap (mg, w(map1,i), 1, w(map1,ip(i)), 1)
  130       continue
c
c           apply householder transformations to constraint matrix.
c
            if (krank.gt.0 .and. krank.lt.n) then
               do 140 i = krank,1,-1
                  call dh12 (2, i, krank+1, n, w(i,1), mdw, ws(n1+i-1),
     +                      w(map1,1), mdw, 1, mg)
  140          continue
            endif
c
c           compute permuted inequality constraint matrix times r-inv.
c
            do 160 i = map1,m
               do 150 j = 1,krank
                  w(i,j) = (w(i,j)-ddot(j-1,w(1,j),1,w(i,1),mdw))/w(j,j)
  150          continue
  160       continue
         endif
c
c        solve the reduced problem with dlpdp algorithm,
c        the least projected distance problem.
c
         call dlpdp(w(map1,1), mdw, mg, krank, n-krank, prgopt, x,
     +             xnorm, mdlpdp, ws(n2), ip(n+1))
c
c        compute solution in original coordinates.
c
         if (mdlpdp.eq.1) then
            do 170 i = krank,1,-1
               x(i) = (x(i)-ddot(krank-i,w(i,i+1),mdw,x(i+1),1))/w(i,i)
  170       continue
c
c           apply householder transformation to solution vector.
c
            if (krank.lt.n) then
               do 180 i = 1,krank
                  call dh12 (2, i, krank+1, n, w(i,1), mdw, ws(n1+i-1),
     +                      x, 1, 1, 1)
  180          continue
            endif
c
c           repermute variables to their input order.
c
            if (minman.gt.0) then
               do 190 i = minman,1,-1
                  call dswap (1, x(i), 1, x(ip(i)), 1)
  190          continue
c
c              variables are now in original coordinates.
c              add solution of unconstrained problem.
c
               do 200 i = 1,n
                  x(i) = x(i) + ws(i)
  200          continue
c
c              compute the residual vector norm.
c
               rnorm = sqrt(rnorm**2+xnorm**2)
            endif
         else
            mode = 2
         endif
      else
         call dcopy (n, ws, 1, x, 1)
      endif
c
c     compute covariance matrix based on the orthogonal decomposition
c     from dhfti( ).
c
      if (.not.cov .or. krank.le.0) go to 370
      krm1 = krank - 1
      krp1 = krank + 1
c
c     copy diagonal terms to working array.
c
      call dcopy (krank, w, mdw+1, ws(n2), 1)
c
c     reciprocate diagonal terms.
c
      do 210 j = 1,krank
         w(j,j) = 1.d0/w(j,j)
  210 continue
c
c     invert the upper triangular qr factor on itself.
c
      if (krank.gt.1) then
         do 230 i = 1,krm1
            do 220 j = i+1,krank
               w(i,j) = -ddot(j-i,w(i,i),mdw,w(i,j),1)*w(j,j)
  220       continue
  230    continue
      endif
c
c     compute the inverted factor times its transpose.
c
      do 250 i = 1,krank
         do 240 j = i,krank
            w(i,j) = ddot(krank+1-j,w(i,j),mdw,w(j,j),mdw)
  240    continue
  250 continue
c
c     zero out lower trapezoidal part.
c     copy upper triangular to lower triangular part.
c
      if (krank.lt.n) then
         do 260 j = 1,krank
            call dcopy (j, w(1,j), 1, w(j,1), mdw)
  260    continue
c
         do 270 i = krp1,n
            call dcopy (i, 0.d0, 0, w(i,1), mdw)
  270    continue
c
c        apply right side transformations to lower triangle.
c
         n3 = n2 + krp1
         do 330 i = 1,krank
            l = n1 + i
            k = n2 + i
            rb = ws(l-1)*ws(k-1)
c
c           if rb.ge.0.d0, transformation can be regarded as zero.
c
            if (rb.lt.0.d0) then
               rb = 1.d0/rb
c
c              store unscaled rank one householder update in work array.
c
               call dcopy (n, 0.d0, 0, ws(n3), 1)
               l = n1 + i
               k = n3 + i
               ws(k-1) = ws(l-1)
c
               do 280 j = krp1,n
                  ws(n3+j-1) = w(i,j)
  280          continue
c
               do 290 j = 1,n
                  ws(j) = rb*(ddot(j-i,w(j,i),mdw,ws(n3+i-1),1)+
     +                    ddot(n-j+1,w(j,j),1,ws(n3+j-1),1))
  290          continue
c
               l = n3 + i
               gam = 0.5d0*rb*ddot(n-i+1,ws(l-1),1,ws(i),1)
               call daxpy (n-i+1, gam, ws(l-1), 1, ws(i), 1)
               do 320 j = i,n
                  do 300 l = 1,i-1
                     w(j,l) = w(j,l) + ws(n3+j-1)*ws(l)
  300             continue
c
                  do 310 l = i,j
                     w(j,l) = w(j,l) + ws(j)*ws(n3+l-1)+ws(l)*ws(n3+j-1)
  310             continue
  320          continue
            endif
  330    continue
c
c        copy lower triangle to upper triangle to symmetrize the
c        covariance matrix.
c
         do 340 i = 1,n
            call dcopy (i, w(i,1), mdw, w(1,i), 1)
  340    continue
      endif
c
c     repermute rows and columns.
c
      do 350 i = minman,1,-1
         k = ip(i)
         if (i.ne.k) then
            call dswap (1, w(i,i), 1, w(k,k), 1)
            call dswap (i-1, w(1,i), 1, w(1,k), 1)
            call dswap (k-i-1, w(i,i+1), mdw, w(i+1,k), 1)
            call dswap (n-k, w(i, k+1), mdw, w(k, k+1), mdw)
         endif
  350 continue
c
c     put in normalized residual sum of squares scale factor
c     and symmetrize the resulting covariance matrix.
c
      do 360 j = 1,n
         call dscal (j, fac, w(1,j), 1)
         call dcopy (j, w(1,j), 1, w(j,1), mdw)
  360 continue
c
  370 ip(1) = krank
      ip(2) = n + max(m,n) + (mg+2)*(n+7)
      return
      end
