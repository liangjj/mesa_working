*deck dlpdp
      subroutine dlpdp (a, mda, m, n1, n2, prgopt, x, wnorm, mode, ws,
     +   is)
c***begin prologue  dlpdp
c***subsidiary
c***purpose  subsidiary to dlsei
c***library   slatec
c***type      double precision (lpdp-s, dlpdp-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c  **** double precision version of lpdp ****
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
c     called by subprogram dlsi( ).
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
c***see also  dlsei
c***routines called  dcopy, ddot, dnrm2, dscal, dwnnls
c***revision history  (yymmdd)
c   790701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated the author section.  (wrb)
c***end prologue  dlpdp
c
      integer i, is(*), iw, ix, j, l, m, mda, mode, modew, n, n1, n2,
     *     np1
      double precision a(mda,*), ddot, dnrm2, fac, one,
     *     prgopt(*), rnorm, sc, wnorm, ws(*), x(*), ynorm, zero
      save zero, one, fac
      data zero,one /0.0d0,1.0d0/, fac /0.1d0/
c***first executable statement  dlpdp
      n = n1 + n2
      mode = 1
      if (m .gt. 0) go to 20
         if (n .le. 0) go to 10
            x(1) = zero
            call dcopy(n,x,0,x,1)
   10    continue
         wnorm = zero
      go to 200
   20 continue
c        begin block permitting ...exits to 190
            np1 = n + 1
c
c           scale nonzero rows of inequality matrix to have length one.
            do 40 i = 1, m
               sc = dnrm2(n,a(i,1),mda)
               if (sc .eq. zero) go to 30
                  sc = one/sc
                  call dscal(np1,sc,a(i,1),mda)
   30          continue
   40       continue
c
c           scale rt.-side vector to have length one (or zero).
            ynorm = dnrm2(m,a(1,np1),1)
            if (ynorm .eq. zero) go to 50
               sc = one/ynorm
               call dscal(m,sc,a(1,np1),1)
   50       continue
c
c           scale cols of matrix h.
            j = n1 + 1
   60       if (j .gt. n) go to 70
               sc = dnrm2(m,a(1,j),1)
               if (sc .ne. zero) sc = one/sc
               call dscal(m,sc,a(1,j),1)
               x(j) = sc
               j = j + 1
            go to 60
   70       continue
            if (n1 .le. 0) go to 130
c
c              copy transpose of (h g y) to work array ws(*).
               iw = 0
               do 80 i = 1, m
c
c                 move col of transpose of h into work array.
                  call dcopy(n2,a(i,n1+1),mda,ws(iw+1),1)
                  iw = iw + n2
c
c                 move col of transpose of g into work array.
                  call dcopy(n1,a(i,1),mda,ws(iw+1),1)
                  iw = iw + n1
c
c                 move component of vector y into work array.
                  ws(iw+1) = a(i,np1)
                  iw = iw + 1
   80          continue
               ws(iw+1) = zero
               call dcopy(n,ws(iw+1),0,ws(iw+1),1)
               iw = iw + n
               ws(iw+1) = one
               iw = iw + 1
c
c              solve eu=f subject to (transpose of h)u=0, u.ge.0.  the
c              matrix e = transpose of (g y), and the (n+1)-vector
c              f = transpose of (0,...,0,1).
               ix = iw + 1
               iw = iw + m
c
c              do not check lengths of work arrays in this usage of
c              dwnnls( ).
               is(1) = 0
               is(2) = 0
               call dwnnls(ws,np1,n2,np1-n2,m,0,prgopt,ws(ix),rnorm,
     *                     modew,is,ws(iw+1))
c
c              compute the components of the soln denoted above by w.
               sc = one - ddot(m,a(1,np1),1,ws(ix),1)
               if (one + fac*abs(sc) .eq. one .or. rnorm .le. zero)
     *            go to 110
                  sc = one/sc
                  do 90 j = 1, n1
                     x(j) = sc*ddot(m,a(1,j),1,ws(ix),1)
   90             continue
c
c                 compute the vector q=y-gw.  overwrite y with this
c                 vector.
                  do 100 i = 1, m
                     a(i,np1) = a(i,np1) - ddot(n1,a(i,1),mda,x,1)
  100             continue
               go to 120
  110          continue
                  mode = 2
c        .........exit
                  go to 190
  120          continue
  130       continue
            if (n2 .le. 0) go to 180
c
c              copy transpose of (h q) to work array ws(*).
               iw = 0
               do 140 i = 1, m
                  call dcopy(n2,a(i,n1+1),mda,ws(iw+1),1)
                  iw = iw + n2
                  ws(iw+1) = a(i,np1)
                  iw = iw + 1
  140          continue
               ws(iw+1) = zero
               call dcopy(n2,ws(iw+1),0,ws(iw+1),1)
               iw = iw + n2
               ws(iw+1) = one
               iw = iw + 1
               ix = iw + 1
               iw = iw + m
c
c              solve rv=s subject to v.ge.0.  the matrix r =(transpose
c              of (h q)), where q=y-gw.  the (n2+1)-vector s =(transpose
c              of (0,...,0,1)).
c
c              do not check lengths of work arrays in this usage of
c              dwnnls( ).
               is(1) = 0
               is(2) = 0
               call dwnnls(ws,n2+1,0,n2+1,m,0,prgopt,ws(ix),rnorm,modew,
     *                     is,ws(iw+1))
c
c              compute the components of the soln denoted above by z.
               sc = one - ddot(m,a(1,np1),1,ws(ix),1)
               if (one + fac*abs(sc) .eq. one .or. rnorm .le. zero)
     *            go to 160
                  sc = one/sc
                  do 150 j = 1, n2
                     l = n1 + j
                     x(l) = sc*ddot(m,a(1,l),1,ws(ix),1)*x(l)
  150             continue
               go to 170
  160          continue
                  mode = 2
c        .........exit
                  go to 190
  170          continue
  180       continue
c
c           account for scaling of rt.-side vector in solution.
            call dscal(n,ynorm,x,1)
            wnorm = dnrm2(n1,x,1)
  190    continue
  200 continue
      return
      end
