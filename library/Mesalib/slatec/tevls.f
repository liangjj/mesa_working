*deck tevls
      subroutine tevls (n, d, e2, ierr)
c***begin prologue  tevls
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (tevls-s)
c***author  (unknown)
c***description
c
c     this subroutine finds the eigenvalues of a symmetric tridiagonal
c     matrix by the rational ql method.
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e2 contains the subdiagonal elements of the input matrix
c           in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        e2 has been destroyed,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c***see also  blktri
c***references  c. h. reinsch, eigenvalues of a real, symmetric, tri-
c                 diagonal matrix, algorithm 464, communications of the
c                 acm 16, 11 (november 1973), pp. 689.
c***routines called  (none)
c***common blocks    cblkt
c***revision history  (yymmdd)
c   801001  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920528  description revised and references section added.  (wrb)
c***end prologue  tevls
c
      integer         i          ,j          ,l          ,m          ,
     1                n          ,ii         ,l1         ,mml        ,
     2                ierr
      real            d(*)       ,e2(*)
      real            b          ,c          ,f          ,g          ,
     1                h          ,p          ,r          ,s          ,
     2                machep
c
      common /cblkt/  npp        ,k          ,machep     ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  tevls
      ierr = 0
      if (n .eq. 1) go to 115
c
      do 101 i=2,n
         e2(i-1) = e2(i)*e2(i)
  101 continue
c
      f = 0.0
      b = 0.0
      e2(n) = 0.0
c
      do 112 l=1,n
         j = 0
         h = machep*(abs(d(l))+sqrt(e2(l)))
         if (b .gt. h) go to 102
         b = h
         c = b*b
c
c     ********** look for small squared sub-diagonal element **********
c
  102    do 103 m=l,n
            if (e2(m) .le. c) go to 104
c
c     ********** e2(n) is always zero, so there is no exit
c                through the bottom of the loop **********
c
  103    continue
c
  104    if (m .eq. l) go to 108
  105    if (j .eq. 30) go to 114
         j = j+1
c
c     ********** form shift **********
c
         l1 = l+1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1)-g)/(2.0*s)
         r = sqrt(p*p+1.0)
         d(l) = s/(p+sign(r,p))
         h = g-d(l)
c
         do 106 i=l1,n
            d(i) = d(i)-h
  106    continue
c
         f = f+h
c
c     ********** rational ql transformation **********
c
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m-l
c
c     ********** for i=m-1 step -1 until l do -- **********
c
         do 107 ii=1,mml
            i = m-ii
            p = g*h
            r = p+e2(i)
            e2(i+1) = s*r
            s = e2(i)/r
            d(i+1) = h+s*(h+d(i))
            g = d(i)-e2(i)/g
            if (g .eq. 0.0) g = b
            h = g*p/r
  107    continue
c
         e2(l) = s*g
         d(l) = h
c
c     ********** guard against underflowed h **********
c
         if (h .eq. 0.0) go to 108
         if (abs(e2(l)) .le. abs(c/h)) go to 108
         e2(l) = h*e2(l)
         if (e2(l) .ne. 0.0) go to 105
  108    p = d(l)+f
c
c     ********** order eigenvalues **********
c
         if (l .eq. 1) go to 110
c
c     ********** for i=l step -1 until 2 do -- **********
c
         do 109 ii=2,l
            i = l+2-ii
            if (p .ge. d(i-1)) go to 111
            d(i) = d(i-1)
  109    continue
c
  110    i = 1
  111    d(i) = p
  112 continue
c
      if (abs(d(n)) .ge. abs(d(1))) go to 115
      nhalf = n/2
      do 113 i=1,nhalf
         ntop = n-i
         dhold = d(i)
         d(i) = d(ntop+1)
         d(ntop+1) = dhold
  113 continue
      go to 115
c
c     ********** set error -- no convergence to an
c                eigenvalue after 30 iterations **********
c
  114 ierr = l
  115 return
c
c     ********** last card of tqlrat **********
c
      end
