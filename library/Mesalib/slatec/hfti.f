*deck hfti
      subroutine hfti (a, mda, m, n, b, mdb, nb, tau, krank, rnorm, h,
     +   g, ip)
c***begin prologue  hfti
c***purpose  solve a linear least squares problems by performing a qr
c            factorization of the matrix using householder
c            transformations.
c***library   slatec
c***category  d9
c***type      single precision (hfti-s, dhfti-d)
c***keywords  curve fitting, linear least squares, qr factorization
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c***description
c
c     dimension a(mda,n),(b(mdb,nb) or b(m)),rnorm(nb),h(n),g(n),ip(n)
c
c     this subroutine solves a linear least squares problem or a set of
c     linear least squares problems having the same matrix but different
c     right-side vectors.  the problem data consists of an m by n matrix
c     a, an m by nb matrix b, and an absolute tolerance parameter tau
c     whose usage is described below.  the nb column vectors of b
c     represent right-side vectors for nb distinct linear least squares
c     problems.
c
c     this set of problems can also be written as the matrix least
c     squares problem
c
c                       ax = b,
c
c     where x is the n by nb solution matrix.
c
c     note that if b is the m by m identity matrix, then x will be the
c     pseudo-inverse of a.
c
c     this subroutine first transforms the augmented matrix (a b) to a
c     matrix (r c) using premultiplying householder transformations with
c     column interchanges.  all subdiagonal elements in the matrix r are
c     zero and its diagonal elements satisfy
c
c                       abs(r(i,i)).ge.abs(r(i+1,i+1)),
c
c                       i = 1,...,l-1, where
c
c                       l = min(m,n).
c
c     the subroutine will compute an integer, krank, equal to the number
c     of diagonal terms of r that exceed tau in magnitude. then a
c     solution of minimum euclidean length is computed using the first
c     krank rows of (r c).
c
c     to be specific we suggest that the user consider an easily
c     computable matrix norm, such as, the maximum of all column sums of
c     magnitudes.
c
c     now if the relative uncertainty of b is eps, (norm of uncertainty/
c     norm of b), it is suggested that tau be set approximately equal to
c     eps*(norm of a).
c
c     the user must dimension all arrays appearing in the call list..
c     a(mda,n),(b(mdb,nb) or b(m)),rnorm(nb),h(n),g(n),ip(n).  this
c     permits the solution of a range of problems in the same array
c     space.
c
c     the entire set of parameters for hfti are
c
c     input..
c
c     a(*,*),mda,m,n    the array a(*,*) initially contains the m by n
c                       matrix a of the least squares problem ax = b.
c                       the first dimensioning parameter of the array
c                       a(*,*) is mda, which must satisfy mda.ge.m
c                       either m.ge.n or m.lt.n is permitted.  there
c                       is no restriction on the rank of a.  the
c                       condition mda.lt.m is considered an error.
c
c     b(*),mdb,nb       if nb = 0 the subroutine will perform the
c                       orthogonal decomposition but will make no
c                       references to the array b(*).  if nb.gt.0
c                       the array b(*) must initially contain the m by
c                       nb matrix b of the least squares problem ax =
c                       b.  if nb.ge.2 the array b(*) must be doubly
c                       subscripted with first dimensioning parameter
c                       mdb.ge.max(m,n).  if nb = 1 the array b(*) may
c                       be either doubly or singly subscripted.  in
c                       the latter case the value of mdb is arbitrary
c                       but it should be set to some valid integer
c                       value such as mdb = m.
c
c                       the condition of nb.gt.1.and.mdb.lt. max(m,n)
c                       is considered an error.
c
c     tau               absolute tolerance parameter provided by user
c                       for pseudorank determination.
c
c     h(*),g(*),ip(*)   arrays of working space used by hfti.
c
c     output..
c
c     a(*,*)            the contents of the array a(*,*) will be
c                       modified by the subroutine. these contents
c                       are not generally required by the user.
c
c     b(*)              on return the array b(*) will contain the n by
c                       nb solution matrix x.
c
c     krank             set by the subroutine to indicate the
c                       pseudorank of a.
c
c     rnorm(*)          on return, rnorm(j) will contain the euclidean
c                       norm of the residual vector for the problem
c                       defined by the j-th column vector of the array
c                       b(*,*) for j = 1,...,nb.
c
c     h(*),g(*)         on return these arrays respectively contain
c                       elements of the pre- and post-multiplying
c                       householder transformations used to compute
c                       the minimum euclidean length solution.
c
c     ip(*)             array in which the subroutine records indices
c                       describing the permutation of column vectors.
c                       the contents of arrays h(*),g(*) and ip(*)
c                       are not generally required by the user.
c
c***references  c. l. lawson and r. j. hanson, solving least squares
c                 problems, prentice-hall, inc., 1974, chapter 14.
c***routines called  h12, r1mach, xermsg
c***revision history  (yymmdd)
c   790101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   901005  replace usage of diff with usage of r1mach.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hfti
      dimension a(mda,*),b(mdb,*),h(*),g(*),rnorm(*)
      integer ip(*)
      double precision sm,dzero
      save releps
      data releps /0.e0/
c***first executable statement  hfti
      if (releps.eq.0) releps = r1mach(4)
      szero=0.
      dzero=0.d0
      factor=0.001
c
      k=0
      ldiag=min(m,n)
      if (ldiag.le.0) go to 270
      if (.not.mda.lt.m) go to 5
      nerr=1
      iopt=2
      call xermsg ('slatec', 'hfti', 'mda.lt.m, probable error.',
     +   nerr, iopt)
      return
    5 continue
c
      if (.not.(nb.gt.1.and.max(m,n).gt.mdb)) go to 6
      nerr=2
      iopt=2
      call xermsg ('slatec', 'hfti',
     +   'mdb.lt.max(m,n).and.nb.gt.1. probable error.', nerr, iopt)
      return
    6 continue
c
          do 80 j=1,ldiag
          if (j.eq.1) go to 20
c
c     update squared column lengths and find lmax
c    ..
          lmax=j
              do 10 l=j,n
              h(l)=h(l)-a(j-1,l)**2
              if (h(l).gt.h(lmax)) lmax=l
   10         continue
          if (factor*h(lmax) .gt. hmax*releps) go to 50
c
c     compute squared column lengths and find lmax
c    ..
   20     lmax=j
              do 40 l=j,n
              h(l)=0.
                  do 30 i=j,m
   30             h(l)=h(l)+a(i,l)**2
              if (h(l).gt.h(lmax)) lmax=l
   40         continue
          hmax=h(lmax)
c    ..
c     lmax has been determined
c
c     do column interchanges if needed.
c    ..
   50     continue
          ip(j)=lmax
          if (ip(j).eq.j) go to 70
              do 60 i=1,m
              tmp=a(i,j)
              a(i,j)=a(i,lmax)
   60         a(i,lmax)=tmp
          h(lmax)=h(j)
c
c     compute the j-th transformation and apply it to a and b.
c    ..
   70     call h12 (1,j,j+1,m,a(1,j),1,h(j),a(1,j+1),1,mda,n-j)
   80     call h12 (2,j,j+1,m,a(1,j),1,h(j),b,1,mdb,nb)
c
c     determine the pseudorank, k, using the tolerance, tau.
c    ..
          do 90 j=1,ldiag
          if (abs(a(j,j)).le.tau) go to 100
   90     continue
      k=ldiag
      go to 110
  100 k=j-1
  110 kp1=k+1
c
c     compute the norms of the residual vectors.
c
      if (nb.le.0) go to 140
          do 130 jb=1,nb
          tmp=szero
          if (kp1.gt.m) go to 130
              do 120 i=kp1,m
  120         tmp=tmp+b(i,jb)**2
  130     rnorm(jb)=sqrt(tmp)
  140 continue
c                                           special for pseudorank = 0
      if (k.gt.0) go to 160
      if (nb.le.0) go to 270
          do 150 jb=1,nb
              do 150 i=1,n
  150         b(i,jb)=szero
      go to 270
c
c     if the pseudorank is less than n compute householder
c     decomposition of first k rows.
c    ..
  160 if (k.eq.n) go to 180
          do 170 ii=1,k
          i=kp1-ii
  170     call h12 (1,i,kp1,n,a(i,1),mda,g(i),a,mda,1,i-1)
  180 continue
c
c
      if (nb.le.0) go to 270
          do 260 jb=1,nb
c
c     solve the k by k triangular system.
c    ..
              do 210 l=1,k
              sm=dzero
              i=kp1-l
              if (i.eq.k) go to 200
              ip1=i+1
                  do 190 j=ip1,k
  190             sm=sm+a(i,j)*dble(b(j,jb))
  200         sm1=sm
  210         b(i,jb)=(b(i,jb)-sm1)/a(i,i)
c
c     complete computation of solution vector.
c    ..
          if (k.eq.n) go to 240
              do 220 j=kp1,n
  220         b(j,jb)=szero
              do 230 i=1,k
  230         call h12 (2,i,kp1,n,a(i,1),mda,g(i),b(1,jb),1,mdb,1)
c
c      re-order the solution vector to compensate for the
c      column interchanges.
c    ..
  240         do 250 jj=1,ldiag
              j=ldiag+1-jj
              if (ip(j).eq.j) go to 250
              l=ip(j)
              tmp=b(l,jb)
              b(l,jb)=b(j,jb)
              b(j,jb)=tmp
  250         continue
  260     continue
c    ..
c     the solution vectors, x, are now
c     in the first  n  rows of the array b(,).
c
  270 krank=k
      return
      end
