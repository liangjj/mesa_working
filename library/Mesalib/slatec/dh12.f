*deck dh12
      subroutine dh12 (mode, lpivot, l1, m, u, iue, up, c, ice, icv,
     +   ncv)
c***begin prologue  dh12
c***subsidiary
c***purpose  subsidiary to dhfti, dlsei and dwnnls
c***library   slatec
c***type      double precision (h12-s, dh12-d)
c***author  (unknown)
c***description
c
c      *** double precision version of h12 ******
c
c     c.l.lawson and r.j.hanson, jet propulsion laboratory, 1973 jun 12
c     to appear in 'solving least squares problems', prentice-hall, 1974
c
c     construction and/or application of a single
c     householder transformation..     q = i + u*(u**t)/b
c
c     mode    = 1 or 2   to select algorithm  h1  or  h2 .
c     lpivot is the index of the pivot element.
c     l1,m   if l1 .le. m   the transformation will be constructed to
c            zero elements indexed from l1 through m.   if l1 gt. m
c            the subroutine does an identity transformation.
c     u(),iue,up    on entry to h1 u() contains the pivot vector.
c                   iue is the storage increment between elements.
c                                       on exit from h1 u() and up
c                   contain quantities defining the vector u of the
c                   householder transformation.   on entry to h2 u()
c                   and up should contain quantities previously computed
c                   by h1.  these will not be modified by h2.
c     c()    on entry to h1 or h2 c() contains a matrix which will be
c            regarded as a set of vectors to which the householder
c            transformation is to be applied.  on exit c() contains the
c            set of transformed vectors.
c     ice    storage increment between elements of vectors in c().
c     icv    storage increment between vectors in c().
c     ncv    number of vectors in c() to be transformed. if ncv .le. 0
c            no operations will be done on c().
c
c***see also  dhfti, dlsei, dwnnls
c***routines called  daxpy, ddot, dswap
c***revision history  (yymmdd)
c   790101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   900911  added ddot to double precision statement.  (wrb)
c***end prologue  dh12
      integer i, i2, i3, i4, ice, icv, incr, iue, j, kl1, kl2, klp,
     *     l1, l1m1, lpivot, m, mml1p2, mode, ncv
      double precision b, c, cl, clinv, one, ul1m1, sm, u, up, ddot
      dimension u(iue,*), c(*)
c     begin block permitting ...exits to 140
c***first executable statement  dh12
         one = 1.0d0
c
c     ...exit
         if (0 .ge. lpivot .or. lpivot .ge. l1 .or. l1 .gt. m) go to 140
         cl = abs(u(1,lpivot))
         if (mode .eq. 2) go to 40
c           ****** construct the transformation. ******
            do 10 j = l1, m
               cl = max(abs(u(1,j)),cl)
   10       continue
            if (cl .gt. 0.0d0) go to 20
c     .........exit
               go to 140
   20       continue
            clinv = one/cl
            sm = (u(1,lpivot)*clinv)**2
            do 30 j = l1, m
               sm = sm + (u(1,j)*clinv)**2
   30       continue
            cl = cl*sqrt(sm)
            if (u(1,lpivot) .gt. 0.0d0) cl = -cl
            up = u(1,lpivot) - cl
            u(1,lpivot) = cl
         go to 50
   40    continue
c        ****** apply the transformation  i+u*(u**t)/b  to c. ******
c
         if (cl .gt. 0.0d0) go to 50
c     ......exit
            go to 140
   50    continue
c     ...exit
         if (ncv .le. 0) go to 140
         b = up*u(1,lpivot)
c        b  must be nonpositive here.  if b = 0., return.
c
         if (b .lt. 0.0d0) go to 60
c     ......exit
            go to 140
   60    continue
         b = one/b
         mml1p2 = m - l1 + 2
         if (mml1p2 .le. 20) go to 80
            l1m1 = l1 - 1
            kl1 = 1 + (l1m1 - 1)*ice
            kl2 = kl1
            klp = 1 + (lpivot - 1)*ice
            ul1m1 = u(1,l1m1)
            u(1,l1m1) = up
            if (lpivot .ne. l1m1) call dswap(ncv,c(kl1),icv,c(klp),icv)
            do 70 j = 1, ncv
               sm = ddot(mml1p2,u(1,l1m1),iue,c(kl1),ice)
               sm = sm*b
               call daxpy(mml1p2,sm,u(1,l1m1),iue,c(kl1),ice)
               kl1 = kl1 + icv
   70       continue
            u(1,l1m1) = ul1m1
c     ......exit
            if (lpivot .eq. l1m1) go to 140
            kl1 = kl2
            call dswap(ncv,c(kl1),icv,c(klp),icv)
         go to 130
   80    continue
            i2 = 1 - icv + ice*(lpivot - 1)
            incr = ice*(l1 - lpivot)
            do 120 j = 1, ncv
               i2 = i2 + icv
               i3 = i2 + incr
               i4 = i3
               sm = c(i2)*up
               do 90 i = l1, m
                  sm = sm + c(i3)*u(1,i)
                  i3 = i3 + ice
   90          continue
               if (sm .eq. 0.0d0) go to 110
                  sm = sm*b
                  c(i2) = c(i2) + sm*up
                  do 100 i = l1, m
                     c(i4) = c(i4) + sm*u(1,i)
                     i4 = i4 + ice
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
      return
      end
