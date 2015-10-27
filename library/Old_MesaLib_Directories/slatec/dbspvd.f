*deck dbspvd
      subroutine dbspvd (t, k, nderiv, x, ileft, ldvnik, vnikx, work)
c***begin prologue  dbspvd
c***purpose  calculate the value and all derivatives of order less than
c            nderiv of all basis functions which do not vanish at x.
c***library   slatec
c***category  e3, k6
c***type      double precision (bspvd-s, dbspvd-d)
c***keywords  differentiation of b-spline, evaluation of b-spline
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract    **** a double precision routine ****
c
c         dbspvd is the bsplvd routine of the reference.
c
c         dbspvd calculates the value and all derivatives of order
c         less than nderiv of all basis functions which do not
c         (possibly) vanish at x.  ileft is input such that
c         t(ileft) .le. x .lt. t(ileft+1).  a call to intrv(t,n+1,x,
c         ilo,ileft,mflag) will produce the proper ileft.  the output of
c         dbspvd is a matrix vnikx(i,j) of dimension at least (k,nderiv)
c         whose columns contain the k nonzero basis functions and
c         their nderiv-1 right derivatives at x, i=1,k, j=1,nderiv.
c         these basis functions have indices ileft-k+i, i=1,k,
c         k .le. ileft .le. n.  the nonzero part of the i-th basis
c         function lies in (t(i),t(i+k)), i=1,n).
c
c         if x=t(ileft+1) then vnikx contains left limiting values
c         (left derivatives) at t(ileft+1).  in particular, ileft = n
c         produces left limiting values at the right end point
c         x=t(n+1).  to obtain left limiting values at t(i), i=k+1,n+1,
c         set x= next lower distinct knot, call intrv to get ileft,
c         set x=t(i), and then call dbspvd.
c
c     description of arguments
c         input      t,x are double precision
c          t       - knot vector of length n+k, where
c                    n = number of b-spline basis functions
c                    n = sum of knot multiplicities-k
c          k       - order of the b-spline, k .ge. 1
c          nderiv  - number of derivatives = nderiv-1,
c                    1 .le. nderiv .le. k
c          x       - argument of basis functions,
c                    t(k) .le. x .le. t(n+1)
c          ileft   - largest integer such that
c                    t(ileft) .le. x .lt.  t(ileft+1)
c          ldvnik  - leading dimension of matrix vnikx
c
c         output     vnikx,work are double precision
c          vnikx   - matrix of dimension at least (k,nderiv) contain-
c                    ing the nonzero basis functions at x and their
c                    derivatives columnwise.
c          work    - a work vector of length (k+1)*(k+2)/2
c
c     error conditions
c         improper input is a fatal error
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  dbspvn, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbspvd
c
      integer i,ideriv,ileft,ipkmd,j,jj,jlow,jm,jp1mid,k,kmd, kp1, l,
     1 ldummy, m, mhigh, nderiv
      double precision factor, fkmd, t, v, vnikx, work, x
c     dimension t(ileft+k), work((k+1)*(k+2)/2)
c     a(i,j) = work(i+j*(j+1)/2),  i=1,j+1  j=1,k-1
c     a(i,k) = w0rk(i+k*(k-1)/2)  i=1.k
c     work(1) and work((k+1)*(k+2)/2) are not used.
      dimension t(*), vnikx(ldvnik,*), work(*)
c***first executable statement  dbspvd
      if(k.lt.1) go to 200
      if(nderiv.lt.1 .or. nderiv.gt.k) go to 205
      if(ldvnik.lt.k) go to 210
      ideriv = nderiv
      kp1 = k + 1
      jj = kp1 - ideriv
      call dbspvn(t, jj, k, 1, x, ileft, vnikx, work, iwork)
      if (ideriv.eq.1) go to 100
      mhigh = ideriv
      do 20 m=2,mhigh
        jp1mid = 1
        do 10 j=ideriv,k
          vnikx(j,ideriv) = vnikx(jp1mid,1)
          jp1mid = jp1mid + 1
   10   continue
        ideriv = ideriv - 1
        jj = kp1 - ideriv
        call dbspvn(t, jj, k, 2, x, ileft, vnikx, work, iwork)
   20 continue
c
      jm = kp1*(kp1+1)/2
      do 30 l = 1,jm
        work(l) = 0.0d0
   30 continue
c     a(i,i) = work(i*(i+3)/2) = 1.0       i = 1,k
      l = 2
      j = 0
      do 40 i = 1,k
        j = j + l
        work(j) = 1.0d0
        l = l + 1
   40 continue
      kmd = k
      do 90 m=2,mhigh
        kmd = kmd - 1
        fkmd = kmd
        i = ileft
        j = k
        jj = j*(j+1)/2
        jm = jj - j
        do 60 ldummy=1,kmd
          ipkmd = i + kmd
          factor = fkmd/(t(ipkmd)-t(i))
          do 50 l=1,j
            work(l+jj) = (work(l+jj)-work(l+jm))*factor
   50     continue
          i = i - 1
          j = j - 1
          jj = jm
          jm = jm - j
   60   continue
c
        do 80 i=1,k
          v = 0.0d0
          jlow = max(i,m)
          jj = jlow*(jlow+1)/2
          do 70 j=jlow,k
            v = work(i+jj)*vnikx(j,m) + v
            jj = jj + j + 1
   70     continue
          vnikx(i,m) = v
   80   continue
   90 continue
  100 return
c
c
  200 continue
      call xermsg ('slatec', 'dbspvd', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  205 continue
      call xermsg ('slatec', 'dbspvd',
     +   'nderiv does not satisfy 1.le.nderiv.le.k', 2, 1)
      return
  210 continue
      call xermsg ('slatec', 'dbspvd',
     +   'ldvnik does not satisfy ldvnik.ge.k', 2, 1)
      return
      end
