*deck dppqad
      subroutine dppqad (ldc, c, xi, lxi, k, x1, x2, pquad)
c***begin prologue  dppqad
c***purpose  compute the integral on (x1,x2) of a k-th order b-spline
c            using the piecewise polynomial (pp) representation.
c***library   slatec
c***category  h2a2a1, e3, k6
c***type      double precision (ppqad-s, dppqad-d)
c***keywords  b-spline, data fitting, interpolation, quadrature, splines
c***author  amos, d. e., (snla)
c***description
c
c     abstract    **** a double precision routine ****
c         dppqad computes the integral on (x1,x2) of a k-th order
c         b-spline using the piecewise polynomial representation
c         (c,xi,lxi,k).  here the taylor expansion about the left
c         end point xi(j) of the j-th interval is integrated and
c         evaluated on subintervals of (x1,x2) which are formed by
c         included break points.  integration outside (xi(1),xi(lxi+1))
c         is permitted.
c
c     description of arguments
c         input      c,xi,x1,x2 are double precision
c           ldc    - leading dimension of matrix c, ldc .ge. k
c           c(i,j) - right taylor derivatives at xi(j), i=1,k , j=1,lxi
c           xi(*)  - break point array of length lxi+1
c           lxi    - number of polynomial pieces
c           k      - order of b-spline, k .ge. 1
c           x1,x2  - end points of quadrature interval, normally in
c                    xi(1) .le. x .le. xi(lxi+1)
c
c         output     pquad is double precision
c           pquad  - integral of the pp representation over (x1,x2)
c
c     error conditions
c         improper input is a fatal error
c
c***references  d. e. amos, quadrature subroutines for splines and
c                 b-splines, report sand79-1825, sandia laboratories,
c                 december 1979.
c***routines called  dintrv, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dppqad
c
      integer i, ii, il, ilo, il1, il2, im, k, ldc, left, lxi, mf1, mf2
      double precision a,aa,bb,c,dx,flk,pquad,q,s,ss,ta,tb,x,xi,x1,x2
      dimension xi(*), c(ldc,*), ss(2)
c
c***first executable statement  dppqad
      pquad = 0.0d0
      if(k.lt.1) go to 100
      if(lxi.lt.1) go to 105
      if(ldc.lt.k) go to 110
      aa = min(x1,x2)
      bb = max(x1,x2)
      if (aa.eq.bb) return
      ilo = 1
      call dintrv(xi, lxi, aa, ilo, il1, mf1)
      call dintrv(xi, lxi, bb, ilo, il2, mf2)
      q = 0.0d0
      do 40 left=il1,il2
        ta = xi(left)
        a = max(aa,ta)
        if (left.eq.1) a = aa
        tb = bb
        if (left.lt.lxi) tb = xi(left+1)
        x = min(bb,tb)
        do 30 ii=1,2
          ss(ii) = 0.0d0
          dx = x - xi(left)
          if (dx.eq.0.0d0) go to 20
          s = c(k,left)
          flk = k
          im = k - 1
          il = im
          do 10 i=1,il
            s = s*dx/flk + c(im,left)
            im = im - 1
            flk = flk - 1.0d0
   10     continue
          ss(ii) = s*dx
   20     continue
          x = a
   30   continue
        q = q + (ss(1)-ss(2))
   40 continue
      if (x1.gt.x2) q = -q
      pquad = q
      return
c
c
  100 continue
      call xermsg ('slatec', 'dppqad', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'dppqad', 'lxi does not satisfy lxi.ge.1',
     +   2, 1)
      return
  110 continue
      call xermsg ('slatec', 'dppqad', 'ldc does not satisfy ldc.ge.k',
     +   2, 1)
      return
      end
