*deck dbspvn
      subroutine dbspvn (t, jhigh, k, index, x, ileft, vnikx, work,
     +   iwork)
c***begin prologue  dbspvn
c***purpose  calculate the value of all (possibly) nonzero basis
c            functions at x.
c***library   slatec
c***category  e3, k6
c***type      double precision (bspvn-s, dbspvn-d)
c***keywords  evaluation of b-spline
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract    **** a double precision routine ****
c         dbspvn is the bsplvn routine of the reference.
c
c         dbspvn calculates the value of all (possibly) nonzero basis
c         functions at x of order max(jhigh,(j+1)*(index-1)), where t(k)
c         .le. x .le. t(n+1) and j=iwork is set inside the routine on
c         the first call when index=1.  ileft is such that t(ileft) .le.
c         x .lt. t(ileft+1).  a call to dintrv(t,n+1,x,ilo,ileft,mflag)
c         produces the proper ileft.  dbspvn calculates using the basic
c         algorithm needed in dbspvd.  if only basis functions are
c         desired, setting jhigh=k and index=1 can be faster than
c         calling dbspvd, but extra coding is required for derivatives
c         (index=2) and dbspvd is set up for this purpose.
c
c         left limiting values are set up as described in dbspvd.
c
c     description of arguments
c
c         input      t,x are double precision
c          t       - knot vector of length n+k, where
c                    n = number of b-spline basis functions
c                    n = sum of knot multiplicities-k
c          jhigh   - order of b-spline, 1 .le. jhigh .le. k
c          k       - highest possible order
c          index   - index = 1 gives basis functions of order jhigh
c                          = 2 denotes previous entry with work, iwork
c                              values saved for subsequent calls to
c                              dbspvn.
c          x       - argument of basis functions,
c                    t(k) .le. x .le. t(n+1)
c          ileft   - largest integer such that
c                    t(ileft) .le. x .lt.  t(ileft+1)
c
c         output     vnikx, work are double precision
c          vnikx   - vector of length k for spline values.
c          work    - a work vector of length 2*k
c          iwork   - a work parameter.  both work and iwork contain
c                    information necessary to continue for index = 2.
c                    when index = 1 exclusively, these are scratch
c                    variables and can be used for other purposes.
c
c     error conditions
c         improper input is a fatal error.
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbspvn
c
      integer ileft, imjp1, index, ipj, iwork, jhigh, jp1, jp1ml, k, l
      double precision t, vm, vmprev, vnikx, work, x
c     dimension t(ileft+jhigh)
      dimension t(*), vnikx(*), work(*)
c     content of j, deltam, deltap is expected unchanged between calls.
c     work(i) = deltap(i), work(k+i) = deltam(i), i = 1,k
c***first executable statement  dbspvn
      if(k.lt.1) go to 90
      if(jhigh.gt.k .or. jhigh.lt.1) go to 100
      if(index.lt.1 .or. index.gt.2) go to 105
      if(x.lt.t(ileft) .or. x.gt.t(ileft+1)) go to 110
      go to (10, 20), index
   10 iwork = 1
      vnikx(1) = 1.0d0
      if (iwork.ge.jhigh) go to 40
c
   20 ipj = ileft + iwork
      work(iwork) = t(ipj) - x
      imjp1 = ileft - iwork + 1
      work(k+iwork) = x - t(imjp1)
      vmprev = 0.0d0
      jp1 = iwork + 1
      do 30 l=1,iwork
        jp1ml = jp1 - l
        vm = vnikx(l)/(work(l)+work(k+jp1ml))
        vnikx(l) = vm*work(l) + vmprev
        vmprev = vm*work(k+jp1ml)
   30 continue
      vnikx(jp1) = vmprev
      iwork = jp1
      if (iwork.lt.jhigh) go to 20
c
   40 return
c
c
   90 continue
      call xermsg ('slatec', 'dbspvn', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  100 continue
      call xermsg ('slatec', 'dbspvn',
     +   'jhigh does not satisfy 1.le.jhigh.le.k', 2, 1)
      return
  105 continue
      call xermsg ('slatec', 'dbspvn', 'index is not 1 or 2', 2, 1)
      return
  110 continue
      call xermsg ('slatec', 'dbspvn',
     +   'x does not satisfy t(ileft).le.x.le.t(ileft+1)', 2, 1)
      return
      end
