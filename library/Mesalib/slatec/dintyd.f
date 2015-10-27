*deck dintyd
      subroutine dintyd (t, k, yh, nyh, dky, iflag)
c***begin prologue  dintyd
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (intyd-s, dintyd-d)
c***author  watts, h. a., (snla)
c***description
c
c   dintyd approximates the solution and derivatives at t by polynomial
c   interpolation. must be used in conjunction with the integrator
c   package ddebdf.
c ----------------------------------------------------------------------
c dintyd computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.
c this routine is called by ddebdf with k = 0,1 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in lsode usage documentation.)
c ----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c ----------------------------------------------------------------------
c
c***see also  ddebdf
c***routines called  (none)
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dintyd
c
      integer i, ic, ier, iflag, iownd, iowns, j, jb, jb2, jj, jj1,
     1      jp1, jstart, k, kflag, l, maxord, meth, miter, n, nfe,
     2      nje, nq, nqu, nst, nyh
      double precision c, dky, el0, h, hmin, hmxi, hu, r, rownd,
     1      rowns, s, t, tn, tp, uround, yh
      dimension yh(nyh,*),dky(*)
      common /ddebd1/ rownd,rowns(210),el0,h,hmin,hmxi,hu,tn,uround,
     1                iownd(14),iowns(6),ier,jstart,kflag,l,meth,miter,
     2                maxord,n,nq,nst,nfe,nje,nqu
c
c     begin block permitting ...exits to 130
c***first executable statement  dintyd
         iflag = 0
         if (k .lt. 0 .or. k .gt. nq) go to 110
            tp = tn - hu*(1.0d0 + 100.0d0*uround)
            if ((t - tp)*(t - tn) .le. 0.0d0) go to 10
               iflag = -2
c     .........exit
               go to 130
   10       continue
c
            s = (t - tn)/h
            ic = 1
            if (k .eq. 0) go to 30
               jj1 = l - k
               do 20 jj = jj1, nq
                  ic = ic*jj
   20          continue
   30       continue
            c = ic
            do 40 i = 1, n
               dky(i) = c*yh(i,l)
   40       continue
            if (k .eq. nq) go to 90
               jb2 = nq - k
               do 80 jb = 1, jb2
                  j = nq - jb
                  jp1 = j + 1
                  ic = 1
                  if (k .eq. 0) go to 60
                     jj1 = jp1 - k
                     do 50 jj = jj1, j
                        ic = ic*jj
   50                continue
   60             continue
                  c = ic
                  do 70 i = 1, n
                     dky(i) = c*yh(i,jp1) + s*dky(i)
   70             continue
   80          continue
c     .........exit
               if (k .eq. 0) go to 130
   90       continue
            r = h**(-k)
            do 100 i = 1, n
               dky(i) = r*dky(i)
  100       continue
         go to 120
  110    continue
c
            iflag = -1
  120    continue
  130 continue
      return
c     ----------------------- end of subroutine dintyd
c     -----------------------
      end
