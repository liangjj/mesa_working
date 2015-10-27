*deck xcon
      subroutine xcon (x, ix, ierror)
c***begin prologue  xcon
c***purpose  to provide single-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      single precision (xcon-s, dxcon-d)
c***keywords  extended-range single-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c     real x
c     integer ix
c
c                  converts (x,ix) = x*radix**ix
c                  to decimal form in preparation for
c                  printing, so that (x,ix) = x*10**ix
c                  where 1/10 .le. abs(x) .lt. 1
c                  is returned, except that if
c                  (abs(x),ix) is between radix**(-2l)
c                  and radix**(2l) then the reduced
c                  form with ix = 0 is returned.
c
c***see also  xset
c***references  (none)
c***routines called  xadj, xc210, xred
c***common blocks    xblk2
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  xcon
      real x
      integer ix
c
c   the conditions imposed on l and kmax by this subroutine
c are
c    (1) 4 .le. l .le. 2**nbits - 1 - kmax
c
c    (2) kmax .le. ((2**nbits)-2)/log10r - l
c
c these conditions must be met by appropriate coding
c in subroutine xset.
c
      real radix, radixl, rad2l, dlg10r
      integer l, l2, kmax
      common /xblk2/ radix, radixl, rad2l, dlg10r, l, l2, kmax
      save /xblk2/, ispace
c
      real a, b, z
c
      data ispace /1/
c   the parameter ispace is the increment used in form-
c ing the auxiliary index of the decimal extended-range
c form. the returned value of ix will be an integer mult-
c iple of ispace. ispace must satisfy 1 .le. ispace .le.
c l/2. if a value greater than 1 is taken, the returned
c value of x will satisfy 10**(-ispace) .le. abs(x) .le. 1
c when (abs(x),ix) .lt. radix**(-2l) and 1/10 .le. abs(x)
c .lt. 10**(ispace-1) when (abs(x),ix) .gt. radix**(2l).
c
c***first executable statement  xcon
      ierror=0
      call xred(x, ix,ierror)
      if (ierror.ne.0) return
      if (ix.eq.0) go to 150
      call xadj(x, ix,ierror)
      if (ierror.ne.0) return
c
c case 1 is when (x,ix) is less than radix**(-2l) in magnitude,
c case 2 is when (x,ix) is greater than radix**(2l) in magnitude.
      itemp = 1
      icase = (3+sign(itemp,ix))/2
      go to (10, 20), icase
   10 if (abs(x).lt.1.0) go to 30
      x = x/radixl
      ix = ix + l
      go to 30
   20 if (abs(x).ge.1.0) go to 30
      x = x*radixl
      ix = ix - l
   30 continue
c
c at this point, radix**(-l) .le. abs(x) .lt. 1.0     in case 1,
c                      1.0 .le. abs(x) .lt. radix**l  in case 2.
      i = log10(abs(x))/dlg10r
      a = radix**i
      go to (40, 60), icase
   40 if (a.le.radix*abs(x)) go to 50
      i = i - 1
      a = a/radix
      go to 40
   50 if (abs(x).lt.a) go to 80
      i = i + 1
      a = a*radix
      go to 50
   60 if (a.le.abs(x)) go to 70
      i = i - 1
      a = a/radix
      go to 60
   70 if (abs(x).lt.radix*a) go to 80
      i = i + 1
      a = a*radix
      go to 70
   80 continue
c
c at this point i is such that
c radix**(i-1) .le. abs(x) .lt. radix**i      in case 1,
c     radix**i .le. abs(x) .lt. radix**(i+1)  in case 2.
      itemp = ispace/dlg10r
      a = radix**itemp
      b = 10.0**ispace
   90 if (a.le.b) go to 100
      itemp = itemp - 1
      a = a/radix
      go to 90
  100 if (b.lt.a*radix) go to 110
      itemp = itemp + 1
      a = a*radix
      go to 100
  110 continue
c
c at this point itemp is such that
c radix**itemp .le. 10**ispace .lt. radix**(itemp+1).
      if (itemp.gt.0) go to 120
c itemp = 0 if, and only if, ispace = 1 and radix = 16.0
      x = x*radix**(-i)
      ix = ix + i
      call xc210(ix, z, j,ierror)
      if (ierror.ne.0) return
      x = x*z
      ix = j
      go to (130, 140), icase
  120 continue
      i1 = i/itemp
      x = x*radix**(-i1*itemp)
      ix = ix + i1*itemp
c
c at this point,
c radix**(-itemp) .le. abs(x) .lt. 1.0        in case 1,
c           1.0 .le. abs(x) .lt. radix**itemp in case 2.
      call xc210(ix, z, j,ierror)
      if (ierror.ne.0) return
      j1 = j/ispace
      j2 = j - j1*ispace
      x = x*z*10.0**j2
      ix = j1*ispace
c
c at this point,
c  10.0**(-2*ispace) .le. abs(x) .lt. 1.0                in case 1,
c           10.0**-1 .le. abs(x) .lt. 10.0**(2*ispace-1) in case 2.
      go to (130, 140), icase
  130 if (b*abs(x).ge.1.0) go to 150
      x = x*b
      ix = ix - ispace
      go to 130
  140 if (10.0*abs(x).lt.b) go to 150
      x = x/b
      ix = ix + ispace
      go to 140
  150 return
      end
