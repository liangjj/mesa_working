*deck dxset
      subroutine dxset (irad, nradpl, dzero, nbits, ierror)
c***begin prologue  dxset
c***purpose  to provide double-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      double precision (xset-s, dxset-d)
c***keywords  extended-range double-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c
c   subroutine  dxset  must be called prior to calling any other
c extended-range subroutine. it calculates and stores several
c machine-dependent constants in common blocks. the user must
c supply four constants that pertain to his particular computer.
c the constants are
c
c          irad = the internal base of double-precision
c                 arithmetic in the computer.
c        nradpl = the number of radix places carried in
c                 the double-precision representation.
c         dzero = the smallest of 1/dmin, dmax, dmaxln where
c                 dmin = the smallest positive double-precision
c                 number or an upper bound to this number,
c                 dmax = the largest double-precision number
c                 or a lower bound to this number,
c                 dmaxln = the largest double-precision number
c                 such that log10(dmaxln) can be computed by the
c                 fortran system (on most systems dmaxln = dmax).
c         nbits = the number of bits (exclusive of sign) in
c                 an integer computer word.
c
c alternatively, any or all of the constants can be given
c the value 0 (0.0d0 for dzero). if a constant is zero, dxset tries
c to assign an appropriate value by calling i1mach
c (see p.a.fox, a.d.hall, n.l.schryer, algorithm 528 framework
c for a portable library, acm transactions on math software,
c v.4, no.2, june 1978, 177-188).
c
c   this is the setting-up subroutine for a package of subroutines
c that facilitate the use of extended-range arithmetic. extended-range
c arithmetic on a particular computer is defined on the set of numbers
c of the form
c
c               (x,ix) = x*radix**ix
c
c where x is a double-precision number called the principal part,
c ix is an integer called the auxiliary index, and radix is the
c internal base of the double-precision arithmetic.  obviously,
c each real number is representable without error by more than one
c extended-range form.  conversions between  different forms are
c essential in carrying out arithmetic operations.  with the choice
c of radix we have made, and the subroutines we have written, these
c conversions are performed without error (at least on most computers).
c (see smith, j.m., olver, f.w.j., and lozier, d.w., extended-range
c arithmetic and normalized legendre polynomials, acm transactions on
c mathematical software, march 1981).
c
c   an extended-range number  (x,ix)  is said to be in adjusted form if
c x and ix are zero or
c
c           radix**(-l) .le. abs(x) .lt. radix**l
c
c is satisfied, where l is a computer-dependent integer defined in this
c subroutine. two extended-range numbers in adjusted form can be added,
c subtracted, multiplied or divided (if the divisor is nonzero) without
c causing overflow or underflow in the principal part of the result.
c with proper use of the extended-range subroutines, the only overflow
c that can occur is integer overflow in the auxiliary index. if this
c is detected, the software calls xerror (a general error-handling
c fortran subroutine package).
c
c   multiplication and division is performed by setting
c
c                 (x,ix)*(y,iy) = (x*y,ix+iy)
c or
c                 (x,ix)/(y,iy) = (x/y,ix-iy).
c
c pre-adjustment of the operands is essential to avoid
c overflow or  underflow of the principal part. subroutine
c dxadj (see below) may be called to transform any extended-
c range number into adjusted form.
c
c   addition and subtraction require the use of subroutine dxadd
c (see below).  the input operands need not be in adjusted form.
c however, the result of addition or subtraction is returned
c in adjusted form.  thus, for example, if (x,ix),(y,iy),
c (u,iu),  and (v,iv) are in adjusted form, then
c
c                 (x,ix)*(y,iy) + (u,iu)*(v,iv)
c
c can be computed and stored in adjusted form with no explicit
c calls to dxadj.
c
c   when an extended-range number is to be printed, it must be
c converted to an extended-range form with decimal radix.  subroutine
c dxcon is provided for this purpose.
c
c   the subroutines contained in this package are
c
c     subroutine dxadd
c usage
c                  call dxadd(x,ix,y,iy,z,iz,ierror)
c                  if (ierror.ne.0) return
c description
c                  forms the extended-range sum  (z,iz) =
c                  (x,ix) + (y,iy).  (z,iz) is adjusted
c                  before returning. the input operands
c                  need not be in adjusted form, but their
c                  principal parts must satisfy
c                  radix**(-2l).le.abs(x).le.radix**(2l),
c                  radix**(-2l).le.abs(y).le.radix**(2l).
c
c     subroutine dxadj
c usage
c                  call dxadj(x,ix,ierror)
c                  if (ierror.ne.0) return
c description
c                  transforms (x,ix) so that
c                  radix**(-l) .le. abs(x) .lt. radix**l.
c                  on most computers this transformation does
c                  not change the mantissa of x provided radix is
c                  the number base of double-precision arithmetic.
c
c     subroutine dxc210
c usage
c                  call dxc210(k,z,j,ierror)
c                  if (ierror.ne.0) return
c description
c                  given k this subroutine computes j and z
c                  such that  radix**k = z*10**j, where z is in
c                  the range 1/10 .le. z .lt. 1.
c                  the value of z will be accurate to full
c                  double-precision provided the number
c                  of decimal places in the largest
c                  integer plus the number of decimal
c                  places carried in double-precision does not
c                  exceed 60. dxc210 is called by subroutine
c                  dxcon when necessary. the user should
c                  never need to call dxc210 directly.
c
c     subroutine dxcon
c usage
c                  call dxcon(x,ix,ierror)
c                  if (ierror.ne.0) return
c description
c                  converts (x,ix) = x*radix**ix
c                  to decimal form in preparation for
c                  printing, so that (x,ix) = x*10**ix
c                  where 1/10 .le. abs(x) .lt. 1
c                  is returned, except that if
c                  (abs(x),ix) is between radix**(-2l)
c                  and radix**(2l) then the reduced
c                  form with ix = 0 is returned.
c
c     subroutine dxred
c usage
c                  call dxred(x,ix,ierror)
c                  if (ierror.ne.0) return
c description
c                  if
c                  radix**(-2l) .le. (abs(x),ix) .le. radix**(2l)
c                  then dxred transforms (x,ix) so that ix=0.
c                  if (x,ix) is outside the above range,
c                  then dxred takes no action.
c                  this subroutine is useful if the
c                  results of extended-range calculations
c                  are to be used in subsequent ordinary
c                  double-precision calculations.
c
c***references  smith, olver and lozier, extended-range arithmetic and
c                 normalized legendre polynomials, acm trans on math
c                 softw, v 7, n 1, march 1981, pp 93--105.
c***routines called  i1mach, xermsg
c***common blocks    dxblk1, dxblk2, dxblk3
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c           calls to xerror changed to calls to xermsg.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxset
      integer irad, nradpl, nbits
      double precision dzero, dzerox
      common /dxblk1/ nbitsf
      save /dxblk1/
      double precision radix, radixl, rad2l, dlg10r
      integer l, l2, kmax
      common /dxblk2/ radix, radixl, rad2l, dlg10r, l, l2, kmax
      save /dxblk2/
      integer nlg102, mlg102, lg102
      common /dxblk3/ nlg102, mlg102, lg102(21)
      save /dxblk3/
      integer iflag
      save iflag
c
      dimension log102(20), lgtemp(20)
      save log102
c
c   log102 contains the first 60 digits of log10(2) for use in
c conversion of extended-range numbers to base 10 .
      data log102 /301,029,995,663,981,195,213,738,894,724,493,026,768,
     * 189,881,462,108,541,310,428/
c
c following coding prevents dxset from being executed more than once.
c this is important because some subroutines (such as dxnrmp and
c dxlegf) call dxset to make sure extended-range arithmetic has
c been initialized. the user may want to pre-empt this call, for
c example when i1mach is not available. see coding below.
      data iflag /0/
c***first executable statement  dxset
      ierror=0
      if (iflag .ne. 0) return
      iradx = irad
      nrdplc = nradpl
      dzerox = dzero
      iminex = 0
      imaxex = 0
      nbitsx = nbits
c following 5 statements should be deleted if i1mach is
c not available or not configured to return the correct
c machine-dependent values.
      if (iradx .eq. 0) iradx = i1mach (10)
      if (nrdplc .eq. 0) nrdplc = i1mach (14)
      if (dzerox .eq. 0.0d0) iminex = i1mach (15)
      if (dzerox .eq. 0.0d0) imaxex = i1mach (16)
      if (nbitsx .eq. 0) nbitsx = i1mach (8)
      if (iradx.eq.2) go to 10
      if (iradx.eq.4) go to 10
      if (iradx.eq.8) go to 10
      if (iradx.eq.16) go to 10
      call xermsg ('slatec', 'dxset', 'improper value of irad', 201, 1)
      ierror=201
      return
   10 continue
      log2r=0
      if (iradx.eq.2) log2r = 1
      if (iradx.eq.4) log2r = 2
      if (iradx.eq.8) log2r = 3
      if (iradx.eq.16) log2r = 4
      nbitsf=log2r*nrdplc
      radix = iradx
      dlg10r = log10(radix)
      if (dzerox .ne. 0.0d0) go to 14
      lx = min ((1-iminex)/2, (imaxex-1)/2)
      go to 16
   14 lx = 0.5d0*log10(dzerox)/dlg10r
c radix**(2*l) should not overflow, but reduce l by 1 for further
c protection.
      lx=lx-1
   16 l2 = 2*lx
      if (lx.ge.4) go to 20
      call xermsg ('slatec', 'dxset', 'improper value of dzero', 202, 1)
      ierror=202
      return
   20 l = lx
      radixl = radix**l
      rad2l = radixl**2
c    it is necessary to restrict nbits (or nbitsx) to be less than some
c upper limit because of binary-to-decimal conversion. such conversion
c is done by dxc210 and requires a constant that is stored to some fixed
c precision. the stored constant (log102 in this routine) provides
c for conversions accurate to the last decimal digit when the integer
c word length does not exceed 63. a lower limit of 15 bits is imposed
c because the software is designed to run on computers with integer word
c length of at least 16 bits.
      if (15.le.nbitsx .and. nbitsx.le.63) go to 30
      call xermsg ('slatec', 'dxset', 'improper value of nbits', 203, 1)
      ierror=203
      return
   30 continue
      kmax = 2**(nbitsx-1) - l2
      nb = (nbitsx-1)/2
      mlg102 = 2**nb
      if (1.le.nrdplc*log2r .and. nrdplc*log2r.le.120) go to 40
      call xermsg ('slatec', 'dxset', 'improper value of nradpl', 204,
     +             1)
      ierror=204
      return
   40 continue
      nlg102 = nrdplc*log2r/nb + 3
      np1 = nlg102 + 1
c
c   after completion of the following loop, ic contains
c the integer part and lgtemp contains the fractional part
c of log10(iradx) in radix 1000.
      ic = 0
      do 50 ii=1,20
        i = 21 - ii
        it = log2r*log102(i) + ic
        ic = it/1000
        lgtemp(i) = mod(it,1000)
   50 continue
c
c   after completion of the following loop, lg102 contains
c log10(iradx) in radix mlg102. the radix point is
c between lg102(1) and lg102(2).
      lg102(1) = ic
      do 80 i=2,np1
        lg102x = 0
        do 70 j=1,nb
          ic = 0
          do 60 kk=1,20
            k = 21 - kk
            it = 2*lgtemp(k) + ic
            ic = it/1000
            lgtemp(k) = mod(it,1000)
   60     continue
          lg102x = 2*lg102x + ic
   70   continue
        lg102(i) = lg102x
   80 continue
c
c check special conditions required by subroutines...
      if (nrdplc.lt.l) go to 90
      call xermsg ('slatec', 'dxset', 'nradpl .ge. l', 205, 1)
      ierror=205
      return
   90 if (6*l.le.kmax) go to 100
      call xermsg ('slatec', 'dxset', '6*l .gt. kmax', 206, 1)
      ierror=206
      return
  100 continue
      iflag = 1
      return
      end
