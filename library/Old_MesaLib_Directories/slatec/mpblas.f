*deck mpblas
      subroutine mpblas (i1)
c***begin prologue  mpblas
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpblas-a)
c***author  (unknown)
c***description
c
c     this subroutine is called to set up brent's 'mp' package
c     for use by the extended precision inner products from the blas.
c
c     in the slatec library we require the extended precision mp number
c     to have a mantissa twice as long as double precision numbers.
c     the calculation of mpt (and mpmxr which is the actual array size)
c     in this routine will give 2x (or slightly more) on the machine
c     that we are running on.  the integer array size of 30 was chosen
c     to be slightly longer than the longest integer array needed on
c     any machine that we are currently aware of.
c
c***see also  dqdota, dqdoti
c***references  r. p. brent, a fortran multiple-precision arithmetic
c                 package, acm transactions on mathematical software 4,
c                 1 (march 1978), pp. 57-70.
c               r. p. brent, mp, a fortran multiple-precision arithmetic
c                 package, algorithm 524, acm transactions on mathema-
c                 tical software 4, 1 (march 1978), pp. 71-81.
c***routines called  i1mach, xermsg
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920501  reformatted the references section.  (wrb)
c   930124  increased array size in mpcon for sun -r8, and calculate
c               size for quad precision for 2x dp.  (rwc)
c***end prologue  mpblas
      common /mpcom/ mpb, mpt, mpm, mplun, mpmxr, mpr(30)
c***first executable statement  mpblas
      i1 = 1
c
c     for full extended precision accuracy, mpb should be as large as
c     possible, subject to the restrictions in brent's paper.
c
c     statements below are for an integer wordlength of  48, 36, 32,
c     24, 18, and 16.  pick one, or generate a new one.
c       48     mpb = 4194304
c       36     mpb =   65536
c       32     mpb =   16384
c       24     mpb =    1024
c       18     mpb =     128
c       16     mpb =      64
c
      mpbexp = i1mach(8)/2-2
      mpb = 2**mpbexp
c
c     set up remaining parameters
c                  unit for error messages
      mplun = i1mach(4)
c                  number of mp digits
      mpt = (2*i1mach(14)+mpbexp-1)/mpbexp
c                  dimension of r
      mpmxr = mpt+4
c
      if (mpmxr.gt.30) then
         call xermsg('slatec', 'mpblas',
     *      'array space not sufficient for quad precision 2x ' //
     *      'double precision, proceeding.', 1, 1)
         mpt = 26
         mpmxr = 30
      endif
c                  exponent range
      mpm = min(32767,i1mach(9)/4-1)
      return
      end
