*deck beskes
      subroutine beskes (xnu, x, nin, bke)
c***begin prologue  beskes
c***purpose  compute a sequence of exponentially scaled modified bessel
c            functions of the third kind of fractional order.
c***library   slatec (fnlib)
c***category  c10b3
c***type      single precision (beskes-s, dbskes-d)
c***keywords  exponentially scaled, fnlib, fractional order,
c             modified bessel function, sequence of bessel functions,
c             special functions, third kind
c***author  fullerton, w., (lanl)
c***description
c
c beskes computes a sequence of exponentially scaled
c (i.e., multipled by exp(x)) modified bessel
c functions of the third kind of order xnu + i at x, where x .gt. 0,
c xnu lies in (-1,1), and i = 0, 1, ... , nin - 1, if nin is positive
c and i = 0, -1, ... , nin + 1, if nin is negative.  on return, the
c vector bke(.) contains the results at x for order starting at xnu.
c
c***references  (none)
c***routines called  r1mach, r9knus, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  beskes
      dimension bke(*)
      save alnbig
      data alnbig / 0. /
c***first executable statement  beskes
      if (alnbig.eq.0.) alnbig = log (r1mach(2))
c
      v = abs(xnu)
      n = abs(nin)
c
      if (v .ge. 1.) call xermsg ('slatec', 'beskes',
     +   'abs(xnu) must be lt 1', 2, 2)
      if (x .le. 0.) call xermsg ('slatec', 'beskes', 'x is le 0', 3,
     +   2)
      if (n .eq. 0) call xermsg ('slatec', 'beskes',
     +   'n the number in the sequence is 0', 4, 2)
c
      call r9knus (v, x, bke(1), bknu1, iswtch)
      if (n.eq.1) return
c
      vincr = sign (1.0, real(nin))
      direct = vincr
      if (xnu.ne.0.) direct = vincr*sign(1.0,xnu)
      if (iswtch .eq. 1 .and. direct .gt. 0.) call xermsg ('slatec',
     +   'beskes', 'x so small bessel k-sub-xnu+1 overflows', 5, 2)
      bke(2) = bknu1
c
      if (direct.lt.0.) call r9knus (abs(xnu+vincr), x, bke(2), bknu1,
     1  iswtch)
      if (n.eq.2) return
c
      vend = abs(xnu+nin) - 1.0
      if ((vend-0.5)*log(vend)+0.27-vend*(log(x)-.694) .gt. alnbig)
     1call xermsg ( 'slatec', 'beskes',
     2'x so small or abs(nu) so big that bessel k-sub-nu overflows',
     35, 2)
c
      v = xnu
      do 10 i=3,n
        v = v + vincr
        bke(i) = 2.0*v*bke(i-1)/x + bke(i-2)
 10   continue
c
      return
      end
