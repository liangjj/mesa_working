*deck dbskes
      subroutine dbskes (xnu, x, nin, bke)
c***begin prologue  dbskes
c***purpose  compute a sequence of exponentially scaled modified bessel
c            functions of the third kind of fractional order.
c***library   slatec (fnlib)
c***category  c10b3
c***type      double precision (beskes-s, dbskes-d)
c***keywords  exponentially scaled, fnlib, fractional order,
c             modified bessel function, sequence of bessel functions,
c             special functions, third kind
c***author  fullerton, w., (lanl)
c***description
c
c dbskes(xnu,x,nin,bke) computes a double precision sequence
c of exponentially scaled modified bessel functions
c of the third kind of order xnu + i at x, where x .gt. 0,
c xnu lies in (-1,1), and i = 0, 1, ... , nin - 1, if nin is positive
c and i = 0, -1, ... , nin + 1, if nin is negative.  on return, the
c vector bke(.) contains the results at x for order starting at xnu.
c xnu, x, and bke are double precision.  nin is integer.
c
c***references  (none)
c***routines called  d1mach, d9knus, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbskes
      double precision xnu, x, bke(*), bknu1, v, vincr, vend, alnbig,
     1  d1mach, direct
      save alnbig
      data alnbig / 0.d0 /
c***first executable statement  dbskes
      if (alnbig.eq.0.d0) alnbig = log (d1mach(2))
c
      v = abs(xnu)
      n = abs(nin)
c
      if (v .ge. 1.d0) call xermsg ('slatec', 'dbskes',
     +   'abs(xnu) must be lt 1', 2, 2)
      if (x .le. 0.d0) call xermsg ('slatec', 'dbskes', 'x is le 0', 3,
     +   2)
      if (n .eq. 0) call xermsg ('slatec', 'dbskes',
     +   'n the number in the sequence is 0', 4, 2)
c
      call d9knus (v, x, bke(1), bknu1, iswtch)
      if (n.eq.1) return
c
      vincr = sign (1.0, real(nin))
      direct = vincr
      if (xnu.ne.0.d0) direct = vincr*sign(1.d0, xnu)
      if (iswtch .eq. 1 .and. direct .gt. 0.) call xermsg ('slatec',
     +   'dbskes', 'x so small bessel k-sub-xnu+1 overflows', 5, 2)
      bke(2) = bknu1
c
      if (direct.lt.0.) call d9knus (abs(xnu+vincr), x, bke(2), bknu1,
     1  iswtch)
      if (n.eq.2) return
c
      vend = abs (xnu+nin) - 1.0d0
      if ((vend-.5d0)*log(vend)+0.27d0-vend*(log(x)-.694d0) .gt.
     +   alnbig) call xermsg ('slatec', 'dbskes',
     +      'x so small or abs(nu) so big that bessel k-sub-nu ' //
     +      'overflows', 5, 2)
c
      v = xnu
      do 10 i=3,n
        v = v + vincr
        bke(i) = 2.0d0*v*bke(i-1)/x + bke(i-2)
 10   continue
c
      return
      end
