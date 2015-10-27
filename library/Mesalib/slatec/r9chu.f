*deck r9chu
      function r9chu (a, b, z)
c***begin prologue  r9chu
c***subsidiary
c***purpose  evaluate for large z  z**a * u(a,b,z) where u is the
c            logarithmic confluent hypergeometric function.
c***library   slatec (fnlib)
c***category  c11
c***type      single precision (r9chu-s, d9chu-d)
c***keywords  fnlib, logarithmic confluent hypergeometric function,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate for large z  z**a * u(a,b,z)  where u is the logarithmic
c confluent hypergeometric function.  a rational approximation due to y.
c l. luke is used.  when u is not in the asymptotic region, i.e., when a
c or b is large compared with z, considerable significance loss occurs.
c a warning is provided when the computed result is less than half
c precision.
c
c***references  (none)
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   770801  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9chu
      dimension aa(4), bb(4)
      logical first
      save eps, sqeps, first
      data first /.true./
c***first executable statement  r9chu
      if (first) then
         eps = 4.0*r1mach(4)
         sqeps = sqrt (r1mach(4))
      endif
      first = .false.
c
      bp = 1.0 + a - b
      ab = a*bp
      ct2 = 2.0*(z-ab)
      sab = a + bp
c
      bb(1) = 1.0
      aa(1) = 1.0
c
      ct3 = sab + 1.0 + ab
      bb(2) = 1.0 + 2.0*z/ct3
      aa(2) = 1.0 + ct2/ct3
c
      anbn = ct3 + sab + 3.0
      ct1 = 1.0 + 2.0*z/anbn
      bb(3) = 1.0 + 6.0*ct1*z/ct3
      aa(3) = 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3
c
      do 30 i=4,300
        x2i1 = 2*i - 3
        ct1 = x2i1/(x2i1-2.0)
        anbn = anbn + x2i1 + sab
        ct2 = (x2i1 - 1.0) / anbn
        c2 = x2i1*ct2 - 1.0
        d1z = x2i1*2.0*z/anbn
c
        ct3 = sab*ct2
        g1 = d1z + ct1*(c2+ct3)
        g2 = d1z - c2
        g3 = ct1*(1.0 - ct3 - 2.0*ct2)
c
        bb(4) = g1*bb(3) + g2*bb(2) + g3*bb(1)
        aa(4) = g1*aa(3) + g2*aa(2) + g3*aa(1)
        if (abs(aa(4)*bb(1)-aa(1)*bb(4)).lt.eps*abs(bb(4)*bb(1)))
     1    go to 40
c
c if overflows or underflows prove to be a problem, the statements
c below could be altered to incorporate a dynamically adjusted scale
c factor.
c
        do 20 j=1,3
          bb(j) = bb(j+1)
          aa(j) = aa(j+1)
 20     continue
 30   continue
      call xermsg ('slatec', 'r9chu', 'no convergence in 300 terms', 1,
     +   2)
c
 40   r9chu = aa(4)/bb(4)
c
      if (r9chu .lt. sqeps .or. r9chu .gt. 1.0/sqeps) call xermsg
     +   ('slatec', 'r9chu', 'answer less than half precision', 2, 1)
c
      return
      end
