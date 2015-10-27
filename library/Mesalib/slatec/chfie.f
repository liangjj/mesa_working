*deck chfie
      real function chfie (x1, x2, f1, f2, d1, d2, a, b)
c***begin prologue  chfie
c***subsidiary
c***purpose  evaluates integral of a single cubic for pchia
c***library   slatec (pchip)
c***type      single precision (chfie-s, dchfie-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c          chfie:  cubic hermite function integral evaluator.
c
c     called by  pchia  to evaluate the integral of a single cubic (in
c     hermite form) over an arbitrary interval (a,b).
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        real  x1, x2, f1, f2, d1, d2, a, b
c        real  value, chfie
c
c        value = chfie (x1, x2, f1, f2, d1, d2, a, b)
c
c   parameters:
c
c     value -- (output) value of the requested integral.
c
c     x1,x2 -- (input) endpoints if interval of definition of cubic.
c
c     f1,f2 -- (input) function values at the ends of the interval.
c
c     d1,d2 -- (input) derivative values at the ends of the interval.
c
c     a,b -- (input) endpoints of interval of integration.
c
c***see also  pchia
c***routines called  (none)
c***revision history  (yymmdd)
c   820730  date written
c   820805  converted to slatec library version.
c   870813  minor cosmetic changes.
c   890411  1. added save statements (vers. 3.2).
c           2. added six to real declaration.
c   890411  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated author section in prologue.  (wrb)
c   930503  corrected to set value=0 when ierr.ne.0.  (fnf)
c   930504  eliminated ierr and changed name from chfiv to chfie.  (fnf)
c***end prologue  chfie
c
c  programming notes:
c  1. there is no error return from this routine because zero is
c     indeed the mathematically correct answer when x1.eq.x2 .
c**end
c
c  declare arguments.
c
      real  x1, x2, f1, f2, d1, d2, a, b
c
c  declare local variables.
c
      real  dterm, four, fterm, h, half, phia1, phia2, phib1, phib2,
     *      psia1, psia2, psib1, psib2, six, ta1, ta2, tb1, tb2, three,
     *      two, ua1, ua2, ub1, ub2
      save half, two, three, four, six
c
c  initialize.
c
      data  half /0.5/,  two /2./,  three /3./,  four /4./,  six /6./
c
c  validity check input.
c
c***first executable statement  chfie
      if (x1 .eq. x2)  then
         chfie = 0
      else
         h = x2 - x1
         ta1 = (a - x1) / h
         ta2 = (x2 - a) / h
         tb1 = (b - x1) / h
         tb2 = (x2 - b) / h
c
         ua1 = ta1**3
         phia1 = ua1 * (two - ta1)
         psia1 = ua1 * (three*ta1 - four)
         ua2 = ta2**3
         phia2 =  ua2 * (two - ta2)
         psia2 = -ua2 * (three*ta2 - four)
c
         ub1 = tb1**3
         phib1 = ub1 * (two - tb1)
         psib1 = ub1 * (three*tb1 - four)
         ub2 = tb2**3
         phib2 =  ub2 * (two - tb2)
         psib2 = -ub2 * (three*tb2 - four)
c
         fterm =   f1*(phia2 - phib2) + f2*(phib1 - phia1)
         dterm = ( d1*(psia2 - psib2) + d2*(psib1 - psia1) )*(h/six)
c
         chfie = (half*h) * (fterm + dterm)
      endif
c
      return
c------------- last line of chfie follows ------------------------------
      end
