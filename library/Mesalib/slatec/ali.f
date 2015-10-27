*deck ali
      function ali (x)
c***begin prologue  ali
c***purpose  compute the logarithmic integral.
c***library   slatec (fnlib)
c***category  c5
c***type      single precision (ali-s, dli-d)
c***keywords  fnlib, logarithmic integral, special functions
c***author  fullerton, w., (lanl)
c***description
c
c ali(x) computes the logarithmic integral; i.e., the
c integral from 0.0 to x of (1.0/ln(t))dt.
c
c***references  (none)
c***routines called  ei, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  ali
c***first executable statement  ali
      if (x .le. 0.0) call xermsg ('slatec', 'ali',
     +   'log integral undefined for x le 0', 1, 2)
      if (x .eq. 1.0) call xermsg ('slatec', 'ali',
     +   'log integral undefined for x = 1', 2, 2)
c
      ali = ei (log(x) )
c
      return
      end
