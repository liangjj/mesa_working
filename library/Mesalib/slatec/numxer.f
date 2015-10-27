*deck numxer
      function numxer (nerr)
c***begin prologue  numxer
c***purpose  return the most recent error number.
c***library   slatec (xerror)
c***category  r3c
c***type      integer (numxer-i)
c***keywords  error number, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        numxer returns the most recent error number,
c        in both numxer and the parameter nerr.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save
c***revision history  (yymmdd)
c   790801  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   910411  made user-callable and added keywords section.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  numxer
c***first executable statement  numxer
      nerr = j4save(1,0,.false.)
      numxer = nerr
      return
      end
