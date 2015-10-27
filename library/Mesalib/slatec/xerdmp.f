*deck xerdmp
      subroutine xerdmp
c***begin prologue  xerdmp
c***purpose  print the error tables and then clear them.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerdmp-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xerdmp prints the error tables, then clears them.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  xersve
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  changed call of xersav to xersve.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerdmp
c***first executable statement  xerdmp
      call xersve (' ',' ',' ',0,0,0,kount)
      return
      end
