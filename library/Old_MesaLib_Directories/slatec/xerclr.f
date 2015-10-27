*deck xerclr
      subroutine xerclr
c***begin prologue  xerclr
c***purpose  reset current error number to zero.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerclr-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary in order to determine that a certain
c        error has occurred again since the last time numxer was
c        referenced.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerclr
c***first executable statement  xerclr
      junk = j4save(1,0,.true.)
      return
      end
