*deck xgetf
      subroutine xgetf (kontrl)
c***begin prologue  xgetf
c***purpose  return the current value of the error control flag.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xgetf-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c   abstract
c        xgetf returns the current value of the error control flag
c        in kontrl.  see subroutine xsetf for flag value meanings.
c        (kontrl is an output parameter only.)
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
c***end prologue  xgetf
c***first executable statement  xgetf
      kontrl = j4save(2,0,.false.)
      return
      end
