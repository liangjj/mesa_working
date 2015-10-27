*deck xgetun
      subroutine xgetun (iunit)
c***begin prologue  xgetun
c***purpose  return the (first) output file to which error messages
c            are being sent.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xgetun-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xgetun gets the (first) output file to which error messages
c        are being sent.  to find out if more than one file is being
c        used, one must use the xgetua routine.
c
c     description of parameter
c      --output--
c        iunit - the logical unit number of the  (first) unit to
c                which error messages are being sent.
c                a value of zero means that the default file, as
c                defined by the i1mach routine, is being used.
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
c***end prologue  xgetun
c***first executable statement  xgetun
      iunit = j4save(3,0,.false.)
      return
      end
