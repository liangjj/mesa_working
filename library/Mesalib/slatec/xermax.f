*deck xermax
      subroutine xermax (max)
c***begin prologue  xermax
c***purpose  set maximum number of times any error message is to be
c            printed.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xermax-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xermax sets the maximum number of times any message
c        is to be printed.  that is, non-fatal messages are
c        not to be printed after they have occurred max times.
c        such non-fatal messages may be printed less than
c        max times even if they occur max times, if error
c        suppression mode (kontrl=0) is ever in effect.
c
c     description of parameter
c      --input--
c        max - the maximum number of times any one message
c              is to be printed.
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
c***end prologue  xermax
c***first executable statement  xermax
      junk = j4save(4,max,.true.)
      return
      end
