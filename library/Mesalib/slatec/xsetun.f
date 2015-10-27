*deck xsetun
      subroutine xsetun (iunit)
c***begin prologue  xsetun
c***purpose  set output file to which error messages are to be sent.
c***library   slatec (xerror)
c***category  r3b
c***type      all (xsetun-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xsetun sets the output file to which error messages are to
c        be sent.  only one file will be used.  see xsetua for
c        how to declare more than one file.
c
c     description of parameter
c      --input--
c        iunit - an input parameter giving the logical unit number
c                to which error messages are to be sent.
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
c***end prologue  xsetun
c***first executable statement  xsetun
      junk = j4save(3,iunit,.true.)
      junk = j4save(5,1,.true.)
      return
      end
