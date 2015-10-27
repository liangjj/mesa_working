*deck xerhlt
      subroutine xerhlt (messg)
c***begin prologue  xerhlt
c***subsidiary
c***purpose  abort program execution and print error message.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerhlt-a)
c***keywords  abort program execution, error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        ***note*** machine dependent routine
c        xerhlt aborts the execution of the program.
c        the error message causing the abort is given in the calling
c        sequence, in case one needs it for printing on a dayfile,
c        for example.
c
c     description of parameters
c        messg is as in xermsg.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900206  routine changed from user-callable to subsidiary.  (wrb)
c   900510  changed calling sequence to delete length of character
c           and changed routine name from xerabt to xerhlt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerhlt
      character*(*) messg
c***first executable statement  xerhlt
      stop
      end
