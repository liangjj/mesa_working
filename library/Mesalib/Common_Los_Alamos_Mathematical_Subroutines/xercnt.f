*deck xercnt
      subroutine xercnt (librar, subrou, messg, nerr, level, kontrl)
c***begin prologue  xercnt
c***subsidiary
c***purpose  allow user control over handling of errors.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xercnt-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to xercnt.
c        if the user has provided his own version of xercnt, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        librar - the library that the routine is in.
c        subrou - the subroutine that xermsg is being called from
c        messg  - the first 20 characters of the error message.
c        nerr   - same as in the call to xermsg.
c        level  - same as in the call to xermsg.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
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
c   900510  changed calling sequence to include library and subroutine
c           names, changed routine name from xerctl to xercnt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xercnt
      character*(*) librar, subrou, messg
c***first executable statement  xercnt
      return
      end
