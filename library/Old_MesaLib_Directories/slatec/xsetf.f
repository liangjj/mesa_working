*deck xsetf
      subroutine xsetf (kontrl)
c***begin prologue  xsetf
c***purpose  set the error control flag.
c***library   slatec (xerror)
c***category  r3a
c***type      all (xsetf-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xsetf sets the error control flag value to kontrl.
c        (kontrl is an input parameter only.)
c        the following table shows how each message is treated,
c        depending on the values of kontrl and level.  (see xermsg
c        for description of level.)
c
c        if kontrl is zero or negative, no information other than the
c        message itself (including numeric values, if any) will be
c        printed.  if kontrl is positive, introductory messages,
c        trace-backs, etc., will be printed in addition to the message.
c
c              abs(kontrl)
c        level        0              1              2
c        value
c          2        fatal          fatal          fatal
c
c          1     not printed      printed         fatal
c
c          0     not printed      printed        printed
c
c         -1     not printed      printed        printed
c                                  only           only
c                                  once           once
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save, xermsg
c***revision history  (yymmdd)
c   790801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  change call to xerrwv to xermsg.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xsetf
      character *8 xern1
c***first executable statement  xsetf
      if (abs(kontrl) .gt. 2) then
         write (xern1, '(i8)') kontrl
         call xermsg ('slatec', 'xsetf',
     *      'invalid argument = ' // xern1, 1, 2)
         return
      endif
c
      junk = j4save(2,kontrl,.true.)
      return
      end
