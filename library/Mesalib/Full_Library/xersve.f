*deck xersve
      subroutine xersve (librar, subrou, messg, kflag, nerr, level,
     +   icount)
c***begin prologue  xersve
c***subsidiary
c***purpose  record that an error has occurred.
c***library   slatec (xerror)
c***category  r3
c***type      all (xersve-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c *usage:
c
c        integer  kflag, nerr, level, icount
c        character * (len) librar, subrou, messg
c
c        call xersve (librar, subrou, messg, kflag, nerr, level, icount)
c
c *arguments:
c
c        librar :in    is the library that the message is from.
c        subrou :in    is the subroutine that the message is from.
c        messg  :in    is the message to be saved.
c        kflag  :in    indicates the action to be performed.
c                      when kflag > 0, the message in messg is saved.
c                      when kflag=0 the tables will be dumped and
c                      cleared.
c                      when kflag < 0, the tables will be dumped and
c                      not cleared.
c        nerr   :in    is the error number.
c        level  :in    is the error severity.
c        icount :out   the number of times this message has been seen,
c                      or zero if the table has overflowed and does not
c                      contain this message specifically.  when kflag=0,
c                      icount will not be altered.
c
c *description:
c
c   record that this error occurred and possibly dump and clear the
c   tables.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called
c***revision history  (yymmdd)
c   800319  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900413  routine modified to remove reference to kflag.  (wrb)
c   900510  changed to add library name and subroutine to calling
c           sequence, use if-then-else, make number of saved entries
c           easily changeable, changed routine name from xersav to
c           xersve.  (rwc)
c   910626  added libtab and subtab to save statement.  (bks)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xersve
      parameter (lentab=10)
      character*(*) librar, subrou, messg
      character*8  libtab(lentab), subtab(lentab), lib, sub
      character*20 mestab(lentab), mes
      dimension nertab(lentab), levtab(lentab), kount(lentab)
      common/io/inp, iout
      save libtab, subtab, mestab, nertab, levtab, kount, kountx, nmsg
      data kountx/0/, nmsg/0/
c***first executable statement  xersve
c
      iu=iout
      if (kflag.le.0) then
c
c        dump the table.
c
         if (nmsg.eq.0) return
c
c        print to each unit.
c
c           print the table header.
c
            write (iu,9000)
c
c           print body of table.
c
            do 10 i = 1,nmsg
               write (iu,9010) libtab(i), subtab(i), mestab(i),
     *         nertab(i),levtab(i),kount(i)
   10       continue
c
c           print number of other errors.
c
            if (kountx.ne.0) write (iu,9020) kountx
            write (iu,9030)
c
c        clear the error tables.
c
         if (kflag.eq.0) then
            nmsg = 0
            kountx = 0
         endif
      else
c
c        process a message...
c        search for this messg, or else an empty slot for this messg,
c        or else determine that the error table is full.
c
         lib = librar
         sub = subrou
         mes = messg
         do 30 i = 1,nmsg
            if (lib.eq.libtab(i) .and. sub.eq.subtab(i) .and.
     *         mes.eq.mestab(i) .and. nerr.eq.nertab(i) .and.
     *         level.eq.levtab(i)) then
                  kount(i) = kount(i) + 1
                  icount = kount(i)
                  return
            endif
   30    continue
c
         if (nmsg.lt.lentab) then
c
c           empty slot found for new message.
c
            nmsg = nmsg + 1
            libtab(i) = lib
            subtab(i) = sub
            mestab(i) = mes
            nertab(i) = nerr
            levtab(i) = level
            kount (i) = 1
            icount    = 1
         else
c
c           table is full.
c
            kountx = kountx+1
            icount = 0
         endif
      endif
      return
c
c     formats.
c
 9000 format ('0          error message summary' /
     +   ' library    subroutine message start             nerr',
     +   '     level     count')
 9010 format (1x,a,3x,a,3x,a,3i10)
 9020 format ('0other errors not individually tabulated = ', i10)
 9030 format (1x)
      end
