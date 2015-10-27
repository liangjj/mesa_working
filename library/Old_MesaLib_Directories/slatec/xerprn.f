*deck xerprn
      subroutine xerprn (prefix, npref, messg, nwrap)
c***begin prologue  xerprn
c***subsidiary
c***purpose  print error messages processed by xermsg.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerprn-a)
c***keywords  error messages, printing, xerror
c***author  fong, kirby, (nmfecc at llnl)
c***description
c
c this routine sends one or more lines to each of the (up to five)
c logical units to which error messages are to be sent.  this routine
c is called several times by xermsg, sometimes with a single line to
c print and sometimes with a (potentially very long) message that may
c wrap around into multiple lines.
c
c prefix  input argument of type character.  this argument contains
c         characters to be put at the beginning of each line before
c         the body of the message.  no more than 16 characters of
c         prefix will be used.
c
c npref   input argument of type integer.  this argument is the number
c         of characters to use from prefix.  if it is negative, the
c         intrinsic function len is used to determine its length.  if
c         it is zero, prefix is not used.  if it exceeds 16 or if
c         len(prefix) exceeds 16, only the first 16 characters will be
c         used.  if npref is positive and the length of prefix is less
c         than npref, a copy of prefix extended with blanks to length
c         npref will be used.
c
c messg   input argument of type character.  this is the text of a
c         message to be printed.  if it is a long message, it will be
c         broken into pieces for printing on multiple lines.  each line
c         will start with the appropriate prefix and be followed by a
c         piece of the message.  nwrap is the number of characters per
c         piece; that is, after each nwrap characters, we break and
c         start a new line.  in addition the characters '$$' embedded
c         in messg are a sentinel for a new line.  the counting of
c         characters up to nwrap starts over for each new line.  the
c         value of nwrap typically used by xermsg is 72 since many
c         older error messages in the slatec library are laid out to
c         rely on wrap-around every 72 characters.
c
c nwrap   input argument of type integer.  this gives the maximum size
c         piece into which to break messg for printing on multiple
c         lines.  an embedded '$$' ends a line, and the count restarts
c         at the following character.  if a line break does not occur
c         on a blank (it would split a word) that word is moved to the
c         next line.  values of nwrap less than 16 will be treated as
c         16.  values of nwrap greater than 132 will be treated as 132.
c         the actual line length will be npref + nwrap after npref has
c         been adjusted to fall between 0 and 16 and nwrap has been
c         adjusted to fall between 16 and 132.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  i1mach, xgetua
c***revision history  (yymmdd)
c   880621  date written
c   880708  revised after the slatec cml subcommittee meeting of
c           june 29 and 30 to change the name to xerprn and to rework
c           the handling of the new line sentinel to behave like the
c           slash character in format statements.
c   890706  revised with the help of fred fritsch and reg clemens to
c           streamline the coding and fix a bug that caused extra blank
c           lines to be printed.
c   890721  revised to add a new feature.  a negative value of npref
c           causes len(prefix) to be used as the length.
c   891013  revised to correct error in calculating prefix length.
c   891214  prologue converted to version 4.0 format.  (wrb)
c   900510  added code to break messages between words.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerprn
      character*(*) prefix, messg
      integer npref, nwrap
      character*148 cbuff
      integer iu(5), nunit
      character*2 newlin
      parameter (newlin = '$$')
c***first executable statement  xerprn
      call xgetua(iu,nunit)
c
c       a zero value for a logical unit number means to use the standard
c       error message unit instead.  i1mach(4) retrieves the standard
c       error message unit.
c
      n = i1mach(4)
      do 10 i=1,nunit
         if (iu(i) .eq. 0) iu(i) = n
   10 continue
c
c       lpref is the length of the prefix.  the prefix is placed at the
c       beginning of cbuff, the character buffer, and kept there during
c       the rest of this routine.
c
      if ( npref .lt. 0 ) then
         lpref = len(prefix)
      else
         lpref = npref
      endif
      lpref = min(16, lpref)
      if (lpref .ne. 0) cbuff(1:lpref) = prefix
c
c       lwrap is the maximum number of characters we want to take at one
c       time from messg to print on one line.
c
      lwrap = max(16, min(132, nwrap))
c
c       set lenmsg to the length of messg, ignore any trailing blanks.
c
      lenmsg = len(messg)
      n = lenmsg
      do 20 i=1,n
         if (messg(lenmsg:lenmsg) .ne. ' ') go to 30
         lenmsg = lenmsg - 1
   20 continue
   30 continue
c
c       if the message is all blanks, then print one blank line.
c
      if (lenmsg .eq. 0) then
         cbuff(lpref+1:lpref+1) = ' '
         do 40 i=1,nunit
            write(iu(i), '(a)') cbuff(1:lpref+1)
   40    continue
         return
      endif
c
c       set nextc to the position in messg where the next substring
c       starts.  from this position we scan for the new line sentinel.
c       when nextc exceeds lenmsg, there is no more to print.
c       we loop back to label 50 until all pieces have been printed.
c
c       we look for the next occurrence of the new line sentinel.  the
c       index intrinsic function returns zero if there is no occurrence
c       or if the length of the first argument is less than the length
c       of the second argument.
c
c       there are several cases which should be checked for in the
c       following order.  we are attempting to set lpiece to the number
c       of characters that should be taken from messg starting at
c       position nextc.
c
c       lpiece .eq. 0   the new line sentinel does not occur in the
c                       remainder of the character string.  lpiece
c                       should be set to lwrap or lenmsg+1-nextc,
c                       whichever is less.
c
c       lpiece .eq. 1   the new line sentinel starts at messg(nextc:
c                       nextc).  lpiece is effectively zero, and we
c                       print nothing to avoid producing unnecessary
c                       blank lines.  this takes care of the situation
c                       where the library routine has a message of
c                       exactly 72 characters followed by a new line
c                       sentinel followed by more characters.  nextc
c                       should be incremented by 2.
c
c       lpiece .gt. lwrap+1  reduce lpiece to lwrap.
c
c       else            this last case means 2 .le. lpiece .le. lwrap+1
c                       reset lpiece = lpiece-1.  note that this
c                       properly handles the end case where lpiece .eq.
c                       lwrap+1.  that is, the sentinel falls exactly
c                       at the end of a line.
c
      nextc = 1
   50 lpiece = index(messg(nextc:lenmsg), newlin)
      if (lpiece .eq. 0) then
c
c       there was no new line sentinel found.
c
         idelta = 0
         lpiece = min(lwrap, lenmsg+1-nextc)
         if (lpiece .lt. lenmsg+1-nextc) then
            do 52 i=lpiece+1,2,-1
               if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
                  lpiece = i-1
                  idelta = 1
                  goto 54
               endif
   52       continue
         endif
   54    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      elseif (lpiece .eq. 1) then
c
c       we have a new line sentinel at messg(nextc:nextc+1).
c       don't print a blank line.
c
         nextc = nextc + 2
         go to 50
      elseif (lpiece .gt. lwrap+1) then
c
c       lpiece should be set down to lwrap.
c
         idelta = 0
         lpiece = lwrap
         do 56 i=lpiece+1,2,-1
            if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
               lpiece = i-1
               idelta = 1
               goto 58
            endif
   56    continue
   58    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      else
c
c       if we arrive here, it means 2 .le. lpiece .le. lwrap+1.
c       we should decrement lpiece by one.
c
         lpiece = lpiece - 1
         cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc  = nextc + lpiece + 2
      endif
c
c       print
c
      do 60 i=1,nunit
         write(iu(i), '(a)') cbuff(1:lpref+lpiece)
   60 continue
c
      if (nextc .le. lenmsg) go to 50
      return
      end
