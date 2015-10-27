*deck xermsg
      subroutine xermsg (librar, subrou, messg, nerr, level)
c***begin prologue  xermsg
c***purpose  process error messages for slatec and other libraries.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xermsg-a)
c***keywords  error message, xerror
c***author  fong, kirby, (nmfecc at llnl)
c***description
c
c   xermsg processes a diagnostic message in a manner determined by the
c   value of level and the current value of the library error control
c   flag, kontrl.  see subroutine xsetf for details.
c
c    librar   a character constant (or character variable) with the name
c             of the library.  this will be 'slatec' for the slatec
c             common math library.  the error handling package is
c             general enough to be used by many libraries
c             simultaneously, so it is desirable for the routine that
c             detects and reports an error to identify the library name
c             as well as the routine name.
c
c    subrou   a character constant (or character variable) with the name
c             of the routine that detected the error.  usually it is the
c             name of the routine that is calling xermsg.  there are
c             some instances where a user callable library routine calls
c             lower level subsidiary routines where the error is
c             detected.  in such cases it may be more informative to
c             supply the name of the routine the user called rather than
c             the name of the subsidiary routine that detected the
c             error.
c
c    messg    a character constant (or character variable) with the text
c             of the error or warning message.  in the example below,
c             the message is a character constant that contains a
c             generic message.
c
c                   call xermsg ('slatec', 'mmpy',
c                  *'the order of the matrix exceeds the row dimension',
c                  *3, 1)
c
c             it is possible (and is sometimes desirable) to generate a
c             specific message--e.g., one that contains actual numeric
c             values.  specific numeric values can be converted into
c             character strings using formatted write statements into
c             character variables.  this is called standard fortran
c             internal file i/o and is exemplified in the first three
c             lines of the following example.  you can also catenate
c             substrings of characters to construct the error message.
c             here is an example showing the use of both writing to
c             an internal file and catenating character strings.
c
c                   character*5 charn, charl
c                   write (charn,10) n
c                   write (charl,10) lda
c                10 format(i5)
c                   call xermsg ('slatec', 'mmpy', 'the order'//charn//
c                  *   ' of the matrix exceeds its row dimension of'//
c                  *   charl, 3, 1)
c
c             there are two subtleties worth mentioning.  one is that
c             the // for character catenation is used to construct the
c             error message so that no single character constant is
c             continued to the next line.  this avoids confusion as to
c             whether there are trailing blanks at the end of the line.
c             the second is that by catenating the parts of the message
c             as an actual argument rather than encoding the entire
c             message into one large character variable, we avoid
c             having to know how long the message will be in order to
c             declare an adequate length for that large character
c             variable.  xermsg calls xerprn to print the message using
c             multiple lines if necessary.  if the message is very long,
c             xerprn will break it into pieces of 72 characters (as
c             requested by xermsg) for printing on multiple lines.
c             also, xermsg asks xerprn to prefix each line with ' *  '
c             so that the total line length could be 76 characters.
c             note also that xerprn scans the error message backwards
c             to ignore trailing blanks.  another feature is that
c             the substring '$$' is treated as a new line sentinel
c             by xerprn.  if you want to construct a multiline
c             message without having to count out multiples of 72
c             characters, just use '$$' as a separator.  '$$'
c             obviously must occur within 72 characters of the
c             start of each line to have its intended effect since
c             xerprn is asked to wrap around at 72 characters in
c             addition to looking for '$$'.
c
c    nerr     an integer value that is chosen by the library routine's
c             author.  it must be in the range -99 to 999 (three
c             printable digits).  each distinct error should have its
c             own error number.  these error numbers should be described
c             in the machine readable documentation for the routine.
c             the error numbers need be unique only within each routine,
c             so it is reasonable for each routine to start enumerating
c             errors from 1 and proceeding to the next integer.
c
c    level    an integer value in the range 0 to 2 that indicates the
c             level (severity) of the error.  their meanings are
c
c            -1  a warning message.  this is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.  an attempt is made to only print this
c                message once.
c
c             0  a warning message.  this is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.
c
c             1  a recoverable error.  this is used even if the error is
c                so serious that the routine cannot return any useful
c                answer.  if the user has told the error package to
c                return after recoverable errors, then xermsg will
c                return to the library routine which can then return to
c                the user's routine.  the user may also permit the error
c                package to terminate the program upon encountering a
c                recoverable error.
c
c             2  a fatal error.  xermsg will not return to its caller
c                after it receives a fatal error.  this level should
c                hardly ever be used; it is much better to allow the
c                user a chance to recover.  an example of one of the few
c                cases in which it is permissible to declare a level 2
c                error is a reverse communication library routine that
c                is likely to be called repeatedly until it integrates
c                across some interval.  if there is a serious error in
c                the input such that another step cannot be taken and
c                the library routine is called again without the input
c                error having been corrected by the caller, the library
c                routine will probably be called forever with improper
c                input.  in this case, it is reasonable to declare the
c                error to be fatal.
c
c    each of the arguments to xermsg is input; none will be modified by
c    xermsg.  a routine may make multiple calls to xermsg with warning
c    level messages; however, after a call to xermsg with a recoverable
c    error, the routine should return to the user.  do not try to call
c    xermsg with a second recoverable error after the first recoverable
c    error because the error package saves the error number.  the user
c    can retrieve this error number by calling another entry point in
c    the error handling package and then clear the error number when
c    recovering from the error.  calling xermsg in succession causes the
c    old error number to be overwritten by the latest error number.
c    this is considered harmless for error numbers associated with
c    warning messages but must not be done for error numbers of serious
c    errors.  after a call to xermsg with a recoverable error, the user
c    must be given a chance to call numxer or xerclr to retrieve or
c    clear the error number.
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  fdump, j4save, xercnt, xerhlt, xerprn, xersve
c***revision history  (yymmdd)
c   880101  date written
c   880621  revised as directed at slatec cml meeting of february 1988.
c           there are two basic changes.
c           1.  a new routine, xerprn, is used instead of xerprt to
c               print messages.  this routine will break long messages
c               into pieces for printing on multiple lines.  '$$' is
c               accepted as a new line sentinel.  a prefix can be
c               added to each line to be printed.  xermsg uses either
c               ' ***' or ' *  ' and long messages are broken every
c               72 characters (at most) so that the maximum line
c               length output can now be as great as 76.
c           2.  the text of all messages is now in upper case since the
c               fortran standard document does not admit the existence
c               of lower case.
c   880708  revised after the slatec cml meeting of june 29 and 30.
c           the principal changes are
c           1.  clarify comments in the prologues
c           2.  rename xrprnt to xerprn
c           3.  rework handling of '$$' in xerprn to handle blank lines
c               similar to the way format statements handle the /
c               character for new records.
c   890706  revised with the help of fred fritsch and reg clemens to
c           clean up the coding.
c   890721  revised to use new feature in xerprn to count characters in
c           prefix.
c   891013  revised to correct comments.
c   891214  prologue converted to version 4.0 format.  (wrb)
c   900510  changed test on nerr to be -9999999 < nerr < 99999999, but
c           nerr .ne. 0, and on level to be -2 < level < 3.  added
c           level=-1 logic, changed calls to xersav to xersve, and
c           xerctl to xercnt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xermsg
      character*(*) librar, subrou, messg
      character*8 xlibr, xsubr
      character*72  temp
      character*20  lfirst
c***first executable statement  xermsg
      lkntrl = j4save (2, 0, .false.)
      maxmes = j4save (4, 0, .false.)
c
c       lkntrl is a local copy of the control flag kontrl.
c       maxmes is the maximum number of times any particular message
c          should be printed.
c
c       we print a fatal error message and terminate for an error in
c          calling xermsg.  the error number should be positive,
c          and the level should be between 0 and 2.
c
      if (nerr.lt.-9999999 .or. nerr.gt.99999999 .or. nerr.eq.0 .or.
     *   level.lt.-1 .or. level.gt.2) then
         call xerprn (' ***', -1, 'fatal error in...$$ ' //
     *      'xermsg -- invalid error number or level$$ '//
     *      'job abort due to fatal error.', 72)
         call xersve (' ', ' ', ' ', 0, 0, 0, kdummy)
         call xerhlt (' ***xermsg -- invalid input')
         return
      endif
c
c       record the message.
c
      i = j4save (1, nerr, .true.)
      call xersve (librar, subrou, messg, 1, nerr, level, kount)
c
c       handle print-once warning messages.
c
      if (level.eq.-1 .and. kount.gt.1) return
c
c       allow temporary user override of the control flag.
c
      xlibr  = librar
      xsubr  = subrou
      lfirst = messg
      lerr   = nerr
      llevel = level
      call xercnt (xlibr, xsubr, lfirst, lerr, llevel, lkntrl)
c
      lkntrl = max(-2, min(2,lkntrl))
      mkntrl = abs(lkntrl)
c
c       skip printing if the control flag value as reset in xercnt is
c       zero and the error is not fatal.
c
      if (level.lt.2 .and. lkntrl.eq.0) go to 30
      if (level.eq.0 .and. kount.gt.maxmes) go to 30
      if (level.eq.1 .and. kount.gt.maxmes .and. mkntrl.eq.1) go to 30
      if (level.eq.2 .and. kount.gt.max(1,maxmes)) go to 30
c
c       announce the names of the library and subroutine by building a
c       message in character variable temp (not exceeding 66 characters)
c       and sending it out via xerprn.  print only if control flag
c       is not zero.
c
      if (lkntrl .ne. 0) then
         temp(1:21) = 'message from routine '
         i = min(len(subrou), 16)
         temp(22:21+i) = subrou(1:i)
         temp(22+i:33+i) = ' in library '
         ltemp = 33 + i
         i = min(len(librar), 16)
         temp(ltemp+1:ltemp+i) = librar (1:i)
         temp(ltemp+i+1:ltemp+i+1) = '.'
         ltemp = ltemp + i + 1
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       if lkntrl is positive, print an introductory line before
c       printing the message.  the introductory line tells the choice
c       from each of the following three options.
c       1.  level of the message
c              'informative message'
c              'potentially recoverable error'
c              'fatal error'
c       2.  whether control flag will allow program to continue
c              'prog continues'
c              'prog aborted'
c       3.  whether or not a traceback was requested.  (the traceback
c           may not be implemented at some sites, so this only tells
c           what was requested, not what was delivered.)
c              'traceback requested'
c              'traceback not requested'
c       notice that the line including four prefix characters will not
c       exceed 74 characters.
c       we skip the next block if the introductory line is not needed.
c
      if (lkntrl .gt. 0) then
c
c       the first part of the message tells about the level.
c
         if (level .le. 0) then
            temp(1:20) = 'informative message,'
            ltemp = 20
         elseif (level .eq. 1) then
            temp(1:30) = 'potentially recoverable error,'
            ltemp = 30
         else
            temp(1:12) = 'fatal error,'
            ltemp = 12
         endif
c
c       then whether the program will continue.
c
         if ((mkntrl.eq.2 .and. level.ge.1) .or.
     *       (mkntrl.eq.1 .and. level.eq.2)) then
            temp(ltemp+1:ltemp+14) = ' prog aborted,'
            ltemp = ltemp + 14
         else
            temp(ltemp+1:ltemp+16) = ' prog continues,'
            ltemp = ltemp + 16
         endif
c
c       finally tell whether there should be a traceback.
c
         if (lkntrl .gt. 0) then
            temp(ltemp+1:ltemp+20) = ' traceback requested'
            ltemp = ltemp + 20
         else
            temp(ltemp+1:ltemp+24) = ' traceback not requested'
            ltemp = ltemp + 24
         endif
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       now send out the message.
c
      call xerprn (' *  ', -1, messg, 72)
c
c       if lkntrl is positive, write the error number and request a
c          traceback.
c
      if (lkntrl .gt. 0) then
         write (temp, '(''error number = '', i8)') nerr
         do 10 i=16,22
            if (temp(i:i) .ne. ' ') go to 20
   10    continue
c
   20    call xerprn (' *  ', -1, temp(1:15) // temp(i:23), 72)
         call fdump
      endif
c
c       if lkntrl is not zero, print a blank line and an end of message.
c
      if (lkntrl .ne. 0) then
         call xerprn (' *  ', -1, ' ', 72)
         call xerprn (' ***', -1, 'end of message', 72)
         call xerprn ('    ',  0, ' ', 72)
      endif
c
c       if the error is not fatal or the error is recoverable and the
c       control flag is set for recovery, then return.
c
   30 if (level.le.0 .or. (level.eq.1 .and. mkntrl.le.1)) return
c
c       the program will be stopped due to an unrecovered error or a
c       fatal error.  print the reason for the abort and the error
c       summary if the control flag and the maximum error count permit.
c
      if (lkntrl.gt.0 .and. kount.lt.max(1,maxmes)) then
         if (level .eq. 1) then
            call xerprn
     *         (' ***', -1, 'job abort due to unrecovered error.', 72)
         else
            call xerprn(' ***', -1, 'job abort due to fatal error.', 72)
         endif
         call xersve (' ', ' ', ' ', -1, 0, 0, kdummy)
         call xerhlt (' ')
      else
         call xerhlt (messg)
      endif
      return
      end
