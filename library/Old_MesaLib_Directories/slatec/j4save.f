*deck j4save
      function j4save (iwhich, ivalue, iset)
c***begin prologue  j4save
c***subsidiary
c***purpose  save or recall global variables needed by error
c            handling routines.
c***library   slatec (xerror)
c***type      integer (j4save-i)
c***keywords  error messages, error number, recall, save, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                = 1 refers to current error number.
c                = 2 refers to current error control flag.
c                = 3 refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                = 5 refers to the total number of units to which
c                     each error message is to be written.
c                = 6 refers to the 2nd unit for error messages
c                = 7 refers to the 3rd unit for error messages
c                = 8 refers to the 4th unit for error messages
c                = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c***see also  xermsg
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900205  minor modifications to prologue.  (wrb)
c   900402  added type section.  (wrb)
c   910411  added keywords section.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  j4save
      logical iset
      integer iparam(9)
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***first executable statement  j4save
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end
