*deck xsetua
      subroutine xsetua (iunita, n)
c***begin prologue  xsetua
c***purpose  set logical unit numbers (up to 5) to which error
c            messages are to be sent.
c***library   slatec (xerror)
c***category  r3b
c***type      all (xsetua-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xsetua may be called to declare a list of up to five
c        logical units, each of which is to receive a copy of
c        each error message processed by this package.
c        the purpose of xsetua is to allow simultaneous printing
c        of each error message on, say, a main output file,
c        an interactive terminal, and other files such as graphics
c        communication files.
c
c     description of parameters
c      --input--
c        iunit - an array of up to five unit numbers.
c                normally these numbers should all be different
c                (but duplicates are not prohibited.)
c        n     - the number of unit numbers provided in iunit
c                must have 1 .le. n .le. 5.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save, xermsg
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  change call to xerrwv to xermsg.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xsetua
      dimension iunita(5)
      character *8 xern1
c***first executable statement  xsetua
c
      if (n.lt.1 .or. n.gt.5) then
         write (xern1, '(i8)') n
         call xermsg ('slatec', 'xsetua',
     *      'invalid number of units, n = ' // xern1, 1, 2)
         return
      endif
c
      do 10 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         junk = j4save(index,iunita(i),.true.)
   10 continue
      junk = j4save(5,n,.true.)
      return
      end
