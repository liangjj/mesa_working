*deck @(#)posfil.f	2.1  10/10/91
      function posfil(key,card,file)
c***begin prologue     posfil
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           position, input file
c***author             martin, richard (lanl)
c***source
c***purpose            module to position the input file to the requested key.
c***description
c                      this routine is a logical function which returns .true.
c                      if the "key" is found in the input file.
c                      it rewinds the input file before the search, and
c                      leaves it positioned to the card directly after the
c                      keyword.
c
c                      posfil is a logical function used as:
c                      test=posfil(key,card,file)
c                        key      character string for which to search.
c                        card     character string which acts as a buffer for
c                                 each card read in the input file.  should
c                                 be dimensioned 80 characters in the calling
c                                 routine (or whatever).
c                        file     file to search.
c
c
c***references
c***routines called    captlz(chr), lnkerr(mdutil)
c***end prologue       posfil
      implicit integer(a-z)
      logical posfil
      character*(*) key,card
c
 1000 format(a)
c
c
      posfil=.true.
c
      rewind file
   10 read(file,1000,end=9000) card
         call locase(card,card)
         if(index(card,key).ne.0) return
      goto 10
c
c
 9000 posfil=.false.
c
c
      return
      end
