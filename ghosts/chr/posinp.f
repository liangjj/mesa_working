*deck @(#)posinp.f	2.1  10/10/91
      function posinp(key,card)
c***begin prologue     posinp
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
c                      posinp is a logical function used as:
c                      test=posinp(key,card)
c                        key      character string for which to search.
c                        card     character string which acts as a buffer for
c                                 each card read in the input file.  should
c                                 be dimensioned 80 characters in the calling
c                                 routine (or whatever).
c
c
c***references
c***routines called    captlz(chr), lnkerr(mdutil)
c***end prologue       posinp
      implicit integer(a-z)
      logical posinp
      character*(*) key,card
c
      common/io/inp,iout
c
 1000 format(a)
c
c
      posinp=.true.
c
      rewind inp
   10 read(inp,1000,end=9000) card
         call locase(card,card)
         if(index(card,key).ne.0) return
      goto 10
c
c
 9000 posinp=.false.
c
c
      return
      end
