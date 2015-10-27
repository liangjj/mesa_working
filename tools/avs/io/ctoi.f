*deck @(#)ctoi.f	4.1  7/7/93
      function ctoi(c)
c***begin prologue     ctoi
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, integer, conversion
c***author             saxe, paul (lanl)
c***source             @(#)ctoi.f	4.1   7/7/93
c***purpose            converts a character string into an integer.
c***description
c                      ctoi is an integer function used as:
c                        int=ctoi(c)
c                          c  the input character string.
c
c                      the first unrecognized character other than the sign,
c                      i.e. the first non-digit, terminates the conversion.
c
c***references
c***routines called    (none)
c***end prologue       ctoi
c
      implicit integer (a-z)
      integer ctoi
c
      character*(*) c,digits*10
c
      data digits /'0123456789'/
      save digits
c
c
      temp=0
      if (c(1:1).eq.'-'.or.c(1:1).eq.'+') then
         mn=2
      else
         mn=1
      end if
c
c     note that the first unrecognized character kicks us out
c     of the loop.
      do 1 i=mn,len(c)
         intgr=index(digits,c(i:i))-1
         if(intgr.lt.0) goto 2
         temp=10*temp+intgr
    1 continue
    2 continue
c
      if (c(1:1).eq.'-') temp=-temp
      ctoi=temp
c
c
      return
      end
