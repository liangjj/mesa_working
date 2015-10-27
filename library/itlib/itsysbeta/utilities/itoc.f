*deck @(#)itoc.f	1.1  11/30/90
      function itoc(num)
c***begin prologue     itoc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, integer, conversion
c***author             martin, richard (lanl)
c***source             @(#)itoc.f	1.1   11/30/90
c***purpose            converts an integer into a string.
c***description
c                      itoc is a character function used as:
c                        string=itoc(num)
c                          num  the integer to convert.
c
c                      itoc will work for integers between 10**-15 and 10**16.
c
c***references
c***routines called    (none)
c***end prologue       itoc
      implicit integer(a-z)
      character*(*) itoc
      character digits*10,k*1,str*16
      data digits/'0123456789'/
c
c
      maxsiz=len(itoc)
      str=' '
      n=iabs(num)
      i=0
c
c     generate digits.
   10 i=i+1
         d=mod(n,10)
         str(i:i)=digits(d+1:d+1)
         n=n/10
         if(n.gt.0.and.i.lt.maxsiz) goto 10
c
c     generate sign.
      if(num.lt.0.and.i.lt.maxsiz) then
         i=i+1
         str(i:i)='-'
      endif
c
c     now flip it around.
      halfi=i/2
      do 20 j=1,halfi
         k=str(j:j)
         str(j:j)=str(i+1-j:i+1-j)
         str(i+1-j:i+1-j)=k
   20 continue
c
c
      itoc=str
c
c
      return
      end
