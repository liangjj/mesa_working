!function itoc.f90	5.1  11/6/94
      function itoc(num)
!***begin prologue     itoc
!***date written       850601  (yymmdd)
!***revision date      yymmdd  (yymmdd)
!***keywords           character, integer, conversion
!***author             martin, richard (lanl)
!***source             @(#)itoc.f	5.1   11/6/94
!***purpose            converts an integer into a string.
!***description
!                      itoc is a character function used as:
!                      string=itoc(num)
!                      num  the integer to convert.
!
!                      itoc will work for integers between 10**-15 and 10**16.
!
!***references
!***routines called    (none)
!***end prologue       itoc
      implicit integer(a-z)
      character*(*) itoc
      character digits*10,k*1,str*16
      data digits/'0123456789'/
      save digits
      maxsiz=len(itoc)
      str=' '
      n=iabs(num)
      i=0
!     generate digits.
   10 i=i+1
         d=mod(n,10)
         str(i:i)=digits(d+1:d+1)
         n=n/10
         if(n.gt.0.and.i.lt.maxsiz) goto 10
!     generate sign.
      if(num.lt.0.and.i.lt.maxsiz) then
         i=i+1
         str(i:i)='-'
      endif
!
!     now flip it around.
      halfi=i/2
      do 20 j=1,halfi
         k=str(j:j)
         str(j:j)=str(i+1-j:i+1-j)
         str(i+1-j:i+1-j)=k
   20 continue
      itoc=str
      return
      end
