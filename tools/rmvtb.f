      program rmvtb
      implicit none
      character*80 card
      character*80 newline
      integer cskipb,pos
      integer stdin,stdout,i
c
c     this little code strips out trailing blanks from each
c     line of the standard input, and sends it bakc out to standard
c     ouput. completely blank lines are left in the output.
c
      data stdin /1/
      data stdout/2/
c
 1000 format(a80)
 1010 format(80a1)
c
    5 continue
      do 10 i=1,80
         newline(i:i)=' '
   10 continue
      read(*,1000,end=20) card
      pos=cskipb(card,' ') 
      newline(1:pos)=card(1:pos)
      write(*,1010) (newline(i:i),i=1,pos)
      go to 5
c
c
   20 continue
      stop
      end
*deck @(#)cskipb.f	4.1  7/7/93
      function cskipb(string,substr)
c***begin prologue     cskipb
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, string, skip
c***author             martin, richard (lanl)
c***source             @(#)cskipb.f	4.1   7/7/93
c***purpose            skips backwards over all occurrences of a pattern.
c***description
c                      cskipb is an integer function used as:
c                        pos=cskipb(string,substr)
c                          string  input character string to search.
c                          substr  string to skip over.
c
c                      cskipb returns the index of the last character
c                      which does not match the pattern.  it returns 0
c                      if the string contains only the pattern of interest.
c                      it is generally used to return the length of a string
c                      ignoring trailing blanks, i.e; length=cksipb(str,' ')
c***references
c***routines called    (none)
c***end prologue       cskipb
      implicit integer(a-z)
      integer cskipb
      character*(*) string,substr
c
c
      lenstr=len(string)
      lensub=len(substr)
      cskipb=0
      if(lensub.gt.lenstr) return
c
      do 10 i=lenstr,1,-lensub
         if(string(i+1-lensub:i).ne.substr) then
            cskipb=i
            return
         endif
   10 continue
c
c
      return
      end
