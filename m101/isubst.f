*deck @(#)isubst.f	5.1  11/6/94
      function isubst(char)
c***begin prologue     isubst.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           atomic number, atomic symbol
c***author             martin, richard (lanl)
c***source             @(#)isubst.f	5.1   11/6/94
c***purpose            maps an atomic symbol onto an atomic number.
c***description
c     isubst is an integer function used as:
c     atomno=isubst(char)
c       char    character string containing the atomic symbol.
c               the atomic number is returned as the value of isubst.
c
c***references
c***routines called    streqc(chr), fillel(util), lnkerr(mdutil)
c***end prologue       isubst.f
      implicit none
      integer isubst
c     --- input variables -----
      character*(*) char
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer ni,i
      character*2 keyi(106), atnam
      logical streqc, called
c
      data ni/106/, called/.false./
      save ni,keyi,called
c
      common/io/inp,iout
c
 1000 format(' unrecognized atomic symbol: ',80a1)
c
c     --- a routine to return the atomic number of a center.  the name
c         of the center is passed as a string in "char".  the
c         returned value of this function is the atomic number.
      isubst = -1
      if(streqc('-',char(1:1))) return
      if(.not.called) then
         called = .true.
         call fillel(-1,ni-2,keyi)
      endif
c
c     --- see if the first two characters match an element name.
      atnam=char(1:2)
      do 10 i=1,ni
         isubst=i-2
c
         if(streqc(atnam,keyi(i))) return
c
   10 continue
c
c
c     --- the first two characters in the center name don't match any 
c     known atomic symbol.  it may be a single character symbol followed
c     by an identifier, e.g. n1.  see if the first character matches.
      do 20 i=1,ni
         isubst=i-2
c
         if(streqc(atnam(1:1),keyi(i)(1:1))) then
c
c           --- the first character matches.  make sure that the second
c           character in the element name is a blank.  this is to get
c           around situations in which b1 is confused with be.
c
            if(keyi(i)(2:2).eq.' ') return
c
         endif
   20 continue
c
c     --- unrecognized symbol.
      write(iout,1000) (atnam(i:i),i=1,2)
      call lnkerr(' ')
c
c
      return
      end
