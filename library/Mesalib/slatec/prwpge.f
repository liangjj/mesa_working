*deck prwpge
      subroutine prwpge (key, ipage, lpg, sx, ix)
c***begin prologue  prwpge
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (prwpge-s, dprwpg-d)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     prwpge limits the type of storage to a sequential scheme.
c     virtual memory page read/write subroutine.
c
c     depending on the value of key, subroutine prwpge() performs a page
c     read or write of page ipage. the page has length lpg.
c
c     key       is a flag indicating whether a page read or write is
c               to be performed.
c               if key = 1 data is read.
c               if key = 2 data is written.
c     ipage     is the page number of the matrix to be accessed.
c     lpg       is the length of the page of the matrix to be accessed.
c   sx(*),ix(*) is the matrix to be accessed.
c
c     this subroutine is a modification of the subroutine lrwpge,
c     sandia labs. rept. sand78-0785.
c     modifications by k.l. hiebert and r.j. hanson
c     revised 811130-1000
c     revised yymmdd-hhmm
c
c***see also  splp
c***routines called  prwvir, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  fixed error messages and replaced gotos with
c           if-then-else.  (rwc)
c   910403  updated author and description sections.  (wrb)
c***end prologue  prwpge
      real sx(*)
      dimension ix(*)
c***first executable statement  prwpge
c
c     check if ipage is in range.
c
      if (ipage.lt.1) then
         call xermsg ('slatec', 'prwpge',
     +      'the value of ipage (page number) was not in the range' //
     +      '1.le.ipage.le.maxpge.', 55, 1)
      endif
c
c     check if lpg is positive.
c
      if (lpg.le.0) then
         call xermsg ('slatec', 'prwpge',
     +      'the value of lpg (page length) was nonpositive.', 55, 1)
      endif
c
c     decide if we are reading or writing.
c
      if (key.eq.1) then
c
c        code to do a page read.
c
         call prwvir(key,ipage,lpg,sx,ix)
      else if (key.eq.2) then
c
c        code to do a page write.
c
         call prwvir(key,ipage,lpg,sx,ix)
      else
         call xermsg ('slatec', 'prwpge',
     +      'the value of key (read-write flag) was not 1 or 2.', 55, 1)
      endif
      return
      end
