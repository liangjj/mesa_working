*deck mc20ad
      subroutine mc20ad (nc, maxa, a, inum, jptr, jnum, jdisp)
c***begin prologue  mc20ad
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (mc20as-s, mc20ad-d)
c***author  (unknown)
c***description
c
c     this subprogram is a slight modification of a subprogram
c     from the c. 1979 aere harwell library.  the name of the
c     corresponding harwell code can be obtained by deleting
c     the final letter =d= in the names used here.
c     revised sep. 13, 1979.
c
c     royalties have been paid to aere-uk for use of their codes
c     in the package given here.  any primary usage of the harwell
c     subroutines requires a royalty agreement and payment between
c     the user and aere-uk.  any usage of the sandia written codes
c     dsplp( ) (which uses the harwell subroutines) is permitted.
c
c***see also  dsplp
c***routines called  (none)
c***revision history  (yymmdd)
c   811215  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  mc20ad
      integer inum(*), jnum(*)
      double precision a(*),ace,acep
      dimension jptr(nc)
c***first executable statement  mc20ad
      null = -jdisp
c**      clear jptr
      do 10 j=1,nc
         jptr(j) = 0
   10 continue
c**      count the number of elements in each column.
      do 20 k=1,maxa
         j = jnum(k) + jdisp
         jptr(j) = jptr(j) + 1
   20 continue
c**      set the jptr array
      k = 1
      do 30 j=1,nc
         kr = k + jptr(j)
         jptr(j) = k
         k = kr
   30 continue
c
c**      reorder the elements into column order.  the algorithm is an
c        in-place sort and is of order maxa.
      do 50 i=1,maxa
c        establish the current entry.
         jce = jnum(i) + jdisp
         if (jce.eq.0) go to 50
         ace = a(i)
         ice = inum(i)
c        clear the location vacated.
         jnum(i) = null
c        chain from current entry to store items.
         do 40 j=1,maxa
c        current entry not in correct position.  determine correct
c        position to store entry.
            loc = jptr(jce)
            jptr(jce) = jptr(jce) + 1
c        save contents of that location.
            acep = a(loc)
            icep = inum(loc)
            jcep = jnum(loc)
c        store current entry.
            a(loc) = ace
            inum(loc) = ice
            jnum(loc) = null
c        check if next current entry needs to be processed.
            if (jcep.eq.null) go to 50
c        it does.  copy into current entry.
            ace = acep
            ice = icep
            jce = jcep + jdisp
   40    continue
c
   50 continue
c
c**      reset jptr vector.
      ja = 1
      do 60 j=1,nc
         jb = jptr(j)
         jptr(j) = ja
         ja = jb
   60 continue
      return
      end
