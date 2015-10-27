*deck iploc
      integer function iploc (loc, sx, ix)
c***begin prologue  iploc
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (iploc-s, idloc-d)
c***keywords  relative address determination function, slatec
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c   given a "virtual" location,  iploc returns the relative working
c   address of the vector component stored in sx, ix.  any necessary
c   page swaps are performed automatically for the user in this
c   function subprogram.
c
c   loc       is the "virtual" address of the data to be retrieved.
c   sx ,ix    represent the matrix where the data is stored.
c
c***see also  splp
c***routines called  prwpge, xermsg
c***revision history  (yymmdd)
c   810306  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890606  restructured to match double precision version.  (wrb)
c   890606  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   910731  added code to set iploc to 0 if loc is non-positive.  (wrb)
c***end prologue  iploc
      real sx(*)
      integer ix(*)
c***first executable statement  iploc
      if (loc.le.0) then
         call xermsg ('slatec', 'iploc',
     +     'a value of loc, the first argument, .le. 0 was encountered',
     +     55, 1)
         iploc = 0
         return
      endif
c
c     two cases exist:  (1.le.loc.le.k) .or. (loc.gt.k).
c
      k = ix(3) + 4
      lmx = ix(1)
      lmxm1 = lmx - 1
      if (loc.le.k) then
         iploc = loc
         return
      endif
c
c     compute length of the page, starting address of the page, page
c     number and relative working address.
c
      lpg = lmx-k
      itemp = loc - k - 1
      ipage = itemp/lpg + 1
      iploc = mod(itemp,lpg) + k + 1
      np = abs(ix(lmxm1))
c
c     determine if a page fault has occurred.  if so, write page np
c     and read page ipage.  write the page only if it has been
c     modified.
c
      if (ipage.ne.np) then
         if (sx(lmx).eq.1.0) then
            sx(lmx) = 0.0
            key = 2
            call prwpge (key, np, lpg, sx, ix)
         endif
         key = 1
         call prwpge (key, ipage, lpg, sx, ix)
      endif
      return
      end
