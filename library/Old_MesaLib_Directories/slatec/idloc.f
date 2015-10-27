*deck idloc
      integer function idloc (loc, sx, ix)
c***begin prologue  idloc
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (iploc-s, idloc-d)
c***keywords  relative address determination function, slatec
c***author  boland, w. robert, (lanl)
c           nicol, tom, (university of british columbia)
c***description
c
c   given a "virtual" location,  idloc returns the relative working
c   address of the vector component stored in sx, ix.  any necessary
c   page swaps are performed automatically for the user in this
c   function subprogram.
c
c   loc       is the "virtual" address of the data to be retrieved.
c   sx ,ix    represent the matrix where the data is stored.
c
c***see also  dsplp
c***routines called  dprwpg, xermsg
c***revision history  (yymmdd)
c   890606  date written
c   890606  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   910731  added code to set idloc to 0 if loc is non-positive.  (wrb)
c***end prologue  idloc
      double precision sx(*)
      integer ix(*)
c***first executable statement  idloc
      if (loc.le.0) then
         call xermsg ('slatec', 'idloc',
     +     'a value of loc, the first argument, .le. 0 was encountered',
     +     55, 1)
         idloc = 0
         return
      endif
c
c     two cases exist:  (1.le.loc.le.k) .or. (loc.gt.k).
c
      k = ix(3) + 4
      lmx = ix(1)
      lmxm1 = lmx - 1
      if (loc.le.k) then
         idloc = loc
         return
      endif
c
c     compute length of the page, starting address of the page, page
c     number and relative working address.
c
      lpg = lmx-k
      itemp = loc - k - 1
      ipage = itemp/lpg + 1
      idloc = mod(itemp,lpg) + k + 1
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
            call dprwpg (key, np, lpg, sx, ix)
         endif
         key = 1
         call dprwpg (key, ipage, lpg, sx, ix)
      endif
      return
      end
