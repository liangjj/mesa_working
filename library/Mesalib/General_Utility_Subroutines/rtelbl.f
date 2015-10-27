*deck @(#)rtelbl.f	5.1  11/6/94
      function rtelbl(ncards,route,label)
c***begin prologue     rtelbl
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           route, label
c***author             martin, richard (lanl)
c***source
c***purpose            searches for a label in the route.
c***description
c     rtelbl is an integer function used as:
c       card=rtelbl(ncards,route,label)
c          ncards  the number of cards in the route.
c          route   the nonstandard route.
c          label   the label for which to search.
c
c     if the label is not found, rtelbl=ncards+1
c***references
c***routines called    (none)
c***end prologue       rtelbl
      implicit integer(a-z)
      integer rtelbl
      character*(*) route(ncards),label
c
c
      rtelbl=ncards+1
      lenlbl=len(label)
      do 10 i=1,ncards
         rtstrt=cskipf(route(i),' ')
         rtend=rtstrt+lenlbl-1
         if(index(route(i)(rtstrt:rtend),label).ne.0) then
            rtelbl=i
            return
         endif
   10 continue
c
c
      return
      end
