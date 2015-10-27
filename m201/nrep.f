*deck @(#)nrep.f	5.1  11/6/94
      function nrep(item,list,length)
c***begin prologue     nrep
c***date written       850601  (yymmdd)
c***revision date      11/6/94 
c***keywords           repetitions, vector
c***author             martin, richard (lanl)
c***source             @(#)nrep.f	5.1 11/6/94
c***purpose            finds the number of repetitions of an item
c                      in an array.
c***description
c
c     nrep is an integer function used as:
c                      number=nrep(item,list,length)
c                        item     item to be counted in the list.
c                        list     input vector; integer(length).
c                        length   length of array.
c
c     return as the value of the function, nrep, the number of times
c     that item appears in the list.  the first length elements of list
c     are to be searched.  note that for a match only the absolute value
c     of item need be equal to the absolute of an element of the list.
c
c***references
c***routines called    
c***end prologue       nrep.f
      implicit none
c     --- input variables -----
      integer nrep
      integer item,length
c     --- input arrays (unmodified) ---
      integer list(length)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ltest,i
      integer abs
c
c
      ltest=abs(item)
      nrep=0
      if(length.gt.0) then
         do 10 i=1,length
            if(ltest.eq.abs(list(i))) nrep=nrep+1
   10    continue
      endif
c
c
      return
      end
