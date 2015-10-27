*deck @(#)lsubst.f	5.1  11/6/94
      function lsubst(namcnt,name,n)
c***begin prologue     lsubst.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           atomic number, atomic symbol
c***author             martin, richard (lanl)
c***source             @(#)lsubst.f	5.1   11/6/94
c***purpose            returns the sequential position of a given center in
c                      list.
c***description
c     lsubst is an integer function used as:
c     poscnt=lsubst(namcnt,name,n)
c       namcnt  the list of knownknown names.
c       name    the specific name for which to search.
c       n       the number of names in the array namcnt.
c
c     the sequential number of the center is returned in lsubst.
c***references
c***routines called    lnkerr(mdutil)
c***end prologue       lsubst.f
      implicit none
      integer lsubst
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      character*(*) namcnt(*),name*16
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
c     --- module to return the sequential number of the center whose
c         name is 'name'.  the list of n known names is in 'namcnt',
c         each of which occupy len(name) characters.
      lsubst=0
      do 10 i=1,n
         if(namcnt(i).eq.name) lsubst=i
   10 continue
      if(lsubst.eq.0)
     $   call lnkerr(' unknown center: '//name)
c
c
      return
      end
