*deck @(#)crowd.f	5.1  11/6/94
      subroutine crowd(natoms,dis,an,tooclo)
c***begin prologue     crowd.f
c***date written       870320  yymmdd
c***revision date      11/6/94
c***keywords           distance matrix, coordinates, crowd
c***author             binkley, et al., (g82)
c***source             @(#)crowd.f	5.1   11/6/94
c***purpose            checks the distance matrix for atoms too close together.
c***description
c     call crowd(natoms,dis,an,tooclo)
c       natoms ... number of atoms.
c       dis    ... distance matrix (natoms,natoms).
c       an     ... atomic numbers (natoms).
c       tooclo...  minimum acceptable distance (excludes dummies).
c
c***references
c***routines called    lnkerr(mdutil)
c***end prologue       crowd.f
      implicit none
c     --- input variables -----
      integer natoms
      real*8 tooclo
c     --- input arrays (unmodified) ---
      integer an(natoms)
      real*8 dis(natoms,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j
      character itoc*4
c
      do 20 i=1,natoms
         do 10 j=1,i-1
            if(dis(i,j).le.tooclo) then
               if(an(i).gt.0.and.an(j).gt.0)
     $            call lnkerr('atoms too close together:'
     $                        //itoc(i)//itoc(j))
               endif
   10    continue
   20 continue
c
c
      return
      end
