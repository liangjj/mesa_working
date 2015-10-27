*deck @(#)crowd.f	1.1  11/30/90
      subroutine crowd(natoms,dis,an,tooclo)
c***begin prologue     crowd
c***date written       870320  yymmdd
c***revision date      870320  yymmdd
c***keywords           distance matrix, coordinates, crowd
c***author             gauss82
c***source             @(#)crowd.f	1.1   11/30/90
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
c***end prologue       crowd
      implicit integer(a-z)
      integer an(natoms)
      real*8 dis(natoms,natoms),tooclo
      character itoc*4
      common/io/inp, iout
c
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
