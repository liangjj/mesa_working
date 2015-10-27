*deck @(#)atmset.f	5.1  11/6/94
      subroutine atmset(symat,atprmt,natoms,nop,nset,maxset)
c
c***begin prologue     atmset
c***date written       870715   (yymmdd)
c***revision date      871004   (yymmdd)
c
c    4 october 1987  bug fixed (?)  pws at lanl
c        changing line near end from 'maxset=max(maxset,nset)'.
c        replaced 'nset' with 'ninset' since 'maxset'=='mcu'
c        is surely the largest number of symmetry related atoms.
c
c***keywords           symmetry related atoms
c***author             saxe, paul (lanl)
c***source             @(#)atmset.f	5.1   11/6/94
c
c***purpose            enumerate the sets of symmetry related atoms.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       atmset
c
      implicit integer (a-z)
c
      integer symat(natoms)
      integer atprmt(natoms,nop)
c
      call izero(symat,natoms)
c
c     ----- search through atoms, enumerating related set if the
c           current atoms is not in a set
c
      nset=0
      maxset=0
      do 50 iatom=1,natoms
         if (symat(iatom).eq.0) then
            nset=nset+1
            do 10 op=1,nop
               symat(atprmt(iatom,op))=nset
 10         continue
            ninset=0
            do 20 jatom=1,natoms
               if (symat(jatom).eq.symat(iatom)) ninset=ninset+1
 20         continue
            maxset=max(maxset,ninset)
         end if
 50   continue
c
c
      return
      end
