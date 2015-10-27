*deck @(#)relate.f	5.1  11/6/94
      subroutine relate(symat,natoms,nsymat,ns,relatm,mcu)
c
c***begin prologue     relate
c***date written       870715   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           symmetry related atoms
c***author             saxe, paul (lanl)
c***source             @(#)relate.f	5.1   11/6/94
c
c***purpose            determine the number, and lists, of symmetry
c                      related atoms.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       relate
c
      implicit integer (a-z)
c
      integer symat(natoms)
      integer nsymat(ns)
      integer relatm(mcu,ns)
c
      do 10 set=1,ns
         n=0
         do 5 i=1,natoms
            if (symat(i).eq.set) then
               n=n+1
               relatm(n,set)=i
            end if
 5       continue
         nsymat(set)=n
 10   continue
c
c
      return
      end
