*deck @(#)filchg.f	5.1  11/6/94
      subroutine filchg(natoms,ian,atmchg)
c***begin prologue     filchg
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           atomic charges, vector, float
c***author             martin, richard (lanl)
c***source
c***purpose            fills an array with default atomic charges.
c***description
c                      module to install default atomic charges into the array
c                      atmchg.  in this routine, these are merely the atomic
c                      numbers(ian).  later on, other routines might modify
c                      this data (e.g. a pseudo-potential).
c                      call filchg(natoms,ian,atmchg)
c                        natoms   the number of atoms.
c                        ian      integer atomic numbers, dimensioned natoms.
c                        atmchg   nuclear charges, dimensioned natoms.
c
c                      note that this routine could, and probably should,
c                      be replaced by a generic vector float routine.
c
c***references
c***routines called    float
c***end prologue       filchg
      implicit integer(a-z)
      real*8 atmchg(natoms)
      dimension ian(natoms)
c
c
      do 10 i=1,natoms
   10 atmchg(i)=float(ian(i))
c
c
      return
      end
