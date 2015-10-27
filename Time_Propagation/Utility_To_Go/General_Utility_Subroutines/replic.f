*deck @(#)replic.f	5.1  11/6/94
      subroutine replic(old,new,ptprim,ptcont,nprim,ncont,atomz,
     $                  natoms,ntypes)
c***begin prologue     replic
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           basis set
c***author             saxe, paul (lanl)
c***source
c***purpose            replicates the bookkeeping arrays associated with a
c                      basis set.
c***description
c                      call replic(old,new,ptprim,ptcont,nprim,ncont,atomz,
c                                  natoms,ntypes)
c
c***references
c***routines called    (none)
c***end prologue       replic
      implicit integer (a-z)
c
      real*8 atomz(natoms)
      integer ptprim(natoms,ntypes),ptcont(natoms,ntypes)
      integer nprim(natoms,ntypes),ncont(natoms,ntypes)
c
      atomz(new)=atomz(old)
      do 1 type=1,ntypes
         ptprim(new,type)=ptprim(old,type)
         ptcont(new,type)=ptcont(old,type)
         nprim(new,type)=nprim(old,type)
         ncont(new,type)=ncont(old,type)
    1 continue
c
      return
      end
