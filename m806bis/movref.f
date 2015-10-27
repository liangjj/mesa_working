*deck @(#)movref.f	1.2  7/30/91
      subroutine movref(bfcode,nbf,nrefs,orbcod,norbs,iout)
c
c***begin prologue  movref
c***date written   850102   (yymmdd)
c***revision date  yymmdd   (yymmdd)
c***keywords  drt, distinct row table
c
c***author  saxe, paul,    (lanl)
c***purpose  to shift the codes for a reference from bfcode in
c            scf ordering to a temporary array orbcod in drt
c            ordering
c***description
c
c
c***references
c
c***routines called  (none)
c***end prologue  movref
c
      implicit integer (a-z)
c
      integer bfcode(nbf,nrefs),orbcod(norbs,nrefs),iout(nbf)
c
      do 2 ref=1,nrefs
         do 1 bf=1,nbf
            if (iout(bf).gt.0) orbcod(iout(bf),ref)=bfcode(bf,ref)
    1    continue
    2 continue
c
c
      return
      end
