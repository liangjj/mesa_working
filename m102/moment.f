*deck @(#)moment.f	5.1  11/6/94
      subroutine moment(nprim,momatm,maxmom,natoms,nbtype,maxgrp)
c***begin prologue     moment.f
c***date written       870716   (yymmdd)
c***revision date      11/6/94
c
c***keywords           maxmimum angular momentum
c***author             saxe, paul (lanl)
c***source             @(#)moment.f	5.1   11/6/94
c
c***purpose            find the maximum angula momentum on each atom
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       moment.f
      implicit none
c     --- input variables -----
      integer natoms,nbtype
c     --- input arrays (unmodified) ---
      integer nprim(natoms,nbtype)
      integer maxmom(nbtype)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer momatm(natoms)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxgrp,atom,type
c
c     --- find the maxmimum angular momentum of basis functions
c         on each atom
      maxgrp=-1
      do 20 atom=1,natoms
         momatm(atom)=-1
         do 10 type=1,nbtype
            if (nprim(atom,type).gt.0) then
               momatm(atom)=max(momatm(atom),maxmom(type))
            end if
 10      continue
         maxgrp=max(maxgrp,momatm(atom))
 20   continue
c
c
      return
      end
