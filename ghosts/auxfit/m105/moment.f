*deck @(#)moment.f	1.1  11/20/92
      subroutine moment(nprim,momatm,maxmom,natoms,nbtype,maxgrp)
c
c***begin prologue     moment
c***date written       870716   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           maxmimum angular momentum
c***author             saxe, paul (lanl)
c***source             @(#)moment.f	1.1   11/20/92
c
c***purpose            find the maximum angula momentum on each atom
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       moment
c
      implicit integer (a-z)
c
      integer nprim(natoms,nbtype)
      integer momatm(natoms)
      integer maxmom(nbtype)
c
c     ----- find the maxmimum angular momentum of basis functions
c           on each atom
c
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
