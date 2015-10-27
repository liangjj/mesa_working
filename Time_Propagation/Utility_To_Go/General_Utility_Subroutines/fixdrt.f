*deck @(#)fixdrt.f	5.1  11/6/94
      subroutine fixdrt(orbsym,norbs,ijgrp,nnp,imngrp,imxgrp,
     $     jmngrp,jmxgrp,ngroup)
c
c***begin prologue     fixdrt
c***date written       870806   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           drt
c***author             saxe, paul (lanl)
c***source             @(#)fixdrt.f	5.1   11/6/94
c
c***purpose            to fixup and augment the drt information.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       fixdrt
c
      implicit integer (a-z)
c
      integer orbsym(norbs)
      integer ijgrp(nnp)
      integer imngrp(ngroup)
      integer imxgrp(ngroup)
      integer jmngrp(ngroup)
      integer jmxgrp(ngroup)
c
c     ----- start the symmetry numbering from 0 -----
c
      do 1 i=1,norbs
         orbsym(i)=orbsym(i)-1
    1 continue
c
c     ----- fill in the minimum and maximum i and j value per group
c
      group=1
      imxgrp(group)=norbs
      jmxgrp(group)=norbs
      do 25 i=norbs,1,-1
         do 20 j=i,1,-1
            ij=i*(i-1)/2+j
            if(ijgrp(ij).eq.group) go to 20
            if (i.eq.j.and.i.ne.imxgrp(group)) go to 22
            imngrp(group)=i
            jmngrp(group)=j+1
            go to 24
  22        continue
            imngrp(group)=i+1
            jmngrp(group)=1
  24        continue
            group=group+1
            imxgrp(group)=i
            jmxgrp(group)=j
  20     continue
  25  continue
      imngrp(ngroup)=1
      jmngrp(ngroup)=1
c
c
      return
      end
