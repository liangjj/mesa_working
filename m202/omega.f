*deck @(#)omega.f	5.1  11/6/94
      subroutine omega(natoms,rotmat,trvec,c,cscr,ian,toang,
     $                 usesym,prtsym,pgrp)
c***begin prologue     omega.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           symmetry
c***author             martin, richard (lanl)
c***source             @(#)omega.f	5.1   11/6/94
c***purpose            the final routine in the symmetry package.
c                      communicates symmetry information to the rwf.
c***description
c***references
c***routines called    iosys(io), noones(m202), numdof(m202), numer(symm),
c                      deornt(m202), subgrp(m202),
c                      corpr1(util)
c***end prologue       omega.f
      implicit none
c     --- input variables -----
      integer natoms
      character*4 pgrp
      logical prtsym,usesym
      real*8 toang
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 rotmat(3,3),trvec(3),c(3,natoms),cscr(3,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,at
      real*8 toler,tol2
c
      common/io/inp,iout
      common/tol/toler,tol2
c
c
 1010 format(1x,'rotation matrix:')
 1020 format(1x,22x,3f12.6)
c
c     --- omega is the final routine in the symmetry package.
c     it places the character strings molfor and fwg on r/w file
c     so that they may be read by the archiver.
c     it monitors the framework group during optimizations and turns
c     symmetry off if there is any change from point to point.
c     using the translation vector and rotation matrix provided by
c     ptgrp, it reorints the coordinates in blank common (array c
c     in the calling arguments).
c     it edits the list of operations so that only cartesian two-fold
c     operations remain.  it orders the operations written to the r/w
c     file such that the first nop1 define the largest concise abelian
c     subgroup while the first nop2 (nop2 >= nop1) define the longest
c     abelian subgroup.
c
c
c     --- if symmetry is in use, translate and rotate the molecule to the
c         standard orientation.
      if(usesym)then
c        --- translate.
         do 10 at=1,natoms
            call vadd(c(1,at),c(1,at),trvec,3)
   10    continue
c        --- skip the rotation if the point group is ci.
         if(pgrp.ne.'ci') then
            call ebc(cscr,rotmat,c,3,3,natoms)
            call vclean(cscr,toler,3*natoms)
            call vmove(c,cscr,3*natoms)
         endif
c
c        --- perhaps print the rotation matrix.
         if(prtsym) then
            write(iout,1010)
            do 20 i=1,3
               write(iout,1020) (rotmat(i,j),j=1,3)
   20       continue
         endif
      else
         call deornt(rotmat,0)
      endif
c
c
      return
      end
