*deck @(#)exchng.f	5.1  11/6/94
      subroutine exchng(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,
     #     b,ints,sirow,sjrow,nrows,norbs,nlevs,orbfrm,nsym,nmax,
     $     nwks,nnp,irowsv,jrowsv,segsv,pagesv,iwtsv,jwtsv,traksv,
     #     acoef,bcoef,rbuf,ibuf,lnbuf,cutoff,n,ntotal,diag,
     $     imngrp,imxgrp,jmngrp,jmxgrp,ngroup,levp,unit)
c
c***begin prologue     exchng
c***date written       870807   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci coupling coefficients, exchng
c***author             saxe, paul (lanl)
c***source             @(#)exchng.f	5.1   11/6/94
c
c***purpose            to form one- and two-body loops and the
c     diagonalization tape.
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       exchng
c
c
      implicit integer (a-z)
c
      integer arc(4,nrows)
      integer wt(4,nrows)
      integer nlwks(nrows)
      integer ijadd(nnp)
      integer kadd(norbs,nsym)
      integer ladd(norbs,nsym)
      integer orbsym(norbs)
      integer b(nrows)
      integer irowsv(nlevs)
      integer segsv(nlevs)
      integer pagesv(nlevs)
      integer iwtsv(nlevs)
      integer jwtsv(nlevs)
      integer jrowsv(nlevs)
      integer traksv(nlevs)
      integer ibuf(2,lnbuf)
      integer imngrp(ngroup)
      integer imxgrp(ngroup)
      integer jmngrp(ngroup)
      integer jmxgrp(ngroup)
      real*8 ints(nmax)
      real*8 acoef(nlevs)
      real*8 bcoef(nlevs)
      real*8 rbuf(lnbuf)
      real*8 diag(nwks)
      real*8 cutoff
      character*16 unit
c
c     ----- local arrays -----
c
      real*8 coeffs(20,21)
      real*8 cfs(420)
      integer segmax(22)
      integer trak(228)
      integer jcond(228)
      integer kcond(228)
      integer nxtpag(228)
      integer iarc(228)
      integer jarc(228)
      integer segmin(22)
c
c     ----- local variables -----
c
      real*8 crite,root2,rootn2,toor2,toorn2
      real*8 d,dx,acof,bcof
      real*8 loop
      logical hqp
c
      common /io/ inp,iout
c
      equivalence (coeffs,cfs)
c
      data segmax/16,34,52,63,75,92,102,118,128,137,148,155,162,172,
     a            179,186,193,200,207,214,221,228/
      data (trak(i),i=1,150)
     #          / 1,  3,  1,  3,  1,  2,  9,  1,  1,  7
     a,           2,  9,  1,  7, 10,  9,  0,  4,  4,  3
     b,           0,  2,  4,  9,  0, 10,  0, 10,  4,  0
     c,           9,  3,  0, 10,  0,  4,  4,  9,  0, 10
     d,           0, 10,  4,  3,  0,  2,  4,  0,  3,  9
     e,           0, 10,  0,  0,  0,  0,  0,  0,  0,  0
     f,           0,  0,  0,  0,  0,  0,  0,  0,  1,  0
     g,           0,  0,  0,  0,  0,  0,  5,  5,  0,  0
     h,           0,  3,  5,  0,  3,  0,  0,  5,  0,  0
     i,           0,  0,  0,  0,  0,  1,  1,  0,  0,  0
     j,           1,  0,  0,  0,  0,  1,  1,  0,  0,  0
     k,           1,  0,  1,  0,  0,  1,  1,  0,  0,  0
     l,           1,  0,  1,  0,  0,  0,  1,  0,  0,  0
     m,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     n,           0,  0,  0,  0,  1,  0,  0,  0,  0,  0/
      data (trak(i),i=151,228)
     o           /0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     p,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     q,           0,  0,  0,  0,  0,  0,  0,  8,  0,  0
     r,           0,  0,  0,  0,  8,  0,  0,  0,  0,  0
     s,           0,  6,  0,  0,  0,  0,  0,  0,  6,  0
     t,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     u,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     v,           0,  0,  0,  0,  0,  0,  0,  0/
      data (jcond(i),i=1,150)
     #          /-1,  1, -1,  1,  1,  1,  1,  1, -1,  1
     a,           1,  1, -1,  1,  1,  1, -1,  1,  1,  1
     b,          -1,  1,  1,  1, -1,  1, -1,  1,  1,  1
     c,           1,  1, -1,  1, -1,  1,  1,  1, -1,  1
     d,          -1,  1,  1,  1, -1,  1,  1,  1,  1,  1
     e,          -1,  1,  0,  0,  0,  0,  0,  0,  0,  0
     f,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     g,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     h,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     i,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     j,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     k,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     l,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     m,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     n,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0/
      data (jcond(i),i=151,228)
     o           /0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     p,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     q,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     r,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     s,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     t,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     u,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     v,           0,  0,  0,  0,  0,  0,  0,  0/
      data (kcond(i),i=1,150)
     #          / 0,  1,  0,  1,  0,  1,  0,  0,  0,  1
     a,           1,  0,  0,  1,  1,  0,  0,  0,  0,  0
     b,           0,  1,  0,  0,  0,  1,  0,  1,  0,  1
     c,           0,  0,  0,  1,  0,  0,  0,  0,  0,  1
     d,           0,  1,  0,  0,  0,  1,  0,  1,  0,  0
     e,           0,  1,  0,  1,  1,  1,  0,  1,  1,  0
     f,           1,  1,  0,  0,  1,  1,  1,  0,  0,  1
     g,           1,  0,  1,  1,  0,  0,  1,  1,  1,  1
     h,           0,  0,  1,  1,  0,  1,  0,  1,  1,  1
     i,           1,  0,  0,  1,  0,  1,  1,  0,  0,  1
     j,           1,  0,  0,  1,  1,  1,  1,  0,  0,  1
     k,           1,  0,  1,  0,  1,  1,  1,  0,  0,  1
     l,           1,  0,  1,  0,  1,  0,  1,  0,  0,  1
     m,           0,  1,  0,  1,  1,  1,  0,  0,  1,  0
     n,           0,  1,  0,  0,  1,  1,  1,  0,  0,  0/
      data (kcond(i),i=151,228)
     o           /1,  0,  0,  1,  0,  0,  1,  0,  0,  0
     p,           1,  0,  0,  1,  0,  0,  1,  0,  0,  1
     q,           1,  0,  0,  0,  0,  0,  0,  0,  0,  0
     r,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     s,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     t,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     u,           0,  0,  0,  0,  0,  0,  0,  0,  0,  0
     v,           0,  0,  0,  0,  0,  0,  0,  0/
      data (nxtpag(i),i=1,150)
     #          / 3, 18,  2, 17, 10,  0,  5,  7,  2, 15
     a,           0,  5,  3, 16,  0,  4,  2, 11, 12,  7
     b,           2, 21, 12,  6,  3, 22,  2, 21, 11, 20
     c,           6,  7,  2, 21,  3, 13, 11,  6,  3, 22
     d,           2, 21, 11,  9,  3, 22, 13, 19,  9,  6
     e,           3, 22,  4, 22, 21,  0,  4, 21,  0,  4
     f,          22,  0,  4,  5, 22, 21,  0,  5,  7, 21
     g,           0,  5, 22,  0,  5,  6, 22, 21, 20,  0
     h,           6,  7, 21, 19,  9,  0,  6, 22, 19, 20
     i,           0,  6,  7, 21,  7, 20,  0,  8,  7, 21
     j,          20,  7,  8, 22, 21, 20,  0,  8,  7, 21
     k,          19,  9,  0,  8, 22, 19, 20,  8,  9, 22
     l,          19,  9,  0,  8, 22,  9, 19,  9, 10, 21
     m,          10, 22, 10,  0, 22, 21, 10, 11, 21, 11
     n,          12, 22, 13, 11,  0, 22, 21, 11, 12, 12/
      data (nxtpag(i),i=151,228)
     o          /21, 14, 12, 21, 12, 13, 22, 13, 14, 13
     p,          22, 13, 14, 21, 14, 12, 22, 13, 14, 22
     q,          21, 14, 15, 15,  0, 16, 15,  0, 15, 16
     r,           0, 16, 15, 16,  0, 16, 17, 17,  0, 18
     s,          17,  0, 17, 18,  0, 18, 17, 18,  0, 18
     t,          19,  0, 19, 20,  0, 19, 19, 20,  0, 20
     u,          19, 20,  0, 20, 21, 21,  0, 22, 21,  0
     v,          21, 22,  0, 22, 21, 22,  0, 22/
      data (jarc(i),i=1,150)
     #          / 2,  2,  3,  3,  4,  2,  2,  3,  4,  4
     a,           3,  3,  4,  4,  4,  4,  1,  2,  3,  1
     b,           2,  2,  4,  1,  2,  2,  3,  3,  4,  1
     c,           2,  3,  4,  4,  1,  2,  3,  1,  2,  2
     d,           3,  3,  4,  1,  3,  3,  4,  1,  2,  3
     e,           4,  4,  1,  2,  3,  2,  2,  4,  3,  3
     f,           4,  4,  4,  1,  2,  3,  2,  2,  3,  4
     g,           3,  3,  4,  4,  4,  1,  2,  3,  1,  2
     h,           2,  3,  4,  1,  2,  3,  3,  4,  2,  3
     i,           4,  4,  1,  2,  2,  1,  2,  2,  3,  4
     j,           2,  4,  1,  2,  3,  1,  2,  2,  3,  4
     k,           1,  2,  3,  3,  4,  2,  3,  4,  1,  3
     l,           1,  2,  3,  3,  4,  3,  3,  4,  1,  1
     m,           2,  1,  3,  1,  2,  3,  4,  1,  1,  2
     n,           3,  1,  2,  3,  1,  2,  3,  4,  1,  2/
      data (jarc(i),i=151,228)
     o           /1,  2,  3,  2,  4,  1,  1,  2,  3,  3
     p,           3,  4,  1,  1,  2,  3,  1,  2,  3,  2
     q,           3,  4,  1,  2,  1,  2,  3,  2,  4,  1
     r,           1,  2,  3,  3,  3,  4,  1,  2,  1,  2
     s,           3,  2,  4,  1,  1,  2,  3,  3,  3,  4
     t,           1,  3,  2,  3,  4,  3,  4,  1,  2,  2
     u,           2,  3,  4,  4,  1,  2,  1,  2,  3,  2
     v,           4,  1,  1,  2,  3,  3,  3,  4/
      data (iarc(i),i=1,150)
     #          / 1,  1,  1,  1,  1,  2,  2,  2,  2,  2
     a,           3,  3,  3,  3,  4,  4,  1,  1,  1,  2
     b,           2,  2,  2,  3,  3,  3,  3,  3,  3,  4
     c,           4,  4,  4,  4,  1,  1,  1,  2,  2,  2
     d,           2,  2,  2,  3,  3,  3,  3,  4,  4,  4
     e,           4,  4,  1,  1,  1,  2,  2,  2,  3,  3
     f,           3,  4,  4,  1,  1,  1,  2,  2,  2,  2
     g,           3,  3,  3,  4,  4,  1,  1,  1,  2,  2
     h,           2,  2,  2,  3,  3,  3,  3,  3,  4,  4
     i,           4,  4,  1,  1,  2,  3,  3,  3,  3,  3
     j,           4,  4,  1,  1,  1,  2,  2,  2,  2,  2
     k,           3,  3,  3,  3,  3,  4,  4,  4,  1,  1
     l,           2,  2,  2,  2,  2,  3,  4,  4,  1,  2
     m,           2,  3,  3,  4,  4,  4,  4,  1,  2,  2
     n,           2,  3,  3,  3,  4,  4,  4,  4,  1,  2/
      data (iarc(i),i=151,228)
     o           /3,  3,  3,  4,  4,  1,  2,  2,  2,  3
     p,           4,  4,  1,  2,  2,  2,  3,  3,  3,  4
     q,           4,  4,  1,  2,  3,  3,  3,  4,  4,  1
     r,           2,  2,  2,  3,  4,  4,  1,  2,  3,  3
     s,           3,  4,  4,  1,  2,  2,  2,  3,  4,  4
     t,           1,  1,  2,  2,  2,  3,  4,  1,  1,  2
     u,           3,  3,  3,  4,  1,  2,  3,  3,  3,  4
     v,           4,  1,  2,  2,  2,  3,  4,  4/
      save segmax,trak,jcond,kcond,nxtpag,jarc,iarc
c
c     ----- for second pass on hpq, flip i and j indices -----
c
      hqp=sirow.eq.1.and.sjrow.eq.2
c
c     ----- set up some constants -----
c
      crite=0.00001d+00
      root2=sqrt(2.0d+00)
      rootn2=-root2
      toor2=1.0d+00/root2
      toorn2=-toor2
c
c     ----- build the coefficient table -----
c
      call rzero(coeffs,20*21)
c
      do 140 i=3,20
         a = i-2
         coeffs(i,1) = sqrt(a/(a+1.0d+00))
         coeffs(i,2) = -coeffs(i,1)
         coeffs(i,3) = coeffs(i,1)/sqrt(2.0d+00)
         coeffs(i,4) = -coeffs(i,3)
         coeffs(i,5) = sqrt((a+1.0d+00)/a)
         coeffs(i,6) = -coeffs(i,5)
         coeffs(i,7) = coeffs(i,5)/sqrt(2.0d+00)
         coeffs(i,8) = -coeffs(i,7)
         coeffs(i,9) = sqrt((a+2.0d+00)/(a*2.0d+00))
         coeffs(i,10) = -coeffs(i,9)
         coeffs(i,11) = sqrt(a/(2.0d+00*(a+2.0d+00)))
         coeffs(i,12) = -coeffs(i,11)
         coeffs(i,13) = sqrt(2.0d+00/(a*(a+1.0d+00)))
         coeffs(i,14) = sqrt(a*(a+2.0d+00))/(a+1.0d+00)
         coeffs(i,15) = -sqrt(a*(a+2.0d+00))/(a+1.0d+00)
         coeffs(i,16) = sqrt((a-1.0d+00)*(a+2.0d+00)/(a*(a+1.0d+00)))
         coeffs(i,17) = -coeffs(i,16)
         coeffs(i,18) = -sqrt(2.0d+00/(a*(a+2.0d+00)))
         coeffs(i,19) = 1.0d+00/a
         coeffs(i,20) = -1.0d+00/a
         coeffs(i,21) = -sqrt(2.0d+00)/a
  140 continue
      do 141 i=2,22
         segmin(i)=segmax(i-1)
  141 continue
      segmin(1)=0
c
c     ----- initialize counters -----
c
      call rzero(diag,nwks)
c
c     ----- loop over opening levels of loops, starting at top of graph
c
      call iosys('rewind "guga integrals" on '//unit,0,0,0,' ')
c
      do 2000 group=1,ngroup
         call iosys('read real "guga integrals" from '//unit//
     $        ' without rewinding',nmax,ints,0,' ')
c
         imin=imngrp(group)
         imax=imxgrp(group)
      do 1000 level=imax+1,max(2,imin+1,levp+2),-1
c
c        ----- descend partial walks to the opening level -----
c
         if (level.eq.nlevs) then
c
c           ----- special situation for loops opening at dividing level
c
            hseg=4
            hlev=nlevs
            irow=sirow
            jrow=sjrow
            iwtsv(level)=0
            jwtsv(level)=0
            acoef(level)=1.0d+00
            page=1
            seg=segmin(page)
            segmx=segmax(page)
            go to 200
         end if
c
c        ----- initialize head section ------
c
         irow=sirow
         jrow=sjrow
         hlev=nlevs
         iwtsv(hlev)=0
         jwtsv(hlev)=0
         hseg=0
c
c        ----- start tree search of head sections -----
c
  101    continue
            hseg=hseg+1
            if (hseg.gt.4) then
               hlev=hlev+1
               if (hlev.gt.nlevs) go to 1000
               irow=irowsv(hlev-1)
               jrow=jrowsv(hlev-1)
               hseg=segsv(hlev-1)
               go to 101
            end if
c
c           ----- check that this segment exists on the graph -----
c
            jnxt=arc(hseg,jrow)
            if (jnxt.eq.0) go to 101
c
            inxt=arc(hseg,irow)
            if (inxt.eq.0) go to 101
c
c           ----- hop out of search if descended to desired level -----
c
            if (hlev-1.eq.level) go to 102
c
c           ----- otherwise, descend another level -----
c
            hlev=hlev-1
            irowsv(hlev)=irow
            jrowsv(hlev)=jrow
            segsv(hlev)=hseg
            iwtsv(hlev)=iwtsv(hlev+1)+wt(hseg,irow)
            jwtsv(hlev)=jwtsv(hlev+1)+wt(hseg,jrow)
            hseg=0
            irow=inxt
            jrow=jnxt
         go to 101
c
c     ----- initialize stacks for tree search -----
c
  102    continue
         irowsv(level)=irow
         jrowsv(level)=jrow
         iwtsv(level)=iwtsv(level+1)+wt(hseg,irow)
         jwtsv(level)=jwtsv(level+1)+wt(hseg,jrow)
         irow=inxt
         jrow=jnxt
         acoef(level)=1.0d+00
         page=1
         seg=segmin(page)
         segmx=segmax(page)
c
c        ----- start tree search in downward direction until we have
c              found all of the bottoms of loops eminating from the
c              given rows and page in the table.
c
 200    continue
c
         lev=level
         levm1=lev-1
         i=levm1
         ia=i*(i-1)/2
         is=orbsym(i)
c
         if (i.eq.imax) then
            jmax=jmxgrp(group)
         else
            jmax=i
         end if
         if (i.eq.imin) then
            jmin=jmngrp(group)
         else
            jmin=1
         end if
c
c        ----- initialize the stacks -----
c
         traksv(lev)=0
         acoef(lev)=1.0d+00
c
c        ----- test next segment in this page to see if has correct j arc
c
  201    continue
         seg=seg+1
         if (seg.gt.segmx) then
c
c           ----- have tried all segments at this level, so back up one
c
            lev=lev+1
            if (lev.gt.level) then
               irow=irowsv(level)
               jrow=jrowsv(level)
               go to 101
            end if
c
            levm1=lev-1
            seg=segsv(levm1)
            page=pagesv(levm1)
            segmx=segmax(page)
            irow=irowsv(levm1)
            jrow=jrowsv(levm1)
            go to 201
         end if
c
c        ----- check if arcs lead to valid rows -----
c
         inxt=arc(iarc(seg),irow)
         if (inxt.eq.0) go to 201
         jnxt=arc(jarc(seg),jrow)
         if (jnxt.eq.0) go to 201
c
c        ----- check for j or k level -----
c
         if (jcond(seg).gt.0) then
            if (levm1.gt.jmax) go to 201
            j=levm1
            ij=ijadd(ia+j)
            ijs=xor(is,orbsym(j))
         else if (jcond(seg).lt.0) then
            if (levm1.le.jmin) go to 201
c exchng
            if (levm1.le.levp) go to 201
c end
         end if
c
         if (kcond(seg).gt.0) then
c exchng
            if (levm1.gt.levp) go to 201
c end
            k=levm1
            ijk=ij+kadd(k,ijs+1)
            ijks=xor(ijs,orbsym(k))
         end if
c
c        ----- check for a change in tracks -----
c
         if (trak(seg).ne.0) then
            traksv(levm1)=trak(seg)
         else
            traksv(levm1)=traksv(lev)
         end if
      if (seg.gt.150) go to 2345
c             1   2   3   4   5   6   7   8   9  10
      go to(  1,  1,  1,  1,  1,  1, 44,  1,  3,  3,
     $        1, 45,  4,  4, 50, 51,  1, 40,  1,  1,
     $        6,  6,  7, 46,  9, 54,  2, 52, 41,  5,
     $       47,  8,  2, 53,  1,  1, 42, 48,  2, 52,
     $       36, 55, 43,  1, 11, 11, 12, 10, 13, 49,
     $        2, 53,  1, 77, 77, 77,  1, 79, 77,  1,
     $       80, 78,  1, 71, 67, 68, 87, 75, 83, 69,
     $       68, 76, 70, 82, 71, 71, 67, 68, 67, 87,
     $       75, 83, 69, 68, 83, 68, 76, 70, 69, 70,
     $       82, 71,  1,  6, 16,  6,  6, 17, 16, 74,
     $        8,  1,  1, 18, 19, 18, 18, 22, 24, 20,
     $       19, 24, 19, 23, 21, 20, 21,  1,  1, 11,
     $       11, 27, 11, 28, 81, 27, 13,  1,  1,  3,
     $        2,  4,  2,  1,  2,  2,  1, 71, 63, 72,
     $       84, 65, 85, 73, 29, 66, 64, 71,  1, 30), seg
      call lnkerr('seg 1')
 2345 continue
      go to( 56, 86,  1, 57,  1,  1, 56,  1, 32, 31,
     $       58,  1,  1, 59, 33, 34, 61, 35, 88, 62,
     $       60,  1,  1,  6,  1,  9,  2,  5,  2,  1,
     $        1,  2, 36, 11, 10,  2,  1,  6,  1,  9,
     $        2,  5,  2,  1,  1,  2, 36, 11, 10,  2,
     $        1,  1, 37, 38,  4,  2,  2,  1,  1,  2,
     $       39, 37,  3,  2,  1,  6,  1,  9,  2,  5,
     $        2,  1,  1,  2, 36, 11, 10,  2),  seg-150
      call lnkerr('seg 2')
c
  1      acoef(levm1)=acoef(lev)
         goto 120
  2      acoef(levm1)=-acoef(lev)
         goto 120
   3     junk = b(jrow) + 2
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   4     junk = b(jrow) + 83
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   5     junk = b(jrow) + 82
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   6     junk = b(jrow) + 261
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   7     junk = b(jrow) + 1
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   8     junk = b(jrow) + 102
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
   9     junk = b(jrow) + 362
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  10     junk = b(jrow) + 3
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  11     junk = b(jrow) + 263
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  12     junk = b(jrow) + 84
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  13     junk = b(jrow) + 23
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  16     junk = b(jrow) + 281
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  17     junk = b(jrow) + 402
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  18     junk = b(jrow) + 162
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  19     junk = b(jrow) + 222
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  20     junk = b(jrow) + 143
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  21     junk = b(jrow) + 42
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  22     junk = b(jrow) + 302
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  23     junk = b(jrow) + 303
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  24     junk = b(jrow) + 342
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  27     junk = b(jrow) + 283
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  28     junk=b(jrow)+404
         acoef(levm1)=acoef(lev)*cfs(junk)
         go to 120
  29     acoef(levm1) = acoef(lev) * root2
         go to 120
  30     junk = b(jrow) + 301
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  31     junk = b(jrow) + 304
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  32     junk = b(jrow) + 244
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  33     junk = b(jrow) + 322
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  34     junk = b(jrow) + 243
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  35     junk = b(jrow) + 242
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  36     junk = b(jrow) + 384
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  37     junk = b(jrow) + 262
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  38     junk = b(jrow) + 363
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  39     junk = b(jrow) + 383
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  86     junk = b(jrow) + 241
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  40     junk = b(jrow) + 122
         junk1=junk-61
         acoef(levm1) = acoef(lev) * cfs(junk)
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  41     junk1 = b(jrow) + 162
         acoef(levm1) = acoef(lev) * toorn2
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  42     junk = b(jrow) + 43
         junk1 = junk + 81
         acoef(levm1) = acoef(lev) * cfs(junk)
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  43     junk1 = b(jrow) + 222
         acoef(levm1) = acoef(lev) * toorn2
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  44     junk1=b(jrow)+221
         acoef(levm1) = acoef(lev) * toor2
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  45     junk1 = b(jrow) + 163
         acoef(levm1) = acoef(lev) * toor2
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  46     junk1 = b(jrow) + 162
         acoef(levm1) = acoef(lev) * toor2
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  47     junk = b(jrow) + 122
         junk1 = junk - 81
         acoef(levm1) = acoef(lev) * cfs(junk)
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  48     junk1 = b(jrow) + 222
         acoef(levm1) = acoef(lev) * toor2
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  49     junk = b(jrow) + 43
         junk1 = junk + 101
         acoef(levm1) = acoef(lev) * cfs(junk)
         bcoef(levm1) = acoef(lev) * cfs(junk1)
         go to 120
  50     acoef(levm1) = acoef(lev) + acoef(lev)
         d=0.5d+00
         go to 120
  51     acoef(levm1)=acoef(lev)*root2
         go to 120
  52     acoef(levm1) = -acoef(lev)
         d= -1.0d+00
         go to 120
  53     acoef(levm1) = -acoef(lev) - acoef(lev)
         d = -0.5d+00
         go to 120
  54     junk=b(jrow)+362
         d=1.0d+00/cfs(junk)
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  55     junk = b(jrow) + 384
         d=1.0d+00/cfs(junk)
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
  56     acoef(levm1) = acoef(lev)
         d = -1.0d+00
         go to 120
  57     junk = b(jrow) + 82
         acoef(levm1) = acoef(lev) * cfs(junk)
         d=-1.0d+00
         go to 120
  58     junk = b(jrow) + 3
         acoef(levm1) = acoef(lev) * cfs(junk)
         d=-1.0d+00
         go to 120
  59     junk = b(jrow) + 123
         acoef(levm1) = acoef(lev) * cfs(junk)
         d=-1.0d+00
         go to 120
  60     junk = b(jrow) + 222
         acoef(levm1) = acoef(lev) * cfs(junk)
         d=-1.0d+00
         go to 120
  61     junk = b(jrow) + 62
         acoef(levm1) = acoef(lev) * cfs(junk)
         d=-1.0d+00
         go to 120
  62     junk = b(jrow) + 162
         acoef(levm1) = acoef(lev) * cfs(junk)
         d=-1.0d+00
         go to 120
  63     junk = b(jrow) + 42
         junk1 = junk + 81
         acof = acoef(lev) * cfs(junk)
         bcof = bcoef(lev) * cfs(junk1)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm1) = d
         d = (acof - bcof) / d
         go to 120
  64     junk1 = b(jrow) + 222
         acof = acoef(lev) * toorn2
         bcof = bcoef(lev) * cfs(junk1)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm1) = d
         d = (acof - bcof) / d
         go to 120
  65     junk = b(jrow) + 123
         junk1 = junk - 61
         acof = acoef(lev) * cfs(junk)
         bcof = bcoef(lev) * cfs(junk1)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm1) = d
         d = (acof - bcof) / d
         go to 120
  66     junk1 = b(jrow) + 162
         acof = acoef(lev) * toorn2
         bcof = bcoef(lev) * cfs(junk1)
         d = acof + bcof
         if(abs(d).lt.crite) go to 110
         acoef(levm1) = d
         d = (acof - bcof) / d
         go to 120
  67     junk1 = b(jrow) + 162
         dx=acoef(lev)*toorn2
         d=dx+bcoef(lev)*cfs(junk1)
         if(abs(d).lt.crite) go to 111
         acoef(levm1) = d
         d=-(dx+dx)/d
         go to 120
  68     junk1 = b(jrow) + 222
         dx=acoef(lev)*toorn2
         d=dx+bcoef(lev)*cfs(junk1)
         if(abs(d).lt.crite) go to 111
         acoef(levm1) = d
         d=-(dx+dx)/d
         go to 120
  69     junk = b(jrow) + 62
         junk1 = junk + 81
         dx=acoef(lev)*cfs(junk)
         d=dx+bcoef(lev)*cfs(junk1)
         if(abs(d).lt.crite) go to 111
         acoef(levm1) = d
         d=-(dx+dx)/d
         go to 120
  70     junk = b(jrow) + 143
         junk1 = junk - 101
         dx=acoef(lev)*cfs(junk)
         d=dx+bcoef(lev)*cfs(junk1)
         if(abs(d).lt.crite) go to 111
         acoef(levm1) = d
         d=-(dx+dx)/d
         go to 120
  87     junk1 = b(jrow) + 162
         dx=acoef(lev)*toorn2
         d=dx+bcoef(lev)*cfs(junk1)
         if(abs(d).lt.crite) go to 111
         acoef(levm1)=d
         d=-(dx+dx)/d
         go to 120
  71     acoef(levm1) = acoef(lev)
         bcoef(levm1) = bcoef(lev)
         go to 120
  72     junk1 = b(jrow) + 322
         acoef(levm1) = -acoef(lev)
         bcoef(levm1) = bcoef(lev) * cfs(junk1)
         go to 120
  73     junk1 = b(jrow) + 323
         acoef(levm1) = -acoef(lev)
         bcoef(levm1) = bcoef(lev) * cfs(junk1)
         go to 120
  74     junk=b(jrow)+21
         acoef(levm1)=acoef(lev)*cfs(junk)
         go to 120
  75     junk1 = b(jrow) + 302
         acoef(levm1) = acoef(lev)
         bcoef(levm1) = bcoef(lev) * cfs(junk1)
         go to 120
  76     junk1 = b(jrow) + 303
         acoef(levm1) = acoef(lev)
         bcoef(levm1) = bcoef(lev) * cfs(junk1)
         go to 120
  77     acoef(levm1)=acoef(lev)*toorn2
         d=-2.0d+00
         go to 120
  78     acoef(levm1)=acoef(lev)*rootn2
         d=-2.0d+00
         go to 120
  79     junk=b(jrow)+62
         acoef(levm1)=acoef(lev)*cfs(junk)
         d=-2.0d+00
         go to 120
  80     junk=b(jrow)+143
         acoef(levm1)=acoef(lev)*cfs(junk)
         d=-2.0d+00
         go to 120
  81     junk=b(jrow)+104
         acoef(levm1)=acoef(lev)*cfs(junk)
         go to 120
  82     acoef(levm1) = acoef(lev) * rootn2
         d=-2.0d+00
         go to 120
  83     junk = b(jrow) + 342
         acoef(levm1) = bcoef(lev) * cfs(junk)
         go to 120
  84     junk = b(jrow) + 243
         acoef(levm1) = bcoef(lev) * cfs(junk)
         go to 120
  85     junk = b(jrow) + 242
         acoef(levm1) = bcoef(lev) * cfs(junk)
         go to 120
  88     junk = b(jrow) + 323
         acoef(levm1) = acoef(lev) * cfs(junk)
         go to 120
 110     traksv(levm1)=3
         acoef(levm1)=acof-bcof
         go to 120
111      traksv(levm1) = 2
         acoef(levm1)=-(dx+dx)
  120    continue
c
c        ----- update stacks -----
c
         newpag=nxtpag(seg)
         if (newpag.eq.0) then
            l=levm1
c
            if (orbsym(l).ne.is) go to 201
c
c           ----- have finished a loop, so check if the rows are
c                 the same. (because of 't' values from interacting
c                 spaces, the rows may not be the same, though their
c                 a,b and s values must be.
c
            if (inxt.ne.jnxt) then
c
c              ----- push row, etc on stack and search down for closing -----
c
               tlev=lev-1
               irowsv(tlev)=irow
               jrowsv(tlev)=jrow
               iwtsv(tlev)=iwtsv(tlev+1)+wt(iarc(seg),irow)
               jwtsv(tlev)=jwtsv(tlev+1)+wt(jarc(seg),jrow)
               tarc=0
               itrow=inxt
               jtrow=jnxt
c
  400          continue
               tarc=tarc+1
               if (tarc.gt.4) then
                  tlev=tlev+1
                  if (tlev.ge.lev) go to 201
                  itrow=irowsv(tlev-1)
                  jtrow=jrowsv(tlev-1)
                  tarc=segsv(tlev-1)
                  go to 400
               end if
c
               inxt=arc(tarc,itrow)
               if (inxt.eq.0) go to 400
c
               jnxt=arc(tarc,jtrow)
               if (jnxt.eq.0) go to 400
c
               if (inxt.eq.jnxt) then
                  iwt=iwtsv(tlev)+wt(tarc,itrow)
                  jwt=jwtsv(tlev)+wt(tarc,jtrow)
                  ijkl=ijk+ladd(l,ijks+1)
c
                  go to (411,412,413,414,415,416,417,418,419,420),
     #                 traksv(levm1)
                  call lnkerr('bad track')
c
 411              continue
                     loop=ints(ijkl+1)*acoef(levm1)
                     go to 430
  412             continue
c exchng                     loop=ints(ijkl+2)*acoef(levm1)
                     go to 430
  413             continue
c exchng
                     if (i.ne.j.and.k.ne.l) then
                        loop=ints(ijkl+3)*acoef(levm1)
                     end if
c end
                     go to 430
  414             continue
                     loop=acoef(levm1)*(ints(ijkl+1)+d*ints(ijkl+3))
                     go to 430
  415             continue
c exchng                     loop=acoef(levm1)*(ints(ijkl+3)+d*ints(ijkl+2))
                     if (i.ne.k) then
                        loop=acoef(levm1)*ints(ijkl+3)
                     end if
c end
                     go to 430
  416             continue
c exchng                     loop=acoef(levm1)*(ints(ijkl+3)+ints(ijkl+2))
                     if (i.ne.k) then
                        loop=acoef(levm1)*ints(ijkl+3)
                     end if
c end
                     go to 430
  417             continue
c exchng
                     if (i.ne.k) then
                        loop=acoef(levm1)*(ints(ijkl+3)+ints(ijkl+1))
                     end if
c end
                     go to 430
  418             continue
c exchng                     loop=acoef(levm1)*(ints(ijkl+1)+
c     $                 ints(ijkl+2)+ints(ijkl+3))
                     if (i.ne.k) then
                        loop=acoef(levm1)*(ints(ijkl+1)+ints(ijkl+3))
                     end if
c end
                     go to 430
  419             continue
c exchng                     loop=acoef(levm1)*
c     $                 (ints(ijkl+1)+d*ints(ijkl+2))
                     loop=acoef(levm1)*ints(ijkl+1)
                     go to 430
  420             continue
c exchng                     loop=acoef(levm1)*(ints(ijkl+2)+d*ints(ijkl+1))
                     loop=acoef(levm1)*d*ints(ijkl+1)
  430             continue
                  if (abs(loop).gt.cutoff) then
                     if (iwt.eq.jwt) then
                        do 432 walk=1,nlwks(inxt)
                           diag(iwt+walk)=diag(iwt+walk)+loop
 432                    continue
                     else
                        do 435 walk=1,nlwks(inxt)
                           n=n+1
                           if (n.gt.lnbuf) then
                              ntotal=ntotal+lnbuf
                              call iosys('write integer buffers to '//
     $                             'hamiltonian without rewinding',
     $                             2*lnbuf,ibuf,0,' ')
                              call iosys('write integer buffers to '//
     $                             'hamiltonian without rewinding',
     $                             wptoin(lnbuf),rbuf,0,' ')
                              n=1
                           end if
                           if (hqp) then
                              ibuf(2,n)=jwt+walk
                              ibuf(1,n)=iwt+walk
                           else
                              ibuf(1,n)=jwt+walk
                              ibuf(2,n)=iwt+walk
                           end if
                           rbuf(n)=loop
 435                    continue
                     end if
                  end if
                  go to 400
c
               else
                  tlev=tlev-1
                  if (tlev.le.1) call lnkerr('error in tail segments')
                  irowsv(tlev)=itrow
                  jrowsv(tlev)=jtrow
                  iwtsv(tlev)=iwtsv(tlev+1)+wt(tarc,itrow)
                  jwtsv(tlev)=jwtsv(tlev+1)+wt(tarc,jtrow)
                  segsv(tlev)=tarc
                  tarc=0
                  itrow=inxt
                  jtrow=jnxt
                  go to 400
               end if
c
            else
c------------------------------------------------------------------
c
               iwt=iwtsv(lev)+wt(iarc(seg),irow)
               jwt=jwtsv(lev)+wt(jarc(seg),jrow)
               ijkl=ijk+ladd(l,ijks+1)
c
               go to (211,212,213,214,215,216,217,218,219,220),
     #                                                    traksv(levm1)
               call lnkerr('bad track')
c
c
 211           continue
                  loop=ints(ijkl+1)*acoef(levm1)
                  go to 230
 212           continue
c exchng                     loop=ints(ijkl+2)*acoef(levm1)
                  go to 230
 213           continue
c exchng
                  if (i.ne.j.and.k.ne.l) then
                     loop=ints(ijkl+3)*acoef(levm1)
                  end if
c end
                  go to 230
  214          continue
                  loop=acoef(levm1)*(ints(ijkl+1)+d*ints(ijkl+3))
                  go to 230
  215          continue
c exchng                  loop=acoef(levm1)*(ints(ijkl+3)+d*ints(ijkl+2))
                  if (i.ne.k) then
                     loop=acoef(levm1)*ints(ijkl+3)
                  end if
c end
                  go to 230
  216          continue
c exchng                  loop=acoef(levm1)*(ints(ijkl+3)+ints(ijkl+2))
                  if (i.ne.k) then
                     loop=acoef(levm1)*ints(ijkl+3)
                  end if
c end
                  go to 230
  217          continue
c exchng
                  if (i.ne.k) then
                     loop=acoef(levm1)*(ints(ijkl+3)+ints(ijkl+1))
                  end if
c end
                  go to 230
  218          continue
c exchng                  loop=acoef(levm1)*(ints(ijkl+1)+
c     $              ints(ijkl+2)+ints(ijkl+3))
                  if (i.ne.k) then
                     loop=acoef(levm1)*(ints(ijkl+1)+ints(ijkl+3))
                  end if
c end
                  go to 230
  219          continue
c exchng                  loop=acoef(levm1)*
c     $              (ints(ijkl+1)+d*ints(ijkl+2))
                  loop=acoef(levm1)*ints(ijkl+1)
                  go to 230
  220          continue
c exchng                  loop=acoef(levm1)*(ints(ijkl+2)+d*ints(ijkl+1))
                  loop=acoef(levm1)*d*ints(ijkl+1)
  230          continue
               if (abs(loop).gt.cutoff) then
                  if (iwt.eq.jwt) then
                     do 232 walk=1,nlwks(inxt)
                        diag(iwt+walk)=diag(iwt+walk)+loop
 232                 continue
                  else
                     do 235 walk=1,nlwks(inxt)
                        n=n+1
                        if (n.gt.lnbuf) then
                           ntotal=ntotal+lnbuf
                           call iosys('write integer buffers to '//
     $                          'hamiltonian without rewinding',
     $                          2*lnbuf,ibuf,0,' ')
                           call iosys('write integer buffers to '//
     $                          'hamiltonian without rewinding',
     $                          wptoin(lnbuf),rbuf,0,' ')
                           n=1
                        end if
                        if (hqp) then
                           ibuf(2,n)=jwt+walk
                           ibuf(1,n)=iwt+walk
                        else
                           ibuf(1,n)=jwt+walk
                           ibuf(2,n)=iwt+walk
                        end if
                        rbuf(n)=loop
 235                 continue
                  end if
               end if
c
c           ----- go back up and try next segment -----
c
               go to 201
            end if
         end if
c
c        ----- update stacks and then go down one level -----
c
         if (lev.le.2) go to 201
c
         irowsv(levm1)=irow
         jrowsv(levm1)=jrow
         pagesv(levm1)=page
         iwtsv(levm1)=iwtsv(lev)+wt(iarc(seg),irow)
         jwtsv(levm1)=jwtsv(lev)+wt(jarc(seg),jrow)
         segsv(levm1)=seg
         irow=inxt
         jrow=jnxt
         page=newpag
c
         lev=lev-1
         levm1=lev-1
         seg=segmin(page)
         segmx=segmax(page)
         go to 201
 1000 continue
 2000 continue
c
c
      return
      end
