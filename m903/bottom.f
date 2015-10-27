*deck @(#)bottom.f	5.1  11/6/94
      subroutine bottom(idnarc,idnwt,ia,ib,idnnwk,inrows,
     #                  jdnarc,jdnwt,ja,jb,jdnnwk,jnrows,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #                  irowst,jrowst,levdiv,nlevs,cfs,
     #                  res,ires,wpti4,binsiz,number)
c
c***begin prologue    bottom
c***date written      860910   (yymmdd)
c***revision date     yymmdd   (yymmdd)
c***keywords          ci, coupling coefficients, guga
c***author            saxe, paul,    (lanl)
c***purpose           to evaluate the bottom portions of loops.
c***description
c
c          i and j refer to distinct partial walks on distinct drt's. the
c      drt arrays (dnarc, dnwt, a, b, and dnnwk as well as nrows) are
c      prefixed with  the i or j respectively. irowst and jrowst are the
c      two rows to start with. this routine will then evaluate all bottoms
c      of one-electron loops eminating from these two rows, using a lexical
c      ordering scheme. on return, ires(1,*) and (2,*) contain the i and j
c      partial walk lexical weights, res(4,*) the partial loop coefficient,
c      ires(4,*) the j-index of the orbital (or one-electron integral) involved
c      , and number the number of loop bottoms found.
c
c***references
c
c***routines called
c***end prologue      bottom
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c           these are equivalenced in the calling routine.
      real*8 res(4,binsiz)
      integer ires(wpti4,binsiz)
c
c     ----- external arrays not modified -----
c
      real*8 cfs(420)
      integer idnarc(4,inrows),idnwt(4,inrows),ia(inrows)
      integer ib(inrows),idnnwk(inrows)
      integer jdnarc(4,jnrows),jdnwt(4,jnrows),ja(jnrows)
      integer jb(jnrows),jdnnwk(jnrows)
c
c     ----- external arrays used for scratch -----
c
      real*8 acoef(nlevs)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
c
c
c     ----- local arrays -----
c
      integer segmax(5),segmin(5),nxtpag(39),iarc(39),jarc(39)
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
      equivalence (segmin(2),segmax(1))
c
      data segmin(1) /0/
      data segmax / 11,18,25,32,39 /
      data nxtpag
     #  /  3, 2, 0, 2, 0, 3, 0, 4, 5, 5, 4,
     #     2, 2, 0, 3, 2, 0, 2,
     #     3, 0, 3, 2, 3, 0, 3,
     #     4, 0, 4, 5, 4, 0, 4,
     #     5, 0, 5, 4, 0, 5, 5/
      data   jarc
     #  /  2, 3, 2, 4, 3, 4, 4, 1, 1, 2, 3,
     #     1, 2, 1, 2, 3, 2, 4,
     #     1, 1, 2, 3, 3, 3, 4,
     #     1, 2, 2, 2, 3, 4, 4,
     #     1, 3, 2, 3, 4, 3, 4/
      data   iarc
     #  /  1, 1, 2, 2, 3, 3, 4, 2, 3, 4, 4,
     #     1, 2, 3, 3, 3, 4, 4,
     #     1, 2, 2, 2, 3, 4, 4,
     #     1, 1, 2, 3, 3, 3, 4,
     #     1, 1, 2, 2, 2, 3, 4/
      save segmin,segmax,nxtpag,jarc,iarc
c
      number=0
c
c     ----- initialize stacks for tree search -----
c
      level=levdiv
      irow=irowst
      jrow=jrowst
      iwtsv(level)=0
      jwtsv(level)=0
      acoef(level)=1.0d+00
c
c     ----- determine which page in the table is the start -----
c
      deltaa=ja(jrow)-ia(irow)
      deltab=jb(jrow)-ib(irow)
c
      if (deltaa.eq.0.and.deltab.eq.-1) then
         page=3
      else if (deltaa.eq.-1.and.deltab.eq.1) then
         page=2
      else if (deltaa.eq.0.and.deltab.eq.1) then
         page=4
      else if (deltaa.eq.1.and.deltab.eq.-1) then
         page=5
      else
         call lnkerr('illegal i and j rows in bottom')
      end if
c
      seg=segmin(page)
      segmx=segmax(page)
c
c     ----- start tree search in downward direction until we have
c           found all of the bottoms of loops eminating from the
c           given rows and page in the table.
c
  201 continue
c
         seg=seg+1
         if (seg.gt.segmx) then
c
c           ----- exhausted this page of the table, so back up a level
c
            level=level+1
            if (level.gt.levdiv) then
               go to 9000
            end if
c
            seg=segsv(level-1)
            page=pagesv(level-1)
            segmx=segmax(page)
            irow=irowsv(level-1)
            jrow=jrowsv(level-1)
            go to 201
         end if
c
c        ----- check that this segment is valid on the graph -----
c
         jnxt=jdnarc(jarc(seg),jrow)
         if (jnxt.eq.0) go to 201
c
         inxt=idnarc(iarc(seg),irow)
         if (inxt.eq.0) go to 201
c
c        ----- evaluate this segment -----
c
         go to
     #   ( 1, 1, 1, 3, 1, 4,50, 1, 1, 3, 4,
     #     1, 6, 1, 9, 2, 5, 2,
     #     1, 1, 2,36,11,10, 2,
     #     1, 1, 2,15,16, 3, 2,
     #     1, 1,16,17, 4, 2, 2), seg
c
    1    continue
            acoef(level-1)=acoef(level)
            goto 120
    2    continue
            acoef(level-1)=-acoef(level)
            goto 120
    3    continue
            junk = jb(jrow) + 2
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
    4    continue
            junk = jb(jrow) + 83
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
    5    continue
            junk = jb(jrow) + 82
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
    6    continue
            junk = jb(jrow) + 261
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
    9    continue
            junk = jb(jrow) + 362
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
   10    continue
            junk = jb(jrow) + 3
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
   11    continue
            junk = jb(jrow) + 263
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
   15    continue
            acoef(level-1)=acoef(level)*cfs(jb(jrow)+383)
            go to 120
   16    continue
            acoef(level-1)=acoef(level)*cfs(jb(jrow)+262)
            go to 120
   17    continue
            acoef(level-1)=acoef(level)*cfs(jb(jrow)+363)
            go to 120
   36    continue
            junk = jb(jrow) + 384
            acoef(level-1) = acoef(level) * cfs(junk)
            go to 120
   50    continue
            acoef(level-1) = acoef(level) + acoef(level)
            go to 120
  120    continue
c
c        ----- update stacks -----
c
         newpag=nxtpag(seg)
c
         if (newpag.eq.0) then
c
c           ----- have finished a loop, check for being at bottom
c                 of graph
c
            if (level.eq.2) then
               number=number+1
               if (number.gt.binsiz) then
                  call lnkerr('bin size to small in bottom')
               end if
               ires(1,number)=iwtsv(level)+idnwt(iarc(seg),irow)
               ires(2,number)=jwtsv(level)+jdnwt(jarc(seg),jrow)
               ires(3,number)=level-1
               res(4,number)=acoef(level-1)
               go to 201
            end if
c
c           ----- search for tails to bottom of graph -----
c
            irowsv(level-1)=irow
            jrowsv(level-1)=jrow
            iwtsv(level-1)=iwtsv(level)+idnwt(iarc(seg),irow)
            jwtsv(level-1)=jwtsv(level)+jdnwt(jarc(seg),jrow)
            tlev=level-1
            irow=inxt
            jrow=jnxt
            tseg=0
c
c           ----- search down tail segments to bottom of graph -----
c
  401       continue
               tseg=tseg+1
               if (tseg.gt.4) then
c
c                 ----- exhausted arcs at this level, back up one -----
c
                  tlev=tlev+1
                  if (tlev.ge.level) then
                     irow=irowsv(level-1)
                     jrow=jrowsv(level-1)
                     go to 201
                  end if
c
                  tseg=segsv(tlev-1)
                  irow=irowsv(tlev-1)
                  jrow=jrowsv(tlev-1)
                  go to 401
               end if
c
c              ----- check that this segment leads down on both i and j
c
               jnxt=jdnarc(tseg,jrow)
               if (jnxt.eq.0) go to 401
c
               inxt=idnarc(tseg,irow)
               if (inxt.eq.0) go to 401
c
c              ----- check for reaching bottom of graph -----
c
               if (tlev.eq.2) then
                  number=number+1
                  if (number.gt.binsiz) then
                     call lnkerr('bin size to small in bottom')
                  end if
                  ires(1,number)=iwtsv(tlev)+idnwt(tseg,irow)
                  ires(2,number)=jwtsv(tlev)+jdnwt(tseg,jrow)
                  ires(3,number)=level-1
                  res(4,number)=acoef(level-1)
                  go to 401
               end if
c
c              ----- update stacks and descend -----
c
               irowsv(tlev-1)=irow
               jrowsv(tlev-1)=jrow
               iwtsv(tlev-1)=iwtsv(tlev)+idnwt(tseg,irow)
               jwtsv(tlev-1)=jwtsv(tlev)+jdnwt(tseg,jrow)
               segsv(tlev-1)=tseg
               irow=inxt
               jrow=jnxt
               tlev=tlev-1
               tseg=0
            go to 401
         end if
c
c        ----- update stacks and then go down one level -----
c
         if (level.le.2) go to 201
c
         irowsv(level-1)=irow
         jrowsv(level-1)=jrow
         pagesv(level-1)=page
         iwtsv(level-1)=iwtsv(level)+idnwt(iarc(seg),irow)
         jwtsv(level-1)=jwtsv(level)+jdnwt(jarc(seg),jrow)
         segsv(level-1)=seg
         irow=inxt
         jrow=jnxt
         page=newpag
c
         level=level-1
         seg=segmin(page)
         segmx=segmax(page)
      go to 201
c
c     ----- the return is down here for convenience -----
c
 9000 continue
c
      return
      end
