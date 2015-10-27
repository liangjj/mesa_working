*deck @(#)top.f	5.1  11/6/94
      subroutine top(iuparc,iupwt,ia,ib,iupnwk,inrows,
     #               juparc,jupwt,ja,jb,jupnwk,jnrows,
     #               acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #               irowst,jrowst,levdiv,nlevs,cfs,
     #               res,ires,wpti4,binsiz,number)
c
c***begin prologue    top
c***date written      860910   (yymmdd)
c***revision date     yymmdd   (yymmdd)
c***keywords          ci, coupling coefficients, guga
c***author            saxe, paul,    (lanl)
c***purpose           to evaluate the tops of loops in reverse-lexical order
c***description
c
c          top evaluates the partial one-body coupling coefficient for
c       the top portion of loops eminating from irowst and jrowst at
c       levdiv, using two different drt's (uparc, upwt, a, b, upnwk,
c       and nrows, prefixed by i or j). the reverse-lexical partial
c       walk numbers are returned in ires(1,*) and (2,*), the partial
c       coefficient in res(4,*), and the i'th orbital index in ires(3,*).
c       number is the number of tops of loops found.
c
c***references
c
c***routines called
c***end prologue      top
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c           these are equivalenced in the calling routine.
      real*8 res(4,binsiz)
      integer ires(wpti4,binsiz)
c
c
c     ----- external arrays not modified -----
c
      real*8 cfs(420)
      integer iuparc(4,inrows),iupwt(4,inrows),ia(inrows)
      integer ib(inrows),iupnwk(inrows)
      integer juparc(4,jnrows),jupwt(4,jnrows),ja(jnrows)
      integer jb(jnrows),jupnwk(jnrows)
c
c     ----- external arrays used for scratch -----
c
      real*8 acoef(nlevs)
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs),jwtsv(nlevs)
      integer pagesv(nlevs),segsv(nlevs)
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
     #  /  0, 0, 0, 2, 2, 3, 3, 4, 4, 5, 5,
     #     0, 0, 2, 2, 2, 2, 3,
     #     0, 0, 3, 3, 3, 3, 2,
     #     0, 0, 4, 4, 4, 4, 5,
     #     0, 0, 5, 5, 5, 5, 4/
      data   jarc
     #  /  2, 3, 4, 1, 2, 1, 3, 2, 4, 3, 4,
     #     3, 4, 1, 2, 3, 4, 3,
     #     2, 4, 1, 2, 3, 4, 2,
     #     1, 3, 1, 2, 3, 4, 3,
     #     1, 2, 1, 2, 3, 4, 2/
      data   iarc
     #  /  2, 3, 4, 3, 4, 2, 4, 1, 3, 1, 2,
     #     1, 2, 1, 2, 3, 4, 2,
     #     1, 3, 1, 2, 3, 4, 3,
     #     2, 4, 1, 2, 3, 4, 2,
     #     3, 4, 1, 2, 3, 4, 3/
c
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
      if (deltaa.eq.-1.and.deltab.eq.1) then
         page=2
      else if (deltaa.eq.0.and.deltab.eq.-1) then
         page=3
      else if (deltaa.eq.0.and.deltab.eq.1) then
         page=4
      else if (deltaa.eq.1.and.deltab.eq.-1) then
         page=5
      else
         call lnkerr('illegal i and j rows in top')
      end if
c
      seg=segmin(page)
      segmx=segmax(page)
c
c     ----- start tree search in upward direction until we have
c           found all of the tops of loops eminating from the
c           given rows and page in the table.
c
  201 continue
c
         seg=seg+1
         if (seg.gt.segmx) then
c
c           ----- exhausted this page of the table, so back up a level
c
            level=level-1
            if (level.lt.levdiv) then
               go to 9000
            end if
c
            seg=segsv(level+1)
            page=pagesv(level+1)
            segmx=segmax(page)
            irow=irowsv(level+1)
            jrow=jrowsv(level+1)
            go to 201
         end if
c
c        ----- check that this segment is valid on the graph -----
c
         jnxt=juparc(jarc(seg),jrow)
         if (jnxt.eq.0) go to 201
c
         inxt=iuparc(iarc(seg),irow)
         if (inxt.eq.0) go to 201
c
c        ----- evaluate this segment -----
c
         go to
     #   ( 1, 1,50, 1, 5, 1,10, 1, 3, 1, 4,
     #     1, 3, 1, 6, 2, 2,36,
     #     1, 4, 1, 2,11, 2, 9,
     #     1, 4, 1, 2,16, 2,17,
     #     1, 3, 1,16, 2, 2,15), seg
c
    1    continue
            acoef(level+1)=acoef(level)
            goto 120
    2    continue
            acoef(level+1)=-acoef(level)
            goto 120
    3    continue
            junk = jb(jnxt) + 2
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
    4    continue
            junk = jb(jnxt) + 83
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
    5    continue
            junk = jb(jnxt) + 82
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
    6    continue
            junk = jb(jnxt) + 261
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
    9    continue
            junk = jb(jnxt) + 362
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
   10    continue
            junk = jb(jnxt) + 3
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
   11    continue
            junk = jb(jnxt) + 263
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
   15    continue
            acoef(level+1)=acoef(level)*cfs(jb(jnxt)+383)
            go to 120
   16    continue
            acoef(level+1)=acoef(level)*cfs(jb(jnxt)+262)
            go to 120
   17    continue
            acoef(level+1)=acoef(level)*cfs(jb(jnxt)+363)
            go to 120
   36    continue
            junk = jb(jnxt) + 384
            acoef(level+1) = acoef(level) * cfs(junk)
            go to 120
   50    continue
            acoef(level+1) = acoef(level) + acoef(level)
            go to 120
  120    continue
c
c        ----- update stacks -----
c
         newpag=nxtpag(seg)
         if (newpag.eq.0) then
c
c           ----- have finished a loop, check for being at top
c                 of graph
c
            if (level+1.eq.nlevs) then
               number=number+1
               if (number.gt.binsiz) then
                  call lnkerr('bin size to small in top')
               end if
               ires(1,number)=iwtsv(level)+iupwt(iarc(seg),irow)
               ires(2,number)=jwtsv(level)+jupwt(jarc(seg),jrow)
               ires(3,number)=level
               res(4,number)=acoef(level+1)
               go to 201
            end if
c
c           ----- search for heads to top of graph -----
c
            irowsv(level+1)=irow
            jrowsv(level+1)=jrow
            iwtsv(level+1)=iwtsv(level)+iupwt(iarc(seg),irow)
            jwtsv(level+1)=jwtsv(level)+jupwt(jarc(seg),jrow)
            hlev=level+1
            irow=inxt
            jrow=jnxt
            hseg=0
c
c           ----- search up head segments to top of graph -----
c
  401       continue
               hseg=hseg+1
               if (hseg.gt.4) then
c
c                 ----- exhausted arcs at this level, back up one -----
c
                  hlev=hlev-1
                  if (hlev.le.level) then
                     irow=irowsv(level+1)
                     jrow=jrowsv(level+1)
                     go to 201
                  end if
c
                  hseg=segsv(hlev+1)
                  irow=irowsv(hlev+1)
                  jrow=jrowsv(hlev+1)
                  go to 401
               end if
c
c              ----- check that this segment leads down on both i and j
c
               jnxt=juparc(hseg,jrow)
               if (jnxt.eq.0) go to 401
c
               inxt=iuparc(hseg,irow)
               if (inxt.eq.0) go to 401
c
c              ----- check for reaching top of graph -----
c
               if (hlev+1.eq.nlevs) then
                  number=number+1
                  if (number.gt.binsiz) then
                     call lnkerr('bin size to small in top')
                  end if
                  ires(1,number)=iwtsv(hlev)+iupwt(hseg,irow)
                  ires(2,number)=jwtsv(hlev)+jupwt(hseg,jrow)
                  ires(3,number)=level
                  res(4,number)=acoef(level+1)
                  go to 401
               end if
c
c              ----- update stacks and descend -----
c
               irowsv(hlev+1)=irow
               jrowsv(hlev+1)=jrow
               iwtsv(hlev+1)=iwtsv(hlev)+iupwt(hseg,irow)
               jwtsv(hlev+1)=jwtsv(hlev)+jupwt(hseg,jrow)
               segsv(hlev+1)=hseg
               irow=inxt
               jrow=jnxt
               hlev=hlev+1
               hseg=0
            go to 401
         end if
c
c        ----- update stacks and then go down one level -----
c
         if (level+1.ge.nlevs) go to 201
c
         irowsv(level+1)=irow
         jrowsv(level+1)=jrow
         pagesv(level+1)=page
         iwtsv(level+1)=iwtsv(level)+iupwt(iarc(seg),irow)
         jwtsv(level+1)=jwtsv(level)+jupwt(jarc(seg),jrow)
         segsv(level+1)=seg
         irow=inxt
         jrow=jnxt
         page=newpag
c
         level=level+1
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
