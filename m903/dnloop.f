*deck @(#)dnloop.f	5.1  11/6/94
      subroutine dnloop(idnarc,idnwt,ia,ib,idnnwk,inrows,
     #                  jdnarc,jdnwt,ja,jb,jdnnwk,jnrows,
     #                  acoef,irowsv,jrowsv,iwtsv,jwtsv,pagesv,segsv,
     #                  irowst,jrowst,levdiv,nlevs,cfs,
     #                  res,ires,wpti4,binsiz,number,nij,ijpt)
c
c***begin prologue    dnloop
c***date written      860910   (yymmdd)
c***revision date     yymmdd   (yymmdd)
c***keywords          ci, coupling coefficients, guga
c***author            saxe, paul,    (lanl)
c***purpose           to calculate the one-body coupling coefficients
c                     the lower portion of a lexical drt
c***description
c
c           dnloop generates all the one-electron coupling coefficients
c       eminating from rows irowst and jrowst at level levdiv. irowst
c       and jrowst must be the same except for symmetry. the lexical
c       walk numbers are returned in ires(1,*) and (2,*), the coupling
c       coefficient in res(4,*), and the orbital index pair in ires(3,*).
c       number returns the number of loops found.
c
c***references
c
c***routines called
c***end prologue      dnloop
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
      integer ijpt(nij)
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
c
      save segmin,segmax,nxtpag,jarc,iarc
c
      number=0
c
c     ----- descend partial walks to the opening level -----
c
      do 1000 top=levdiv,2,-1
         if (top.eq.levdiv) then
c
c           ----- special situation for loops opening at dividing level
c
            hseg=4
            hlev=levdiv
            level=top
            irow=irowst
            jrow=jrowst
            iwtsv(top)=0
            jwtsv(top)=0
            acoef(top)=1.0d+00
            page=1
            seg=segmin(page)
            segmx=segmax(page)
            go to 201
         end if
c
c        ----- initialize head section ------
c
         irow=irowst
         jrow=jrowst
         hlev=levdiv
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
               if (hlev.gt.levdiv) go to 1000
               irow=irowsv(hlev-1)
               jrow=jrowsv(hlev-1)
               hseg=segsv(hlev-1)
               go to 101
            end if
c
c           ----- check that this segment exists on the graph -----
c
            jnxt=jdnarc(hseg,jrow)
            if (jnxt.eq.0) go to 101
c
            inxt=idnarc(hseg,irow)
            if (inxt.eq.0) go to 101
c
c           ----- hop out of search if descended to desired level -----
c
            if (hlev-1.eq.top) go to 102
c
c           ----- otherwise, descend another level -----
c
            hlev=hlev-1
            irowsv(hlev)=irow
            jrowsv(hlev)=jrow
            segsv(hlev)=hseg
            iwtsv(hlev)=iwtsv(hlev+1)+idnwt(hseg,irow)
            jwtsv(hlev)=jwtsv(hlev+1)+jdnwt(hseg,jrow)
            hseg=0
            irow=inxt
            jrow=jnxt
         go to 101
c
c     ----- initialize stacks for tree search -----
c
  102    continue
         level=top
         irowsv(level)=irow
         jrowsv(level)=jrow
         iwtsv(level)=iwtsv(level+1)+idnwt(hseg,irow)
         jwtsv(level)=jwtsv(level+1)+jdnwt(hseg,jrow)
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
  201    continue
c
            seg=seg+1
            if (seg.gt.segmx) then
c
c              ----- exhausted this page of the table, so back up a level
c
               level=level+1
               if (level.gt.top) then
                  irow=irowsv(top)
                  jrow=jrowsv(top)
                  go to 101
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
c           ----- check that this segment is valid on the graph -----
c
            jnxt=jdnarc(jarc(seg),jrow)
            if (jnxt.eq.0) go to 201
c
            inxt=idnarc(iarc(seg),irow)
            if (inxt.eq.0) go to 201
c
c           ----- evaluate this segment -----
c
            go to
     #      ( 1, 1, 1, 3, 1, 4,50, 1, 1, 3, 4,
     #        1, 6, 1, 9, 2, 5, 2,
     #        1, 1, 2,36,11,10, 2,
     #        1, 1, 2,15,16, 3, 2,
     #        1, 1,16,17, 4, 2, 2), seg
c
    1       continue
               acoef(level-1)=acoef(level)
               goto 120
    2       continue
               acoef(level-1)=-acoef(level)
               goto 120
    3       continue
               junk = jb(jrow) + 2
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
    4       continue
               junk = jb(jrow) + 83
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
    5       continue
               junk = jb(jrow) + 82
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
    6       continue
               junk = jb(jrow) + 261
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
    9       continue
               junk = jb(jrow) + 362
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
   10       continue
               junk = jb(jrow) + 3
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
   11       continue
               junk = jb(jrow) + 263
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
   15       continue
               acoef(level-1)=acoef(level)*cfs(jb(jrow)+383)
               go to 120
   16       continue
               acoef(level-1)=acoef(level)*cfs(jb(jrow)+262)
               go to 120
   17       continue
               acoef(level-1)=acoef(level)*cfs(jb(jrow)+363)
               go to 120
   36       continue
               junk = jb(jrow) + 384
               acoef(level-1) = acoef(level) * cfs(junk)
               go to 120
   50       continue
               acoef(level-1) = acoef(level) + acoef(level)
               go to 120
  120       continue
c
c           ----- update stacks -----
c
            newpag=nxtpag(seg)
            if (newpag.eq.0) then
c
c              ----- have finished a loop, check for being at bottom
c                    of graph
c
               if (level.eq.2) then
                  number=number+1
                  if (number.gt.binsiz) then
                     call lnkerr('bin size to small in dnloop')
                  end if
                  ires(1,number)=iwtsv(level)+idnwt(iarc(seg),irow)
                  ires(2,number)=jwtsv(level)+jdnwt(jarc(seg),jrow)
                  ires(3,number)=ijpt((top-1)*(top-2)/2+level-1)
                  res(4,number)=acoef(level-1)
                  go to 201
               end if
c
c              ----- search for tails to bottom of graph -----
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
c              ----- search down tail segments to bottom of graph -----
c
  401          continue
                  tseg=tseg+1
                  if (tseg.gt.4) then
c
c                    ----- exhausted arcs at this level, back up one -----
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
c                 ----- check that this segment leads down on both i and j
c
                  jnxt=jdnarc(tseg,jrow)
                  if (jnxt.eq.0) go to 401
c
                  inxt=idnarc(tseg,irow)
                  if (inxt.eq.0) go to 401
c
c                 ----- check for reaching bottom of graph -----
c
                  if (tlev.eq.2) then
                     number=number+1
                     if (number.gt.binsiz) then
                        call lnkerr('bin size to small in dnloop')
                     end if
                     ires(1,number)=iwtsv(tlev)+idnwt(tseg,irow)
                     ires(2,number)=jwtsv(tlev)+jdnwt(tseg,jrow)
                     ires(3,number)=ijpt((top-1)*(top-2)/2+level-1)
                     res(4,number)=acoef(level-1)
                     go to 401
                  end if
c
c                 ----- update stacks and descend -----
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
c           ----- update stacks and then go down one level -----
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
 1000 continue
c
c     ----- the return is down here for convenience -----
c
 9000 continue
c
      return
      end
