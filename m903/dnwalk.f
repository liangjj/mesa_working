*deck @(#)dnwalk.f	5.1  11/6/94
      subroutine dnwalk(idnarc,idnwt,ia,ib,inrows,
     #                  jdnarc,      ja,jb,jnrows,
     #                  irowsv,jrowsv,iwtsv,segsv,
     #                  irowst,jrowst,levdiv,nlevs,
     #                  iwalk,jndn)
c
c***begin prologue    dnwalk
c***date written      860915   (yymmdd)
c***revision date     yymmdd   (yymmdd)
c***keywords          ci, coupling coefficients, guga
c***author            saxe, paul,    (lanl)
c***purpose           to find the walk numbers in drt i that correspond
c                     to all the walks in drt j down from row jrowst
c***description
c
c           dnwalk finds the partial downward walks in drt i corresponding
c       to all the downward walks from row jrowst in drt j. irowst is the
c       row in j corresponding to irowst in i. returned is iwalk, with the
c       i walk number for each lower walk in j.
c
c***references
c
c***routines called
c***end prologue      dnwalk
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      integer iwalk(jndn)
c
c     ----- external arrays not modified -----
c
      integer idnarc(4,inrows),idnwt(4,inrows),ia(inrows)
      integer ib(inrows)
      integer jdnarc(4,jnrows),ja(jnrows),jb(jnrows)
c
c     ----- external arrays used for scratch -----
c
      integer irowsv(nlevs),jrowsv(nlevs),iwtsv(nlevs)
      integer segsv(nlevs)
c
c     ----- common blocks -----
c
      common /io/ inp,iout
c
c     ----- check that the starting rows coincide -----
c
      if (ia(irowst).ne.ja(jrowst).or.ib(irowst).ne.jb(jrowst)) then
         call lnkerr('starting rows in dnwalk do not coincide')
      end if
c
c     ----- descend partial walks -----
c
      jwt=0
      irow=irowst
      jrow=jrowst
      lev=levdiv
      iwtsv(lev)=0
      seg=0
c
c     ----- start tree search of head sections -----
c
  101 continue
         seg=seg+1
         if (seg.gt.4) then
            lev=lev+1
            if (lev.gt.levdiv) go to 1000
            irow=irowsv(lev-1)
            jrow=jrowsv(lev-1)
            seg=segsv(lev-1)
            go to 101
         end if
c
c        ----- check that this segment exists on the graph -----
c
         jnxt=jdnarc(seg,jrow)
         if (jnxt.eq.0) go to 101
c
         inxt=idnarc(seg,irow)
         if (inxt.eq.0) go to 101
c
c        ----- if down to bottom of graph, store information -----
c
         if (lev.eq.2) then
            jwt=jwt+1
            iwalk(jwt)=iwtsv(lev)+idnwt(seg,irow)
            go to 101
         end if
c
c        ----- otherwise, descend another level -----
c
         lev=lev-1
         irowsv(lev)=irow
         jrowsv(lev)=jrow
         segsv(lev)=seg
         iwtsv(lev)=iwtsv(lev+1)+idnwt(seg,irow)
         seg=0
         irow=inxt
         jrow=jnxt
      go to 101
c
c     ----- the return is down here for convenience -----
c
 1000 continue
c
c     ----- check that number of jwalks worked out -----
c
      if (jwt.ne.jndn) then
         call lnkerr('in dnwalk, did not get correct number of walks')
      end if
c
      return
      end
