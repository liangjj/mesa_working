*deck @(#)reord.f	5.1  11/6/94
      subroutine reord(dnarc,upwt,dnwt,dnnwks,offset,nrows,nlevs,
     #                  irowsv,wtsv,dnwtsv,segsv,levdiv,
     #                  dndiag,diag,nwks)
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      real*8 diag(nwks)
c
c     ----- external arrays not modified -----
c
      real*8 dndiag(nwks)
      integer dnarc(4,nrows),upwt(4,nrows),dnwt(4,nrows)
      integer dnnwks(nrows),offset(nrows)
c
c     ----- external arrays used as scratch -----
c
      integer irowsv(nlevs),wtsv(nlevs),dnwtsv(nlevs),segsv(nlevs)
c
c     ----- descend from top of graph to division level -----
c
      hrow=1
      hlev=nlevs
      wtsv(hlev)=0
      dnwtsv(hlev)=0
      harc=0
c
c     ----- start tree search -----
c
  101 continue
         harc=harc+1
         if (harc.gt.4) then
            hlev=hlev+1
            if (hlev.gt.nlevs) go to 9000
            hrow=irowsv(hlev-1)
            harc=segsv(hlev-1)
            go to 101
         end if
c
c        ----- check that this segment exists on the graph -----
c
         hnxt=dnarc(harc,hrow)
         if (hnxt.eq.0) go to 101
c
c        ----- process if at correct level -----
c
         if (hlev-1.eq.levdiv) then
            wt=dnnwks(hnxt)*(wtsv(hlev)+upwt(harc,hnxt))+
     #            offset(hnxt)
            wtdn=dnwtsv(hlev)+dnwt(harc,hrow)
            do 200 i=1,dnnwks(hnxt)
               wt=wt+1
               wtdn=wtdn+1
               diag(wt)=dndiag(wtdn)
  200       continue
         else
            hlev=hlev-1
            irowsv(hlev)=hrow
            segsv(hlev)=harc
            wtsv(hlev)=wtsv(hlev+1)+upwt(harc,hnxt)
            dnwtsv(hlev)=dnwtsv(hlev+1)+dnwt(harc,hrow)
            harc=0
            hrow=hnxt
         end if
c
      go to 101
c
 9000 continue
c
c
      return
      end
