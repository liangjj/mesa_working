*deck @(#)squash.f	5.1  11/6/94
      subroutine squash(a,b,s,arc,nlwks,ptsim,nsim,sim,levpt,levnr,
     #                  a1,b1,s1,arc1,levpt1,levnr1,temp,nrows,
     #                  nrows1,nrowmx,nlevs,nwks,out,mxsim,ntemp)
c
      implicit integer (a-z)
c
      integer a(nrowmx),b(nrowmx),s(nrowmx),ptsim(nrowmx),nsim(nrowmx)
      integer sim(mxsim),levpt(nlevs),levnr(nlevs),a1(nrows1),b1(nrows1)
      integer s1(nrows1),arc1(4,nrows1),arc(4,nrowmx),levpt1(nlevs)
      integer levnr1(nlevs),temp(ntemp),nlwks(nrowmx)
c
c     ----- the top level of the graph is one row -----
c
      pts=0
      a(1)=a1(1)
      b(1)=b1(1)
      s(1)=s1(1)
      ptsim(1)=pts
      nsim(1)=levnr1(nlevs)
      do 1 i=1,nsim(1)
         pts=pts+1
         sim(pts)=i
    1 continue
      levpt(nlevs)=0
      levnr(nlevs)=1
c
c     ----- descend through graph, adding next levels -----
c
      do 100 level=nlevs,2,-1
         pt=levpt(level)
         nrow=levnr(level)
         pt1=levpt1(level)
         pt1m1=levpt1(level-1)
         ptm1=pt+nrow
         levpt(level-1)=ptm1
         nrm1=0
         do 10 row=pt+1,pt+nrow
            p=ptsim(row)
            n=nsim(row)
            do 9 case=1,4
               nm1=0
               do 2 i=1,n
                  row1=sim(p+i)+pt1
                  if (arc1(case,row1).eq.0) go to 2
                  do 83 j=1,nm1
                     if (arc1(case,row1).lt.temp(j)) go to 84
   83             continue
                  nm1=nm1+1
                  if (nm1.gt.ntemp)
     #              call lnkerr('drt: squash--temp too small')
                  temp(nm1)=arc1(case,row1)
                  go to 88
   84             continue
                  t=temp(j)
                  temp(j)=arc1(case,row1)
                  do 85 k=j+1,nm1
                     t1=temp(j)
                     temp(j)=t
                     t=t1
   85             continue
                  nm1=nm1+1
                  if (nm1.gt.ntemp)
     #              call lnkerr('drt: squash--temp too small')
                  temp(nm1)=t
   88             continue
    2          continue
               if (nm1.le.0) then
                  arc(case,row)=0
                  go to 9
               end if
c
c              ----- see if this row already added -----
c
               do 5 rowm1=ptm1+1,ptm1+nrm1
                  if (nsim(rowm1).ne.nm1) go to 5
                  pm1=ptsim(rowm1)
                  do 3 i=1,nm1
                     if (sim(pm1+i).ne.temp(i)) go to 5
    3             continue
c
c                 ----- found the row, so add arc -----
c
                  arc(case,row)=rowm1-ptm1
                  go to 9
    5          continue
c
c              ----- this is a new row, so add -----
c
               nrm1=nrm1+1
               rowm1=ptm1+nrm1
               row1m1=pt1m1+temp(1)
               a(rowm1)=a1(row1m1)
               b(rowm1)=b1(row1m1)
               s(rowm1)=s1(row1m1)
               ptsim(rowm1)=pts
               nsim(rowm1)=nm1
               arc(case,row)=nrm1
               if (pts+nm1.gt.mxsim)
     #            call lnkerr('drt: squash--mxsim too small')
               do 6 i=1,nm1
                  pts=pts+1
                  sim(pts)=temp(i)
    6          continue
    9       continue
   10    continue
         levnr(level-1)=nrm1
  100 continue
c
c     ----- add the bottom row -----
c
      pt=levpt(1)+1
      do 101 case=1,4
         arc(case,pt)=0
  101 continue
      nrows=pt
c
c     ----- eliminate redundant rows -----
c
      do 102 i=1,nrows
         ptsim(i)=0
  102 continue
c
      call pakdrt(nlevs,nrows,a,b,s,ptsim,arc,levpt,levnr,nlwks,.false.)
c
c     ----- number the walks -----
c
      call lwrwks(a,b,s,ptsim,nlwks,arc,levpt,levnr,nlevs,nrows,
     #            nwks,out)
c
c
      return
      end
