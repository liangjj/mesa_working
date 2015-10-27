*deck @(#)pandq.f	5.1  11/6/94
      subroutine pandq(a,b,s,arc,levpt,levnr,nrows,a1,b1,s1,arc1,levpt1,
     #                 levnr1,nrows1,nrowmx,nlevs,refcod,norbs,nrefs,
     #                 torow1,toref,n,p,uoc,bet,alp,doc,temp,levfrm,
     #                 mxref)
c
      implicit integer (a-z)
c
      integer a(nrowmx),b(nrowmx),s(nrowmx),arc(4,nrowmx),levpt(nlevs)
      integer levnr(nlevs),a1(nrows1),b1(nrows1),s1(nrows1)
      integer arc1(4,nrows1),levpt1(nlevs),levnr1(nlevs)
      integer refcod(norbs,nrefs),torow1(nrowmx),toref(*),n(nrowmx)
      integer p(nrowmx),temp(nrefs)
      integer helper(4)
c
      helper(1)=uoc
      helper(2)=alp
      helper(3)=bet
      helper(4)=doc
c
      ptref=0
      levpt(nlevs)=0
      levnr(nlevs)=2
      a(1)=a1(1)
      b(1)=b1(1)
      s(1)=s1(1)
      arc(1,1)=0
      arc(2,1)=0
      arc(3,1)=0
      arc(4,1)=0
      n(1)=-nrefs
      p(1)=ptref
      torow1(1)=1
      do 1 ref=1,nrefs
         ptref=ptref+1
         toref(ptref)=ref
    1 continue
      a(2)=a(1)
      b(2)=b(1)
      s(2)=s(1)
      arc(1,2)=0
      arc(2,2)=0
      arc(3,2)=0
      arc(4,2)=0
      torow1(2)=1
      n(2)=nrefs
      p(2)=ptref
      do 2 ref=1,nrefs
         ptref=ptref+1
         toref(ptref)=ref
    2 continue
c
      do 1000 level=nlevs,2,-1
         orb=level-1
         pt=levpt(level)
         nr=levnr(level)
         pt1=levpt1(level)
         pt1m1=levpt1(level-1)
         ptm1=pt+nr
         levpt(level-1)=ptm1
         nrm1=0
         do 900 row=pt+1,pt+nr
            row1=torow1(row)
            nn=n(row)
            if (nn.eq.0) then
c
c              ----- plain vanilla row in drt -----
c
               do 50 case=1,4
                  if (arc1(case,row1).le.0) go to 50
                  row1m1=arc1(case,row1)+pt1m1
c
c                 ----- decide if the row we get to exists -----
c
                  do 10 rowm1=ptm1+1,ptm1+nrm1
                     if (torow1(rowm1).eq.row1m1.and.n(rowm1).eq.0)
     #                                                     go to 20
   10             continue
c
c                 ----- this row does not exist, so add it -----
c
                  nrm1=nrm1+1
                  rowm1=ptm1+nrm1
                  a(rowm1)=a1(row1m1)
                  b(rowm1)=b1(row1m1)
                  s(rowm1)=s1(row1m1)
                  n(rowm1)=0
                  p(rowm1)=0
                  torow1(rowm1)=row1m1
                  arc(1,rowm1)=0
                  arc(2,rowm1)=0
                  arc(3,rowm1)=0
                  arc(4,rowm1)=0
c
c                 ----- add the arc to the upper row -----
c
   20             continue
                  arc(case,row)=rowm1-ptm1
   50          continue
            else if (nn.lt.0) then
c
c              ----- if we are dowm to the fermi-level, throw the rows out,
c                    since these walks will be counted as references
c
               if (level.le.levfrm) go to 900
c
c              ----- row looks like part of a reference walk, but isn't
c
               pp=p(row)
               do 150 case=1,4
                  if (arc1(case,row1).le.0) go to 150
                  row1m1=arc1(case,row1)+pt1m1
c
c                 ----- loop through references arriving at row,
c                       determining which go where
c
                  nnm1=0
                  do 110 i=1,-nn
                     ref=toref(pp+i)
                     if (refcod(orb,ref).eq.helper(case)) then
                        nnm1=nnm1+1
                        temp(nnm1)=ref
                     end if
  110             continue
c
c                 ----- scan added rows to see if this one exists -----
c
                  do 120 rowm1=ptm1+1,ptm1+nrm1
                     if (torow1(rowm1).ne.row1m1.or.-n(rowm1).ne.nnm1)
     #                                                        go to 120
                     ppm1=p(rowm1)
                     do 115 i=1,nnm1
                        if (toref(ppm1+i).ne.temp(i)) go to 120
  115                continue
                     go to 125
  120             continue
c
c                 ----- must add the row -----
c
                  nrm1=nrm1+1
                  rowm1=ptm1+nrm1
                  if (rowm1.gt.nrowmx) then
                     call lnkerr('drt: pandq--nrowmx too small')
                  end if
                  a(rowm1)=a1(row1m1)
                  b(rowm1)=b1(row1m1)
                  s(rowm1)=s1(row1m1)
                  n(rowm1)=-nnm1
                  p(rowm1)=ptref
                  torow1(rowm1)=row1m1
                  arc(1,rowm1)=0
                  arc(2,rowm1)=0
                  arc(3,rowm1)=0
                  arc(4,rowm1)=0
                  do 123 i=1,nnm1
                     ptref=ptref+1
                     if (ptref.gt.mxref) then
                        call lnkerr('drt: pandq--mxref too small')
                     end if
                     toref(ptref)=temp(i)
  123             continue
c
c                 ----- add the arc to the upper row -----
c
  125             continue
                  arc(case,row)=rowm1-ptm1
  150          continue
            else if (nn.gt.0) then
c
c               ----- and this row is actually on a reference walk -----
c
                pp=p(row)
                do 250 case=1,4
                   nnm1=0
                   do 210 i=1,nn
                      ref=toref(pp+i)
                      if (refcod(orb,ref).eq.helper(case)) then
                         nnm1=nnm1+1
                         temp(nnm1)=ref
                      end if
  210              continue
                   if (nnm1.le.0) go to 250
c
c                  ----- merge with normal external part of drt at
c                        fermi level.
c
                   if (level.eq.levfrm+1) nnm1=0
c
c                  ----- see if this row exists already -----
c
                   row1m1=arc1(case,row1)+pt1m1
                   do 220 rowm1=ptm1+1,ptm1+nrm1
                      if (torow1(rowm1).ne.row1m1.or.n(rowm1).ne.nnm1)
     #                                                       go to 220
                      ppm1=p(rowm1)
                      do 215 i=1,nnm1
                         if (toref(ppm1+i).ne.temp(i)) go to 220
  215                 continue
                      go to 230
  220              continue
c
c                  ----- add a row -----
c
                  nrm1=nrm1+1
                  rowm1=ptm1+nrm1
                  if (rowm1.gt.nrowmx) then
                     call lnkerr('drt: pandq--nrowmx too small 2')
                  end if
                  a(rowm1)=a1(row1m1)
                  b(rowm1)=b1(row1m1)
                  s(rowm1)=s1(row1m1)
                  n(rowm1)=nnm1
                  p(rowm1)=ptref
                  torow1(rowm1)=row1m1
                  arc(1,rowm1)=0
                  arc(2,rowm1)=0
                  arc(3,rowm1)=0
                  arc(4,rowm1)=0
                  do 225 i=1,nnm1
                     ptref=ptref+1
                     if (ptref.gt.mxref) then
                        call lnkerr('drt: pandq--mxref too small 2')
                     end if
                     toref(ptref)=temp(i)
  225             continue
c
c                 ----- add the arc to the upper row -----
c
  230             continue
                  arc(case,row)=rowm1-ptm1
  250          continue
            end if
  900    continue
         levnr(level-1)=nrm1
 1000 continue
c
      nrows=levpt(1)+levnr(1)
c
c     ----- set up x(n) array properly -----
c
      do 1010 row=1,nrows
         if (n(row).lt.0) then
            n(row)=-1
         else if (n(row).gt.0) then
            n(row)=1
         end if
 1010 continue
c
c     ----- pack the similar rows down -----
c
      call pakdrt(nlevs,nrows,a,b,s,n,arc,levpt,levnr,torow1,.false.)
      return
      end
