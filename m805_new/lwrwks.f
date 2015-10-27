*deck @(#)lwrwks.f	5.1  11/6/94
      subroutine lwrwks(a,b,s,x,nlwks,arc,levpt,levnr,nlevs,nrows,nwks,
     #                  out)
c
c***begin prologue  lwrwks
c***date written   841212   (yymmdd)
c***revision date  yymmdd   (yymmdd)
c***keywords  drt,distinct row table,unitary group ci,
c             configuration interaction
c***author  saxe,paul, (lanl)
c***purpose  given a rough distinct row table (drt), this routine will
c            generate the number-of-lower-walks (nlwks) from each row,
c            discard useless rows, and compress the drt.
c***description
c
c   scalar input:
c
c     nlevs   integer
c             the number of levels in the drt
c
c     nrows   integer
c             the number of rows in the drt
c
c     out     integer
c             unit for output messages, warnings, etc.
c
c   array input:
c
c     a       integer (nrows)
c             the 'a' values of the rows
c
c     b       integer (nrows)
c             the 'b' values of the rows
c
c     s       integer (nrows)
c             the symmetry (0, 1, ....) associated with the rows
c
c     x       integer (nrows)
c             the excitation values associated with the rows
c
c     arc     integer (4,nrows)
c             the array of connections between rows
c
c     levpt   integer (nlevs)
c             the offset of the rows for a level
c
c     levnr   integer (nlevs)
c             the number of rows in a level
c
c   scalars returned:
c
c     nrows   integer
c             the updated number of rows in the drt
c
c     nwks    integer
c             the number of walks in the graph
c
c   arrays returned:
c
c     a, b, s, x, arc, levpt, levnr are as on input, but modified to
c                                   reflect the compression of the drt
c
c     nlwks   integer (nrows)
c             the number-of-lower-walks from each row
c
c***routines called  (none)
c***end prologue  lwrwks
c
      implicit integer (a-z)
c
      integer a(nrows),b(nrows),s(nrows),x(nrows),nlwks(nrows)
      integer arc(4,nrows),levpt(nlevs),levnr(nlevs)
c
c     ----- find the one-and-only true bottom to the graph -----
c
      nroot=0
      do 1 row=levpt(1)+1,levpt(1)+levnr(1)
         nlwks(row)=0
         if (a(row).ne.0.or.b(row).ne.0.or.s(row).ne.0) go to 1
         nroot=nroot+1
         nlwks(row)=1
    1 continue
c
      if (nroot.gt.1) then
         write (out,2) nroot
    2    format (//,' ##### drt: lwrwks -- invalid number of bottoms ',
     #           'to the shavitt graph:',i3,//)
         call lnkerr('invalid number of bottoms to graph')
      else if (nroot.le.0) then
         nwks=0
         return
      end if
c
c     ----- generate the number-of-lower-walks from each row -----
c
      do 5 level=2,nlevs
         ptlvm1=levpt(level-1)
         do 4 row=levpt(level)+1,levpt(level)+levnr(level)
            nlwk=0
            do 3 case=1,4
               if (arc(case,row).gt.0) nlwk=nlwk+nlwks(arc(case,row)+
     #                                                   ptlvm1)
    3       continue
            nlwks(row)=nlwk
    4    continue
    5 continue
c
      nwks=nlwks(1)
c
c     ----- remove all rows with zero weights -----
c
      pt=levnr(nlevs)
      do 12 level=nlevs-1,1,-1
         nrow=levnr(level)
         ptlv=levpt(level)
         levpt(level)=pt
         do 11 row=ptlv+1,ptlv+nrow
            if (nlwks(row).le.0) then
c
c              ----- discard the row -----
c
               do 7 rowp1=levpt(level+1)+1,levpt(level+1)+levnr(level+1)
                  do 6 case=1,4
                     if (arc(case,rowp1).eq.row-ptlv) arc(case,rowp1)=0
    6             continue
    7          continue
            else
c
c              ----- move the row -----
c
               pt=pt+1
               a(pt)=a(row)
               b(pt)=b(row)
               s(pt)=s(row)
               x(pt)=x(row)
               nlwks(pt)=nlwks(row)
               arc(1,pt)=arc(1,row)
               arc(2,pt)=arc(2,row)
               arc(3,pt)=arc(3,row)
               arc(4,pt)=arc(4,row)
               do 9 rowp1=levpt(level+1)+1,levpt(level+1)+levnr(level+1)
                  do 8 case=1,4
                     if (arc(case,rowp1).eq.row-ptlv) arc(case,rowp1)=
     #                                                pt-levpt(level)
    8             continue
    9          continue
            end if
   11    continue
         levnr(level)=pt-levpt(level)
   12 continue
c
      nrows=pt
c
c
      return
      end
