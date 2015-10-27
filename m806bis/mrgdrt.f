*deck @(#)mrgdrt.f	1.2  7/30/91
      subroutine mrgdrt(out,nlevs,nrowmx,nrows,nrows1,nrows2,
     #                  a,b,s,x,arc,levpt,levnr,
     #                  a1,b1,s1,x1,arc1,levpt1,levnr1,
     #                  a2,b2,s2,x2,arc2,levpt2,levnr2,
     #                  torow1,torow2,chkx)
c
c***begin prologue  mrgdrt
c***date written   841213   (yymmdd)
c***revision date  yymmdd   (yymmdd)
c***keywords  drt,distinct row table,unitary group ci,
c             configuration interaction
c***author  saxe,paul, (lanl)
c***purpose  to merge to distinct row tables (drts) into a single
c            drt which contains all the configurations in both
c            the original drts.
c***description
c
c   scalars input:
c
c      out     integer
c              unit for output messages, warnings, etc.
c
c      nlevs   integer
c              number of levels in the drts
c
c      nrowmx  integer
c              the maximum nmber of rows permissable in the merged drt
c
c      nrows1  integer
c              the number of rows in the first input drt
c
c      nrows2  integer
c              the number of rows in the secon input drt
c
c   arrays on input:
c
c      a1      integer (nrows1)
c      a2      integer (nrows2)
c              the 'a' values of the rows of the input drts
c
c      b1      integer (nrows1)
c      b2      integer (nrows2)
c              the 'b' values of the rows of the input drts
c
c      s1      integer (nrows1)
c      s2      integer (nrows2)
c              the symmetry associated with the rows of the input drts
c
c      x1      integer (nrows1)
c      x2      integer (nrows2)
c              the excitation value associated with the rows of the
c               input drts.
c
c      arc1    integer (4,nrows1)
c      arc2    integer (4,nrows2)
c              the array of connections between rows on adjacent levels
c
c      levpt1  integer (nlevs)
c      levpt2  integer (nlevs)
c              pointers to the rows at each level of the graph
c
c      levnr1  integer (nlevs)
c      levnr2  integer (nlevs)
c              the number of rows in each level of the input drts
c
c   scalars returned:
c
c      nrows   integer
c              the number of rows in the merged drt
c
c   arrays returned:
c
c      a       integer (nrowmx)
c              the 'a' values of the rows of the merged drt
c
c      b       integer (nrowmx)
c              the 'b' values of the rows of the merged drt
c
c      s       integer (nrowmx)
c              the symmetries associated with the rows of the merged drt
c
c      x       integer (nrowmx)
c              the excitation values associated with the rows of the
c               merged drt
c
c      arc     integer (4,nrowmx)
c              the array of connections between rows of adjacent levels
c               in the merged drt
c
c      levpt   integer (nlevs)
c              the offset of the rows in each level of the merged drt
c
c      levnr   integer (nlevs)
c              the number of rows in each level of the merged drt
c
c   scratch arrays:
c
c      torow1  integer (nrowmx)
c      torow2  integer (nrowmx)
c              used to keep which rows in the input drts a row in the
c               merged drt corresponds to
c
c***routines called  pakdrt
c***end prologue  mrgdrt
c
      implicit integer (a-z)
c
      integer a(nrowmx),b(nrowmx),s(nrowmx),x(nrowmx),arc(4,nrowmx)
      integer levpt(nlevs),levnr(nlevs)
      integer a1(nrows1),b1(nrows1),s1(nrows1),x1(nrows1),arc1(4,nrows1)
      integer levpt1(nlevs),levnr1(nlevs)
      integer a2(nrows2),b2(nrows2),s2(nrows2),x2(nrows2),arc2(4,nrows2)
      integer levpt2(nlevs),levnr2(nlevs)
      integer torow1(nrowmx),torow2(nrowmx)
      logical chkx
c
      nr=0
      do 5 level=nlevs,1,-1
         pt1=levpt1(level)
         nr1=levnr1(level)
         if (level.gt.1) then
            nr1m1=levnr1(level-1)
         else
            nr1m1=0
         end if
         pt2=levpt2(level)
         nr2=levnr2(level)
         levpt(level)=nr
         do 2 row=pt1+1,pt1+nr1
            nr=nr+1
            a(nr)=a1(row)
            b(nr)=b1(row)
            s(nr)=s1(row)
            x(nr)=x1(row)
            arc(1,nr)=arc1(1,row)
            arc(2,nr)=arc1(2,row)
            arc(3,nr)=arc1(3,row)
            arc(4,nr)=arc1(4,row)
    2    continue
         do 4 row=pt2+1,pt2+nr2
            nr=nr+1
            a(nr)=a2(row)
            b(nr)=b2(row)
            s(nr)=s2(row)
            x(nr)=x2(row)
            do 3 case=1,4
               if (arc2(case,row).eq.0) then
                  arc(case,nr)=0
               else
                  arc(case,nr)=arc2(case,row)+nr1m1
               end if
    3       continue
    4    continue
         levnr(level)=nr-levpt(level)
    5 continue
c
      nrows=nr
c
c     ----- eleminate redundant rows -----
c
      call pakdrt(nlevs,nrows,a,b,s,x,arc,levpt,levnr,torow2,chkx)
c
c
      return
      end
