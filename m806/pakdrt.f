*deck @(#)pakdrt.f	5.1  11/6/94
      subroutine pakdrt(nlevs,nrows,a,b,s,x,arc,levpt,levnr,newrow,chk)
c
c
      implicit integer (a-z)
c
      integer a(nrows),b(nrows),s(nrows),x(nrows),arc(4,nrows)
      integer levpt(nlevs),levnr(nlevs),newrow(nrows)
      logical chk,chkx
c
c     ----- pass through each level, sequentially numbering the
c           unique rows
c
      chkx=.not.chk
      do 1 row=1,nrows
         newrow(row)=0
    1 continue
      nr=0
      do 40 row=levpt(1)+levnr(1),levpt(1)+1,-1
         if (newrow(row).lt.0) go to 40
         nr=nr+1
         newrow(row)=nr
         do 39 i=row-1,levpt(1)+1,-1
            if (a(i).eq.a(row).and.b(i).eq.b(row).and.
     #        s(i).eq.s(row)) newrow(i)=-nr
cps     #        s(i).eq.s(row).and.(chkx.or.x(i).eq.x(row))) newrow(i)=-nr
   39    continue
   40 continue
c
      do 6 level=2,nlevs
         pt=levpt(level)
         ptm1=levpt(level-1)
         nrow=levnr(level)
c
c        ----- redo the arcs to the new numbering scheme -----
c
         do 3 row=pt+1,pt+nrow
            do 2 case=1,4
               if (arc(case,row).gt.0) arc(case,row)=
     #                abs(newrow(arc(case,row)+ptm1))
    2       continue
    3    continue
c
c        ----- now find the unique rows -----
c
         do 5 row=pt+nrow,pt+1,-1
            if (newrow(row).lt.0) go to 5
            nr=nr+1
            newrow(row)=nr
            do 4 i=row-1,pt+1,-1
               if (a(i).eq.a(row).and.b(i).eq.b(row).and.
     #             s(i).eq.s(row).and.arc(1,i).eq.arc(1,row).and.
     #             arc(2,i).eq.arc(2,row).and.arc(3,i).eq.arc(3,row)
     #             .and.arc(4,i).eq.arc(4,row)
     #              .and.(chkx.or.x(i).eq.x(row)))
     #               newrow(i)=-nr
    4       continue
    5    continue
    6 continue
c
c     ----- now reverse the new order and compact the drt -----
c
      offset=nr+1
      do 20 level=nlevs,1,-1
         nr=0
         pt=levpt(level)
         nrow=levnr(level)
         do 10 row=pt+1,pt+nrow
            if (newrow(row).gt.0) then
               i=offset-newrow(row)
               nr=nr+1
               a(i)=a(row)
               b(i)=b(row)
               s(i)=s(row)
               x(i)=x(row)
               do 7 case=1,4
                  if (arc(case,row).eq.0) then
                     arc(case,i)=0
                  else
                     arc(case,i)=offset-arc(case,row)
                  end if
    7          continue
            else
               i=offset+newrow(row)
c               x(i)=min(x(i),x(row))
            end if
   10    continue
         levnr(level)=nr
   20 continue
c
c     ----- sort out the pointers -----
c
      nr=0
      do 22 level=nlevs,1,-1
         levpt(level)=nr
         nr=nr+levnr(level)
   22 continue
      nrows=nr
c
c     ----- and the arc array -----
c
      do 25 level=nlevs,2,-1
         pt=levpt(level)
         nrow=levnr(level)
         ptm1=levpt(level-1)
         do 24 row=pt+1,pt+nrow
            do 23 case=1,4
               if (arc(case,row).gt.0) arc(case,row)=arc(case,row)-ptm1
   23       continue
   24    continue
   25 continue
c
c
      return
      end
