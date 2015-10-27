*deck @(#)renumb.f	1.2  7/30/91
      subroutine renumb(levpt,levnr,arc,nlevs,nrows)
c
      implicit integer (a-z)
c
      integer levpt(nlevs),levnr(nlevs),arc(4,nrows)
c
      do 3 lev=2,nlevs
         ptm1=levpt(lev-1)
         do 2 row=levpt(lev)+1,levpt(lev)+levnr(lev)
            do 1 case=1,4
               if (arc(case,row).gt.0) arc(case,row)=arc(case,row)+ptm1
    1       continue
    2    continue
    3 continue
c
      return
      end
