*deck @(#)weight.f	1.2  7/30/91
      subroutine weight(nlevs,nrows,levfrm,levpt,levnr,arc,wt,x,nlwks)
c
c***begin prologue  weight
c***date written   850108   (yymmdd)
c***revision date  yymmdd   (yymmdd)
c***keywords  distinct row table, drt, arc weights
c
c***author  saxe, paul,    (lanl)
c***purpose  to calculate the arc weights for a generalised drt
c
c***description
c
c
c***references
c
c***routines called  (none)
c***end prologue  weight
c
      implicit integer (a-z)
c
      integer levpt(nlevs),levnr(nlevs),arc(4,nrows),wt(4,nrows)
      integer x(nrows),nlwks(nrows)
c
c     ----- form the wieghts for all rows -----
c
      do 4 level=2,nlevs
         pt=levpt(level)
         nr=levnr(level)
         ptm1=levpt(level-1)
         do 3 row=pt+1,pt+nr
            wght=0
            do 2 case=1,4
               rowm1=arc(case,row)
               if (rowm1.gt.0) then
                  wt(case,row)=wght
                  wght=wght+nlwks(ptm1+rowm1)
               else
                  wt(case,row)=0
               end if
    2       continue
    3    continue
    4 continue
c
c     ----- if the configurations are divided into a p and q space,
c           offset the p space to after the q space.
c
c      if (levnr(nlevs).gt.1) then
c         pt=levpt(levfrm+1)
c         nr=levnr(levfrm+1)
c         offset=nlwks(1)
c         do 6 row=pt+1,pt+nr
c            if (x(row).eq.1) then
c               do 5 case=1,4
c                  if (arc(case,row).gt.0) wt(case,row)=wt(case,row)+
c     #                                                      offset
c    5          continue
c            end if
c    6    continue
c      end if
c
c
      return
      end
