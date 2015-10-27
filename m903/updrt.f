*deck @(#)updrt.f	5.1  11/6/94
      subroutine updrt(dnarc,uparc,dnnwks,upnwks,dnwt,upwt,levpt,levnr,
     #                 nlevs,nrows)
c
      implicit integer (a-z)
c
      integer dnarc(4,nrows),uparc(4,nrows),dnnwks(nrows),upnwks(nrows)
      integer dnwt(4,nrows),upwt(4,nrows),levpt(nlevs),levnr(nlevs)
c
      common /io/ inp,iout
c
c     ----- create the upward chaining indices from the downward -----
c
      call izero(uparc,4*nrows)
c
      do 3 level=2,nlevs
         do 2 toprow=levpt(level)+1,levpt(level)+levnr(level)
            do 1 case=1,4
               botrow=dnarc(case,toprow)
               if (botrow.gt.0) then
                  uparc(case,botrow)=toprow
               end if
    1       continue
    2    continue
    3 continue
c
c     ----- form the number of upper walks from each row -----
c
      call izero(upnwks,nrows)
      upnwks(levpt(nlevs)+1)=1
c
      do 6 level=nlevs-1,1,-1
         do 5 row=levpt(level)+1,levpt(level)+levnr(level)
            nupwks=0
            do 4 case=1,4
               uprow=uparc(case,row)
               if (uprow.gt.0) then
                  nupwks=nupwks+upnwks(uprow)
               end if
    4       continue
            upnwks(row)=nupwks
    5    continue
    6 continue
c
c     ----- create the weight arrays for upward direction -----
c
      call izero(upwt,4*nrows)
c
      do 9 level=nlevs-1,1,-1
         do 8 row=levpt(level)+1,levpt(level)+levnr(level)
            wt=0
            do 7 case=1,4
               uprow=uparc(case,row)
               if (uprow.gt.0) then
                  upwt(case,row)=wt
                  wt=wt+upnwks(uprow)
               end if
    7       continue
    8    continue
    9 continue
c
cps      do 910 level=nlevs,1,-1
cps         write (iout,901) level
cps  901    format (/,' level=',i3,t15,'row',t20,'downward arcs',
cps     #           t40,'upward arcs',t55,'lower walks',t65,'upper walks')
cps         do 909 row=levpt(level)+1,levpt(level)+levnr(level)
cps            write (iout,902) row,(dnarc(case,row),case=1,4),
cps     #                       (uparc(case,row),case=1,4),dnnwks(row),
cps     #                       upnwks(row)
cps  902       format (t14,i4,t18,4i4,t38,4i4,t55,i6,t65,i7)
cps  909    continue
cps  910 continue
c
c
cps      do 920 level=nlevs,1,-1
cps         write (iout,911) level
cps  911    format (/,' level=',i3,t15,'row',t20,'downward wt',
cps     #           t40,'upward wt',t55,'lower walks',t65,'upper walks')
cps         do 919 row=levpt(level)+1,levpt(level)+levnr(level)
cps            write (iout,902) row,(dnwt(case,row),case=1,4),
cps     #                       (upwt(case,row),case=1,4),dnnwks(row),
cps     #                       upnwks(row)
cps  919    continue
cps  920 continue
c
c
      return
      end
