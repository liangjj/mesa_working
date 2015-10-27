*deck @(#)pqx.f	1.2  7/30/91
      subroutine pqx(nrows,nlevs,arc,x,prowsv,qrowsv,arcsv,
     #               levpt)
c
      implicit integer (a-z)
c
      integer arc(4,nrows),x(nrows),prowsv(nlevs),qrowsv(nlevs)
      integer arcsv(nlevs),levpt(nlevs)
c
c     ----- loop down p walks, giving similar q parttial walks x values of -1
c
      x(1)=-1
      prow=2
      qrow=1
      level=nlevs
      harc=0
    1 continue
         harc=harc+1
         if (harc.gt.4) then
            level=level+1
            if (level.gt.nlevs) go to 100
            prow=prowsv(level)
            qrow=qrowsv(level)
            harc=arcsv(level)
            go to 1
         end if
c
         pnxt=arc(harc,prow)
         if (pnxt.le.0) go to 1
c
         qnxt=arc(harc,qrow)
         if (qnxt.le.0) go to 1
c
         pnxt=pnxt+levpt(level-1)
         qnxt=qnxt+levpt(level-1)
c
         if (level.le.2) go to 1
         prowsv(level)=prow
         qrowsv(level)=qrow
         arcsv(level)=harc
         x(qnxt)=-1
         level=level-1
         harc=0
         prow=pnxt
         qrow=qnxt
      go to 1
c
c
  100 continue
c
c
      return
      end
