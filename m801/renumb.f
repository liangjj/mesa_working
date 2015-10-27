*deck @(#)renumb.f	5.1  11/6/94
      subroutine renumb(levpt,levnr,arc)
c
      implicit integer (a-z)
      integer numint
c
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
c
      dimension levpt(nlevs),levnr(nlevs),arc(nrows4)
c
c
c     ----- change arc array to incorporate levpt offsets -----
c           so dont have to add levpt(levm1) all the time
c
      do 3 lev=2,nlevs
         pontm1=levpt(lev-1)
         do 2 row=levpt(lev)+1,levpt(lev)+levnr(lev)
            do 1 case=(row-1)*4+1,(row-1)*4+4
               if (arc(case).gt.0) arc(case)=arc(case)+pontm1
    1       continue
    2    continue
    3 continue
      return
      end
