*deck @(#)squeez.f	5.1  11/6/94
      subroutine squeez(bp,sp,arcp,nlwksp,b,s,arc,nlwks)
c
c**********************************************************************
c     pack down the drt arrays now that invalid arcs and points
c     have been eliminate.
c**********************************************************************
c
      implicit integer (a-z)
      integer numint
c
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
c
      dimension bp(nrowsp),sp(nrowsp),arcp(nrows4p),nlwksp(nrowsp)
      dimension b(nrows),s(nrows),arc(nrows4),nlwks(nrows)
c
c
      do 1 i=1,nrows
         b(i)=bp(i)
    1 continue
      do 2 i=1,nrows
         s(i)=sp(i)
    2 continue
      do 3 i=1,nrows4
         arc(i)=arcp(i)
    3 continue
      do 4 i=1,nrows
         nlwks(i)=nlwksp(i)
    4 continue
c
c
      return
      end
