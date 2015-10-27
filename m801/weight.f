*deck @(#)weight.f	5.1  11/6/94
      subroutine weight(a,b,s,arc,levnr,levpt,nlwks,wtab,wtw,wtx,wty
     #,                 wght)
c
c***********************************************************************
c   compute the arc-weight array, starting from the bottom of the      *
c   graph and working up. the algorithm is simple -- the weight of an  *
c   arc is the sum of the weight of the first arc to its left from the *
c   same--same--upper point and the weight (number of lower walks) of  *
c   the point that next arc to the left ends at as a lower point. thus *
c   all vertical arcs (case=1) have a weight of zero. then the weights *
c   of arcs in the external space are transferred to external-weight-  *
c   arrays for use in the external portion of the calculation. these   *
c   arrays are:                                                        *
c       wtab(n)    wt of doubly occupied arc from w point at level n   *
c                  (case=4)                                            *
c       wtw(n,sym) wt of beta walk (case=3) from w point of symmetry   *
c                  sym at level n of graph                             *
c       wtx(n,sym) wt of alpha walk (case=2) from x point of symmetry  *
c                  sym at level n .                                    *
c       wty(n)     wt of alpha walk (case=2) from y point at level n   *
c***********************************************************************
c
      implicit integer (a-z)
      integer numint
c
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /tapes/  out,errout,input,drttap
c
      dimension a(nrows),b(nrows),s(nrows),nlwks(nrows)
      dimension arc(nrows4),levpt(nrows),levnr(nrows)
      dimension wtab(orbfrm),wtw(orbfrm,nsym),wtx(orbfrm,nsym)
      dimension wty(orbfrm),wght(nrows4)
c
c
c     ----- generate the arc-weight array -----
c
      do 32 lev=2,nlevs
         levm1=lev-1
         pontm1=levpt(levm1)
         do 31 row=levpt(lev)+1,levpt(lev)+levnr(lev)
            wt=0
            do 30 case=1,4
               rowm1=arc((row-1)*4+case)
               if (rowm1.gt.0) go to 28
               go to 29
   28          continue
               wght((row-1)*4+case)=wt
               wt=wt+nlwks(pontm1+rowm1)
   29          continue
   30       continue
   31    continue
   32 continue
c
c     ----- generate external-space weight arrays -----
c
      do 72 lev=1,orbfrm
         wtab(lev)=-999999
         wty(lev)=-999999
         do 72 sym=1,nsym
            wtw(lev,sym)=-999999
            wtx(lev,sym)=-999999
   72 continue
c
      if (lev.lt.2) go to 18
      do 17 lev=2,levfrm
         levm1=lev-1
         pontm1=levpt(levm1)
         do 16 row=levpt(lev)+1,levpt(lev)+levnr(lev)
         ia=a(row)
         ib=b(row)
         if (2*ia+ib.gt.2) go to 16
         is=s(row)+1
         do 15 case=1,4
            rowm1=arc((row-1)*4+case)
            if (rowm1.gt.0) go to 13
            go to 14
   13       continue
            if (ia.eq.1.and.case.eq.4) wtab(levm1)=wght((row-1)*4+case)
            if (ia.eq.1.and.case.eq.3) wtw(levm1,is)=wght((row-1)*4+case
     *         )
            if (ib.eq.2.and.case.eq.2) wtx(levm1,is)=wght((row-1)*4+case
     *         )
            if (ib.eq.1.and.case.eq.2) wty(levm1)=wght((row-1)*4+case)
   14       continue
   15    continue
   16    continue
   17 continue
c
   18 continue
      return
      end
