*deck @(#)extern.f	5.1  11/6/94
      subroutine extern(orbfrm,nsym,a,b,s,levpt,levnr,arc,wt,
     #                  wtab,wtw,wtx,wty,nrows,nlevs)
c
      implicit integer (a-z)
c
      integer a(nrows),b(nrows),s(nrows),levpt(nlevs),levnr(nlevs)
      integer arc(4,nrows),wt(4,nrows),wtab(orbfrm),wtw(orbfrm,nsym)
      integer wtx(orbfrm,nsym),wty(orbfrm)
c
c     ----- form external-weight arrays -----
c
      call izero(wtab,orbfrm)
      call izero(wtw,orbfrm*nsym)
      call izero(wtx,orbfrm*nsym)
      call izero(wty,orbfrm)
c
      do 10 lev=2,orbfrm+1
         levm1=lev-1
         ptm1=levpt(levm1)
         do 9 row=levpt(lev)+1,levpt(lev)+levnr(lev)
            ia=a(row)
            ib=b(row)
            if (2*ia+ib.gt.2.or.(ia.eq.0.and.ib.eq.0)) go to 9
            is=s(row)+1
            do 8 case=1,4
               rowm1=arc(case,row)
               if (rowm1.gt.0) then
                  if (ia.eq.1.and.case.eq.4) then
                      wtab(levm1)=wt(case,row)
                  endif
                  if (ia.eq.1.and.case.eq.3) then
                      wtw(levm1,is)=wt(case,row)
                  endif
                  if (ib.eq.2.and.case.eq.2) then
                      wtx(levm1,is)=wt(case,row)
                  endif
                  if (ib.eq.1.and.case.eq.2) then
                      wty(levm1)=wt(case,row)
                  endif
               end if
    8       continue
    9    continue
   10 continue
c
c
      return
      end
