*deck @(#)prtdrt.f	5.1  11/6/94
      subroutine prtdrt(a,b,s,x,nlwks,arc,levpt,levnr,wt,
     #                  nlevs,nrows,out,orbsym,ijgrp,ijadd,kadd,ladd,
     #                  inint,ijxx,ijww,klxx,klww,nklxx,nklww,nnpvir,
     #                  orbfrm,wtab,wtw,wtx,wty,nnp,nsym,norbs)
c
      implicit integer (a-z)
c
      integer a(nrows),b(nrows),s(nrows),x(nrows),arc(4,nrows)
      integer levpt(nlevs),levnr(nlevs),nlwks(nrows),wt(4,nrows)
      integer orbsym(norbs),ijgrp(nnp),ijadd(nnp),kadd(norbs,nsym)
      integer ladd(norbs,nsym),inint(norbs),ijxx(nnp),ijww(nnp)
      integer klxx(nnpvir),klww(nnpvir),nklxx(orbfrm,nsym)
      integer nklww(orbfrm,nsym),wtab(orbfrm),wtw(orbfrm,nsym)
      integer wtx(orbfrm,nsym),wty(orbfrm)
c
      write (out,1)
    1 format (//,'      ***** distinct row table *****',//,
     #   '  row     a  b    s    x  nlwks       arcs       weights',/)
      do 10 level=nlevs,1,-1
         write (out,4)
    4    format (/)
         do 3 row=1,levnr(level)
            i=row+levpt(level)
            write (out,2) row,a(i),b(i),s(i),x(i),nlwks(i),
     #                    (arc(j,i),j=1,4),(wt(k,i),k=1,4)
    2       format (1x,i3,4x,2i3,2i5,i8,3x,4i3,3x,4i5)
    3    continue
   10 continue
c
      call iprint(ijadd,norbs,'ijadd')
      call iprint(ijgrp,norbs,'ijgrp')
c
      write (out,11)
   11 format (/,' inint      kadd and        ladd ',/)
      do 13 i=1,norbs
         write (out,12) i,inint(i),(kadd(i,j),ladd(i,j),j=1,nsym)
   12    format (1x,i3,i8,8x,(2i8,3x))
   13 continue
c
      call iprint(ijxx,norbs,'ijxx')
      call iprint(ijww,norbs,'ijww')
      call iprint(klxx,orbfrm,'klxx')
      call iprint(klww,orbfrm,'klww')
c
      write (out,14)
   14 format (/,'  nklxx and nklww',/)
      do 16 i=1,orbfrm
         write (out,15) i,(nklxx(i,j),nklww(i,j),j=1,nsym)
   15    format (1x,i4,' |',(2i4,3x))
   16 continue
c
      write (out,17)
   17 format (/,' external weight arrays: ab, y, x and w',/)
      do 19 i=1,orbfrm
         write (out,18) i,wtab(i),wty(i),(wtx(i,j),wtw(i,j),j=1,nsym)
   18    format (1x,i3,' | ',2i7,(2i7,3x))
   19 continue
c
      return
      end
