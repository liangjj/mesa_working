*deck @(#)print.f	5.1  11/6/94
      subroutine print(ops,a,b,s,arc,nlwks,wght,wtab,wtw,wtx,wty,
     $     levnr,levpt,kadd,ladd,inint,ijgrp,ijadd,orbtbf,bfcode,bfkey,
     $     bfsym)
c
c     option(1) controls printing in the following fashion--single bits
c               set control each print section as follows
c
c     bit  printed output
c      1    distinct row table, including weights and chaining indices
c      2    external weight arrays
c      4    integral addressing information
c
c
      implicit integer (a-z)
      integer numint
      logical logkey
      logical all
      character*1 multrf,valenc,bfkey
      character*3 codes,words*18
      character*(*) ops
c
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /tapes/  out,errout,input,drttap
c
      dimension orbtbf(norbs),bfcode(nrefs,nbf),bfkey(nbf),bfsym(nbf)
      dimension a(nrows),b(nrows),s(nrows),arc(nrows4),nlwks(nrows)
      dimension levnr(nlevs),levpt(nlevs)
      dimension wtab(orbfrm),wtw(orbfrm,nsym),wtx(orbfrm,nsym)
      dimension wty(orbfrm),kadd(symorb),ladd(symorb),inint(norbs)
      dimension ijgrp(numij),ijadd(numij),wght(nrows4)
c
   92 format(5x,'number of orbitals allowed in ci',i6)
c
c
      write (out,92) norbs
c
      all=logkey(ops,'print=drt=all',.false.,' ')
c
      if (all.or.logkey(ops,'print=drt=table',.false.,' ')) then
      write (out,19)
   19 format (///,t30,'*** distinct row table ***',/)
      do 6 junk=1,nlevs
         lev=nlevs-junk+1
         orb=lev
         if (lev.eq.nlevs) go to 89
         bf=orbtbf(orb)
         write (out,20) orb,bf,bfsym(bf),words(bfcode(1,bf))
   20    format (/,' orbital',i4,' basis function',i4,' symmetry'
     #   ,       i2,4x,a18,/)
   89    continue
         if (lev.ge.levfrm) go to 1
         write (out,21)
   21    format
     #   (' row  a  b s    nlwk  a1  a2  a3  a4 w1      w2      w3     w
     *4')
         go to 2
    1    continue
         write (out,22)
   22    format
     #   (' -----------------------internal space-----------------------
     *-'/,' row  a  b s    nlwk  a1  a2  a3  a4 w1      w2      w3     w
     *4')
    2    continue
         do 5 row=1,levnr(lev)
            i=row+levpt(lev)
            j=(i-1)*4
            if (lev.ge.levfrm) go to 3
            write (out,23) row,a(i),b(i),s(i),nlwks(i),arc(j+1)
     #      ,              arc(j+2),arc(j+3),arc(j+4),wght(j+1),wght(j+2
     *      ),              wght(j+3),wght(j+4)
   23       format (1x,i3,2i3,i2,i8,4i4,i3,3i8)
            go to 4
    3       continue
            write (out,23) row,a(i),b(i),s(i),nlwks(i),arc(j+1)
     #      ,              arc(j+2),arc(j+3),arc(j+4)
     #      ,              wght(j+1),wght(j+2),wght(j+3),wght(j+4)
    4       continue
    5    continue
    6 continue
      endif
c
c
      if (all.or.logkey(ops,'print=drt=ewa',.false.,' ')) then
      write (out,26)
   26 format (///,'   external weight arrays',/)
      do 9 lev=1,orbfrm
         write (out,27) wtab(lev),wty(lev),(wtw(lev,sq),sq=1,nsym)
     #   ,              (wtx(lev,sq),sq=1,nsym)
   27    format (1x,18i5)
    9 continue
      endif
c
c
      if (all.or.logkey(ops,'print=drt=ipa',.false.,' ')) then
      write (out,33)
   33 format (///,'  integrals addressing scheme',//,'   i    group '
     #,       ' nints  ijadd',/)
      do 11 junk=1,norbs
         i=norbs-junk+1
         write (out,28) i,ijgrp(i*(i+1)/2),inint(i),ijadd(i*(i-1)/2+1)
   28    format (4i7)
   11 continue
      write (out,29)
   29 format  (///,' kadd and ladd')
      do 12 i=1,norbs
         write (out,30) i,(kadd(j),ladd(j),j=i,symorb,norbs)
   30    format (i5,8(2i6,3x))
   12 continue
c
c
   13 continue
      endif
c
c
      return
      end
