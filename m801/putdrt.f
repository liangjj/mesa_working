*deck  @(#)putdrt.f	5.1 11/6/94
      subroutine putdrt(bfsym,bfcode,orbtbf,iout,levpt,levnr,a,b,s,arc
     #,                 nlwks,wght,wtab,wtw,wtx,wty,kadd,ladd,ijadd
     #,                 ijgrp,inint,inext,jmnnxt,jmxnxt,ningrp,orbsym
     #,                 ijxx,klxx,nklxx,ijww,klww,nklww,refwt)
c
      implicit integer (a-z)
      integer numint
      real*8 versin
c
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh,val
      common /tapes/  out,errout,input,drttap
      common /dimens/ nbf,nsym,norbs,nrowsp,nrows4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc,nvref
     #,               nijvir
      common /versin/ versin
c
      dimension bfsym(nbf),bfcode(nrefs,nbf),orbtbf(norbs),iout(nbf)
      dimension levpt(nlevs),levnr(nlevs),a(nrows),b(nrows),s(nrows)
      dimension arc(nrows4),nlwks(nrows),wght(nrows4),wtab(orbfrm)
      dimension wtw(orbfrm,nsym),wtx(orbfrm,nsym),wty(orbfrm)
      dimension kadd(symorb),ladd(symorb)
      dimension ijadd(numij),ijgrp(numij),inint(norbs),inext(norbs)
      dimension jmnnxt(norbs),jmxnxt(norbs),ningrp(ngroup),orbsym(norbs)
      dimension ijxx(numij ),klxx(nijvir),nklxx(nsym,orbfrm)
      dimension ijww(numij ),klww(nijvir),nklww(nsym,orbfrm)
      dimension refwt(nvref)
c
c
c     rewind drttap
c     write (drttap) versin
c     write (drttap) label
c     intnum=numint
c     write (drttap) nbf,nsym,norbs,nrows,nrows4
c    #,              nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
c    #,              orbfrm,symorb,numij,ngroup,intnum,nmax,maxb
c    #,              nijvir,nvref
c     write (drttap) kadd
c     write (drttap) ladd
c     write (drttap) ijadd
c     write (drttap) ijgrp
c     write (drttap) inint
c     write (drttap) inext
c     write (drttap) jmnnxt
c     write (drttap) jmxnxt
c     write (drttap) ningrp
c     write (drttap) iout
c     write (drttap) orbsym
c     write (drttap) bfsym
c     write (drttap) bfcode
c     write (drttap) orbtbf
c     write (drttap) levpt
c     write (drttap) levnr
c     write (drttap) a
c     write (drttap) b
c     write (drttap) s
c     write (drttap) arc
c     write (drttap) nlwks
c     write (drttap) wght
c     write (drttap) wtab
c     write (drttap) wtw
c     write (drttap) wtx
c     write (drttap) wty
c     write (drttap) ijxx
c     write (drttap) klxx
c     write (drttap) nklxx
c     write (drttap) ijww
c     write (drttap) klww
c     write (drttap) nklww
c     write (drttap) refwt
c     endfile drttap
      call iosys('write integer nrows to rwf',1,nrows,0,' ')
      call iosys('write integer nlevs to rwf',1,nlevs,0,' ')
      call iosys('write integer norbs to rwf',1,norbs,0,' ')
      call iosys('write integer nrefs to rwf',1,nrefs,0,' ')
      call iosys('write integer nwks to rwf',1,nwks,0,' ')
      call iosys('write integer orbfrm to rwf',1,orbfrm,0,' ')
      call iosys('write integer symorb to rwf',1,symorb,0,' ')
      call iosys('write integer numij to rwf',1,numij,0,' ')
      call iosys('write integer ngroup to rwf',1,ngroup,0,' ')
      call iosys('write integer nmax to rwf',1,nmax,0,' ')
      call iosys('write integer maxb to rwf',1,maxb,0,' ')
      call iosys('write integer nijvir to rwf',1,nijvir,0,' ')
c
      call iosys('write integer kadd to rwf',symorb,kadd,0,' ')
      call iosys('write integer ladd to rwf',symorb,ladd,0,' ')
      call iosys('write integer ijadd to rwf',numij,ijadd,0,' ')
      call iosys('write integer ijgrp to rwf',numij,ijgrp,0,' ')
      call iosys('write integer inint to rwf',norbs,inint,0,' ')
      call iosys('write integer inext to rwf',norbs,inext,0,' ')
      call iosys('write integer jmnnxt to rwf',norbs,jmnnxt,0,' ')
      call iosys('write integer jmxnxt to rwf',norbs,jmxnxt,0,' ')
      call iosys('write integer ningrp to rwf',ngroup,ningrp,0,' ')
      call iosys('write integer iout to rwf',nbf,iout,0,' ')
      call iosys('write integer orbsym to rwf',norbs,orbsym,0,' ')
      call iosys('write integer bfsym to rwf',nbf,bfsym,0,' ')
      call iosys('write integer bfcode to rwf',nrefs*nbf,bfcode,0,' ')
      call iosys('write integer orbtbf to rwf',norbs,orbtbf,0,' ')
      call iosys('write integer levpt to rwf',nlevs,levpt,0,' ')
      call iosys('write integer levnr to rwf',nlevs,levnr,0,' ')
      call iosys('write integer a to rwf',nrows,a,0,' ')
      call iosys('write integer b to rwf',nrows,b,0,' ')
      call iosys('write integer s to rwf',nrows,s,0,' ')
      call iosys('write integer arc to rwf',4*nrows,arc,0,' ')
      call iosys('write integer nlwks to rwf',nrows,nlwks,0,' ')
      call iosys('write integer weight to rwf',4*nrows,wght,0,' ')
      if (orbfrm.gt.0) then
      call iosys('write integer wtab to rwf',orbfrm,wtab,0,' ')
      call iosys('write integer wtw to rwf',orbfrm*nsym,wtw,0,' ')
      call iosys('write integer wtx to rwf',orbfrm*nsym,wtx,0,' ')
      call iosys('write integer wty to rwf',orbfrm,wty,0,' ')
      call iosys('write integer ijxx to rwf',numij,ijxx,0,' ')
      call iosys('write integer klxx to rwf',nijvir,klxx,0,' ')
      call iosys('write integer nklxx to rwf',nsym*orbfrm,nklxx,0,' ')
      call iosys('write integer ijww to rwf',numij,ijww,0,' ')
      call iosys('write integer klww to rwf',nijvir,klww,0,' ')
      call iosys('write integer nklww to rwf',nsym*orbfrm,nklww,0,' ')
      end if
cps      call iosys('write integer ref_weight to rwf',nvref,refwt,0,' ')
c
c
      return
      end
