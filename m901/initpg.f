*deck @(#)initpg.f	5.1  11/6/94
      subroutine initpg(vector,sigma)
c
c
c
      implicit real*8 (a-h,o-z)
c
cvax  extended dummy vector,sigma
c
      dimension vector(nwksmx),sigma(nwksmx)
      integer puwkt,refwlk,bmax,fword,orbfrm
      logical pagein,ieqj
c
      integer nr3,nr6,nw6
      common /statpg/ ir3(100),nr3(100),ir6(100),nr6(100),iw6(100)
     #,               nw6(100),ipgs(100),ipgd(100),ipgo(100),ipgij(100)
     #,               ipgout(100)
      common /level/  ilev
c
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /io/     itape5,itape6
      common /optns/  iguess,irstrt,irooti,irootf,ikl34
      common /loops/nuwk,puwkt,iuwk,juwk,itrak,ipt1,ipt2
      common /all/  acf,d,ccf,ladt,itr1,itr2,ia,ja,itype,isegt,lvfrm1
     *,             nlwki,nlwkj,imax,imin
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /rstrt/iblck1,inxt,lowtri,ndvdit,iblock
      common /pagdat/ieqj,iuwksv,juwksv,nlwki2,nlwkj2,jwalk,iwalk,ii,jj
      iword3=1
      iword4=1
      iword6=1
      ieqj=.true.
      iwalk=0
      ii=1
      jj=1
      nlwki2=intowp(nwks2)
      nwkmx2=nwksmx
      return
c
c-----------------------------------------------------------------pages
c
      entry pages
c
c     write (9,*) ' entry pages '
      ipgs(ilev)=ipgs(ilev)+1
      if(iuwk.ne.juwk) go to 210
c
c-----------------------------------------------------------------paged
c
      entry paged
c
c     write (9,*) ' entry paged '
      ipgd(ilev)=ipgd(ilev)+1
  200 if(nlwki.gt.nwksmx) return
      iuwksv=iuwk
      juwksv=juwk
      ieqj=.true.
      pagein=.true.
      nlwki2=intowp(nlwki)
      nlwkj2=nlwki2
      noffi=puwkt+iuwk-1
      noffj=noffi
      iwalk=intowp(noffi)
      ii=1
      jj=1
      fword=iwalk+iword6
cpws  call wread(itap96,sigma,nlwki2,fword,junk)
      call lnkerr('paged problem')
      ir6(ilev)=ir6(ilev)+1
      nr6(ilev)=nr6(ilev)+nlwki2
      if(ndvdit.eq.0) return
      fword=iwalk+iword3
      call lnkerr('paged problem')
cpws  call wread(itap93,vector,nlwki2,fword,junk)
      ir3(ilev)=ir3(ilev)+1
      nr3(ilev)=nr3(ilev)+nlwki2
      return
  210 continue
      if((nlwki+nlwkj).gt.nwksmx) return
c
c-----------------------------------------------------------------pageo
c
      entry pageo
c
c     write (9,*) ' entry pageo '
      ipgd(ilev)=ipgd(ilev)+1
      iuwksv=iuwk
      juwksv=juwk
      pagein=.true.
      ieqj=.false.
      nlwki2=intowp(nlwki)
      nlwkj2=intowp(nlwkj)
      noffi=puwkt+iuwk-1
      noffj=puwkt+juwk-1
      iwalk=intowp(noffi)
      jwalk=intowp(noffj)
      noffj=noffj-nlwki
      fword=iwalk+iword6
      ii=1
      call lnkerr('pageo error')
cpws  call wread(itap96,sigma,nlwki2,fword,junk)
      ir6(ilev)=ir6(ilev)+1
      nr6(ilev)=nr6(ilev)+nlwki2
      fword=jwalk+iword6
      jj=nlwki+1
cpws  call wread(itap96,sigma(jj),nlwkj2,fword,junk)
      ir6(ilev)=ir6(ilev)+1
      nr6(ilev)=nr6(ilev)+nlwkj2
      fword=iwalk+iword3
cpws  call wread(itap93,vector,nlwki2,fword,junk)
      ir3(ilev)=ir3(ilev)+1
      nr3(ilev)=nr3(ilev)+nlwki2
      fword=jwalk+iword3
cpws  call wread(itap93,vector(jj),nlwkj2,fword,junk)
      ir3(ilev)=ir3(ilev)+1
      nr3(ilev)=nr3(ilev)+nlwkj2
      return
c
c----------------------------------------------------------------pageij
c
      entry pageij
c
c   check that something has changed
c
c     write (9,*) ' entry pageij '
      ipgij(ilev)=ipgij(ilev)+1
      if(iuwksv.eq.iuwk.and.juwksv.eq.juwk) return
      if(ieqj) go to 400
      if(iuwksv.eq.iuwk) go to 310
      if(juwksv.eq.juwk) go to 330
      pagein=.false.
      fword=iwalk+iword6
      call lnkerr('pageij error')
cpws  call wwrit(itap96,sigma(ii),nlwki2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwki2
      fword=jwalk+iword6
cpws  call wwrit(itap96,sigma(jj),nlwkj2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwkj2
      if(iuwk.eq.juwk) go to 200
      go to 210
  310 continue
      fword=jwalk+iword6
      call lnkerr('pageij error')
cpws  call wwrit(itap96,sigma(jj),nlwkj2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwkj2
      if((nlwki+nlwkj).le.nwksmx) go to 320
      fword=iwalk+iword6
cpws  call wwrit(itap96,sigma(ii),nlwki2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwki2
      pagein=.false.
      return
  320 continue
      jj=nlwki+1
      juwksv=juwk
      nlwkj2=intowp(nlwkj)
      noffj=puwkt+juwk-1
      jwalk=intowp(noffj)
      noffj=noffj-nlwki
      fword=jwalk+iword6
      call lnkerr('pageij error')
cpws  call wread(itap96,sigma(jj),nlwkj2,fword,junk)
      ir6(ilev)=ir6(ilev)+1
      nr6(ilev)=nr6(ilev)+nlwkj2
      fword=jwalk+iword3
cpws  call wread(itap93,vector(jj),nlwkj2,fword,junk)
      ir3(ilev)=ir3(ilev)+1
      nr3(ilev)=nr3(ilev)+nlwkj2
      return
  330 continue
      fword=iwalk+iword6
      call lnkerr('pageij error')
cpws  call wwrit(itap96,sigma(ii),nlwki2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwki2
      if(iuwk.ne.juwk) go to 340
      ieqj=.true.
      ii=jj
      iwalk=jwalk
      iuwksv=iuwk
      nlwki2=nlwkj2
      return
  340 continue
      if(nlwki.lt.jj) go to 350
      fword=jwalk+iword6
      call lnkerr('pageij error')
cpws  call wwrit(itap96,sigma(jj),nlwkj2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwkj2
      pagein=.false.
      go to 210
  350 iuwksv=iuwk
      noffi=puwkt+iuwk-1
      iwalk=intowp(noffi)
      nlwki2=intowp(nlwki)
      fword=iwalk+iword6
      ii=1
      call lnkerr('pageij error')
cpws  call wread(itap96,sigma(ii),nlwki2,fword,junk)
      ir6(ilev)=ir6(ilev)+1
      nr6(ilev)=nr6(ilev)+nlwki2
      fword=iwalk+iword3
cpws  call wread(itap93,vector(ii),nlwki2,fword,junk)
      ir3(ilev)=ir3(ilev)+1
      nr3(ilev)=nr3(ilev)+nlwki2
      return
  400 continue
      if(iuwk.ne.iuwksv) go to 420
      ieqj=.false.
      if(ii.eq.1) go to 405
      if(nlwkj.ge.ii) go to 420
      jj=1
      goto 410
  405 continue
      if((nlwki+nlwkj).gt.nwksmx) go to 420
      jj=ii+nlwki
  410 juwksv=juwk
      noffj=puwkt+juwk-1
      jwalk=intowp(noffj)
      noffj=noffj-nlwki
      nlwkj2=intowp(nlwkj)
      fword=jwalk+iword6
      call lnkerr('pageij error')
cpws  call wread(itap96,sigma(jj),nlwkj2,fword,junk)
      ir6(ilev)=ir6(ilev)+1
      nr6(ilev)=nr6(ilev)+nlwkj2
      fword=jwalk+iword3
cpws  call wread(itap93,vector(jj),nlwkj2,fword,junk)
      ir3(ilev)=ir3(ilev)+1
      nr3(ilev)=nr3(ilev)+nlwkj2
      return
  420 continue
      pagein=.false.
      fword=iwalk+iword6
      call lnkerr('pageij error')
cpws  call wwrit(itap96,sigma(ii),nlwki2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwki2
      if(iuwk.eq.juwk) go to 200
      go to 210
c
c----------------------------------------------------------------pageout
c
      entry pageout
c
c     write (9,*) ' entry pageout '
      ipgout(ilev)=ipgout(ilev)+1
      pagein=.false.
      fword=iwalk+iword6
      call lnkerr('pageout error')
cpws  call wwrit(itap96,sigma(ii) ,nlwki2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwki2
      if(ieqj) return
      fword=jwalk+iword6
cpws  call wwrit(itap96,sigma(jj),nlwkj2,fword,junk)
      iw6(ilev)=iw6(ilev)+1
      nw6(ilev)=nw6(ilev)+nlwkj2
      return
      end
