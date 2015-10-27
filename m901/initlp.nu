*deck %W%  %G%
      subroutine initlp
c
c
c
      implicit real*8 (a-h,o-z)
c
      logical pagein
      integer puwk,arr,refwlk,puwkt,symorb,bmax,orbfrm
      integer ijxx(numij ),nklxx(nsym,orbfrm),ijww(numij )
      integer nklww(nsym,orbfrm),klxx(1),klww(1)
      integer kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      integer wtx(orbfrm,nsym),wty(orbfrm),wab(orbfrm),ss(norbs)
      integer itrak(2,10),itype1(4,4)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      real*8 int(nmax),c(nwks,*),s(nwks,*)
c
      common /loops/nuwk,puwkt,iuwk,juwk,itrack,ipt1,ipt2
      common /io/     itape5,itape6
      common /all/  val1,val2,val3,arr,itr1,itr2,ia,ja,itype,iseg
     *,             lvfrm1,nlwki,nlwkj,imax,imin
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf,ikl34
      common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8),ibl(8)
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
c
      data itrak/1,0,2,0,3,0,1,3,3,2,3,2,3,1,0,0,1,2,2,1/
      data itype1/1,4,8,-1,5,2,7,-1,6,9,3,10,11,12,13,17/
c
      save nwkdv2
c
      nwkdv2=nwksmx/2
      return
c
c-----------------------------------------------------------loopin
c
      entry loopin(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
c
c     write (9,*) ' entry loopin'
      puwk=puwkt
      itr1=itrak(1,itrack)
      itr2=itrak(2,itrack)
      if(itr2.eq.0) val2=0.0d+00
      if(itrack.eq.6) val2=1.0d+00
      if(itrack.eq.7) val2=1.0d+00
      if(ndvdit.le.0)  then
         itype=7
         if(iuwk.ne.juwk) return
         ja=puwkt+iuwk
         if(pagein)  then
            ja=ja-noffj
            call diagonal(int,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
         else
            nlwkt=nlwki
            ja=1
   20       if(nlwkt.gt.0) then
               nlwki=min(nlwkt,nwksmx)
               call paged
               call diagonal(int,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
               call pageout
               nlwkt=nlwkt-nlwki
               iuwk=iuwk+nlwki
            go to 20
            endif
         endif
         return
      else
         if(iuwk.eq.juwk) return
         itype=18
         if(itrack.eq.8)itype=19
         ia=puwkt+iuwk
         ja=puwkt+juwk
         if(pagein) then
            ja=ja-noffj
            ia=ia-noffi
            call external(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
         else
            ia=1
            nlwkt=nlwki
   26       if(nlwkt.le.0)return
            nlwki=min(nlwkt,nwkdv2)
            nlwkj=nlwki
            ja=nlwki+1
            call pageo
            call external(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
            call pageout
            iuwk=iuwk+nlwki
            juwk=juwk+nlwki
            nlwkt=nlwkt-nlwki
         go to 26
         endif
      endif
      return
c
c-------------------------------------------------------loopex
c
      entry loopex(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
c
c     write (9,*)' entry loopex'
      itr1=itrak(1,itrack)
      itr2=itrak(2,itrack)
      itype=itype1(ipt2,ipt1)
      if(itrack.eq.6) val2=1.0d+00
      if(itrack.eq.7) val2=1.0d+00
      if(itr2.eq.0) val2=0.0d+00
      if(ndvdit.le.0) then
         if(iuwk.ne.juwk) return
         ja=puwkt+juwk-noffj
         call diagonal(int,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
         return
      else
c        calculate non-diagonal partial externals
         if(ndvdit.le.1) then
            if(iguess.eq.0.and.itype.lt.10) return
         endif
         ia=puwkt+iuwk-noffi
         ja=puwkt+juwk-noffj
         call external(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
      endif
      return
c
c--------------------------------------------------------allext
c
c     calculate all-external contributions
      entry allext(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
c
c     write (9,*) ' entry allext'
      if(ndvdit.le.0) then
         itype=itype+3
         ja=puwkt-noffi
         ia=ja
         call diagonal(int,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      else
         itype=13+itype
         ia=puwkt-noffi
         ja=puwkt-noffj
         call fourx(int,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss
     #,                 ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww
     #,                 klxx,klww)
      endif
c
      return
      end