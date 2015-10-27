*deck @(#)initlp.f	5.1  11/6/94
      subroutine initlp(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $             ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,
     $             klxx,klww)
c
c
c
      implicit real*8 (a-h,o-z)
      integer puwk,arr,refwlk,puwkt,symorb,bmax,orbfrm
      integer ijxx(numij), nklxx(nsym,orbfrm), ijww(numij)
      integer nklww(nsym,orbfrm), klxx(1), klww(1)
      integer kadd(symorb),ladd(symorb),ijadd(numij),wtw(orbfrm,nsym)
      integer wtx(orbfrm,nsym), wty(orbfrm), wab(orbfrm), ss(norbs)
      real*8 ci(1),cj(1),si(1),sj(1),h1(1),h2(1)
      real*8 int(nmax), int1(norbs*norbs),c(nwks,*), s(nwks,*)
      logical pagein
      integer itrak(2,10),itype1(4,4)
      common /loops/nuwk,puwkt,iuwk,juwk,itrack,ipt1,ipt2
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /all/  val1,val2,val3,arr,itr1,itr2,ia,ja,itype,iseg
     *,             lvfrm1,nlwki,nlwkj,imax,imin
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf
      common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
      common /diag/ rep,fzcore,eguess,eci,cnverg,sqcdif,czero
     *,             refwlk,mxiter,icnvg,iter,nroot
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      data itrak/1,0,2,0,3,0,1,3,3,2,3,2,3,1,0,0,1,2,2,1/
      data itype1/1,4,8,-1,5,2,7,-1,6,9,3,10,11,12,13,17/
      save itrak,itype1
c
      save nwkdv2
c
      nwkdv2=nwksmx/2
      return
c
c

c..rlm      entry loopin
      entry loopin(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $             ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,
     $             klxx,klww)
c
      puwk=puwkt
      itr1=itrak(1,itrack)
      itr2=itrak(2,itrack)
      if(itr2.eq.0) val2=0.0d+00
      if(itrack.eq.6) val2=1.0d+00
      if(itrack.eq.7) val2=1.0d+00
c     if(ndvdit.gt.0) go to 25
      itype=7
c     if(iuwk.ne.juwk) return
      if (iuwk.ne.juwk) go to 25
      ja=puwkt+iuwk
      if(pagein) go to 21
      nlwkt=nlwki
      ja=1
   20 if(nlwkt.le.0)return
      nlwki=min(nlwkt,nwksmx)
      call paged
      call diagonal(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      call pageout
      nlwkt=nlwkt-nlwki
      iuwk=iuwk+nlwki
      go to 20
   21 ja=ja-noffj
      call diagonal(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
      return
   25 continue
      itype=18
      if(itrack.eq.8)itype=19
      ia=puwkt+iuwk
      ja=puwkt+juwk
      if(pagein) go to 27
      ia=1
      nlwkt=nlwki
   26 if(nlwkt.le.0)return
      nlwki=min(nlwkt,nwkdv2)
      nlwkj=nlwki
      ja=nlwki+1
      call pageo
      call external(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $              ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      call pageout
      iuwk=iuwk+nlwki
      juwk=juwk+nlwki
      nlwkt=nlwkt-nlwki
      go to 26
   27 ja=ja-noffj
      ia=ia-noffi
      call external(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $              ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      return
c
c
c---------------------------------------------------------loopex
c
      entry loopex(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $             ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,
     $             klxx,klww)
c
c
      itr1=itrak(1,itrack)
      itr2=itrak(2,itrack)
      itype=itype1(ipt2,ipt1)
      if(itrack.eq.6) val2=1.0d+00
      if(itrack.eq.7) val2=1.0d+00
      if(itr2.eq.0) val2=0.0d+00
c     if(ndvdit.gt.0) go to 40
c     if(iuwk.ne.juwk) return
      if (iuwk.ne.juwk) go to 40
      ja=puwkt+juwk-noffj
c..rlm     call diagonal
      call diagonal(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
c     return
c   calculate non-diagonal partial externals
   40 continue
c     if(ndvdit.gt.1) go to 45
c     if(iguess.eq.0.and.itype.lt.10) return
c  45 continue
      ia=puwkt+iuwk-noffi
      ja=puwkt+juwk-noffj
c..rlm     call extern
      call external(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $              ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,
     $              klxx,klww)
      return
c
c-------------------------------------------------------allext
c
c   calculate all-external contributions
c..rlm     entry allext
      entry allext(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $             ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww,
     $             klxx,klww)
c
c     if(ndvdit.gt.0) go to 80
      itype=itype+3
      ja=puwkt-noffi
      ia=ja
      call diagonal(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss)
c..rlm     call diagon
c     return
c
c     ----- following line replaced for dm -----
c..rlm   seems unreachable.
c
c  80 itype=13+itype
   80 continue
      itype=itype+10
c
      ia=puwkt-noffi
      ja=puwkt-noffj
      call fourx(int,int1,c,s,ijadd,kadd,ladd,wtw,wtx,wty,wab,ss,
     $             ci,si,cj,sj,h1,h2,ijxx,nklxx,ijww,nklww)
      return
      end
