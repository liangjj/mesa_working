*deck @(#)initlp.f	2.1  10/10/91
      subroutine initlp
c
c
c
      implicit real*8 (a-h,o-z)
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
      integer puwk,arr,refwlk,puwkt,symorb,bmax,orbfrm
      logical pagein
      dimension itrak(2,10),itype1(4,4)
      data itrak/1,0,2,0,3,0,1,3,3,2,3,2,3,1,0,0,1,2,2,1/
      data itype1/1,4,8,-1,5,2,7,-1,6,9,3,10,11,12,13,17/
      nwkdv2=nwksmx/2
      return
c
c-----------------------------------------------------------loopin
c
      entry loopin
c
c      write (itape6,*) ' entry loopin'
c
c
      puwk=puwkt
      itr1=itrak(1,itrack)
      itr2=itrak(2,itrack)
      if(itr2.eq.0) val2=0.0d+00
      if(itrack.eq.6) val2=1.0d+00
      if(itrack.eq.7) val2=1.0d+00
      if(ndvdit.gt.0) go to 25
      itype=7
      if(iuwk.ne.juwk) return
      ja=puwkt+iuwk
c..bhl
c      if(pagein) go to 21
c      nlwkt=nlwki
c      ja=1
c   20 if(nlwkt.le.0)return
c      nlwki=min(nlwkt,nwksmx)
c      call paged
c      call diagonal
c      call pageout
c      nlwkt=nlwkt-nlwki
c      iuwk=iuwk+nlwki
c      go to 20
c..bhl
   21 ja=ja-noffj
c      write(itape6,*)'  calling diagonal '
c      write(itape6,90012)noffi,noffj,ia,ja,iuwk,juwk,nlwkt,itype,
c     $ itrack
90012  format(10(2x,i6))
      call diagonal
      return
   25 if(iuwk.eq.juwk) return
      itype=18
      if(itrack.eq.8)itype=19
      ia=puwkt+iuwk
      ja=puwkt+juwk
c..bhl
c      if(pagein) go to 27
c      ia=1
c      nlwkt=nlwki
c   26 if(nlwkt.le.0)return
c      nlwki=min(nlwkt,nwkdv2)
c      nlwkj=nlwki
c      ja=nlwki+1
c      call pageo
c      call external
c      call pageout
c      iuwk=iuwk+nlwki
c      juwk=juwk+nlwki
c      nlwkt=nlwkt-nlwki
c      go to 26
c..bhl
   27 ja=ja-noffj
      ia=ia-noffi
c      write(itape6,*)'  calling diagonal '
c      write(itape6,90012)noffi,noffj,ia,ja,iuwk,juwk,nlwkt,itype,
c     $ itrack
c   90012  format(10(2x,i6))
      call external
      return
c
c-------------------------------------------------------loopex
c
      entry loopex
c
c     write (9,*)' entry loopex'
      itr1=itrak(1,itrack)
      itr2=itrak(2,itrack)
      itype=itype1(ipt2,ipt1)
      if(itrack.eq.6) val2=1.0d+00
      if(itrack.eq.7) val2=1.0d+00
      if(itr2.eq.0) val2=0.0d+00
      if(ndvdit.gt.0) go to 40
      if(iuwk.ne.juwk) return
      ja=puwkt+juwk-noffj
      call diagonal
      return
c   calculate non-diagonal partial externals
   40 continue
c..bhl
c 10/20/89      if(ndvdit.gt.1) go to 45
c..bhl
      if(ndvdit.gt.0) go to 45
      if(iguess.eq.0.and.itype.lt.10) return
   45 continue
      ia=puwkt+iuwk-noffi
      ja=puwkt+juwk-noffj
      call external
      return
c
c--------------------------------------------------------allext
c
c   calculate all-external contributions
      entry allext
c
c     write (9,*) ' entry allext'
      if(ndvdit.gt.0) go to 80
      itype=itype+3
      ja=puwkt-noffi
      ia=ja
      call diagonal
      return
   80 itype=13+itype
      ia=puwkt-noffi
      ja=puwkt-noffj
      call fourx
      return
      end
