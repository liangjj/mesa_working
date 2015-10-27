*deck @(#)start.f	5.1  11/6/94
      subroutine start(dvdmat,ints,sigma,vector)
      implicit real*8 (a-h,o-z)
c
      integer fword,fword2,bmax,refwlk,orbfrm,symorb
c
      logical pagein
c
      real*8 dvdmat(lowtri),sigma(nwksmx),vector(nwksmx)
      real*8 ints(nmax)
c
      common /io/     itape5,itape6
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /optns/  iguess,irstrt,irooti,irootf,ikl34
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
c
      data iflag1,iflag2/123456789,0/
      save iflag1,iflag2
c
      save nwkmx2
c
      inxt=1
      iter=0
      ndvdit=0
      iblck1=1
      nwkmx2=intowp(nwksmx)
      if(irstrt.ne.0) itap96=96
cpws  call sfile(itap96,'tape96')
      if(irstrt.le.0) return
c     read(itap95)ndvdit,iter,nroot
c     read(itap95) dvdmat
c     read(itap95) iblck1,inxt
c     read(itap95) iflag
      ksize=iter*(iter-1)/2
      ksize=ndvdit*(ndvdit-1)/2
      write(itape6,1)
    1 format(1x,'restart data'/' ndvdit iter nroot iblck1 inxt')
      write(itape6,2)ndvdit,iter,nroot,iblck1,inxt
    2 format(1x,i4,i6,i5,i7,i6)
      write(itape6,3)
    3 format(1x,'davidson matrix')
      k=0
      do 5 i=1,ndvdit-1
         write(itape6,4)(dvdmat(j+k),j=1,i)
    4    format(5(1x,g15.9))
         k=k+i
    5 continue
      do 10 i=ksize+1,lowtri
   10 dvdmat(i)=0.0d+00
      iword4=1+ndvdit*nwks2
      iword3=iword4-nwks2
      if(inxt.le.norbs) go to 40
      write(itape6,20)
   20 format(1x,'iteration had been completed. simply do davidson.')
      eci=eguess
cpws      call davidson
      inxt=1
      return
   40 continue
      if(iflag.ne.iflag1) go to 70
c   program died while writing 106=>104, iteration spoiled
      write(itape6,41)
   41 format(1x,'program died while writting, iteration restarted.')
      inxt=1
      iblck1=1
      do 50 i=1,nwksmx
   50 sigma(i)=0.0d+00
      fword=iword4
      nword=nwkmx2
      do 60 i=1,npass
         if(i.eq.npass) nword=nwordl
cpws     call wwrit(itap94,sigma,nword,fword,fword)
   60 continue
      return
   70 nword=nwkmx2
      fword=iword4
      fword2=1
      do 80 i=1,npass
         if(i.eq.npass) nword=nwordl
cpws     call wread(itap94,sigma,nword,fword,fword)
cpws     call wwrit(itap96,sigma,nword,fword2,fword2)
   80 continue
cpws  if(nwks.le.nwksmx) call wread(itap93,vector,nwks2,iword3
cpws #,                             junk)
      return
      entry save
c     backspace itap95
c     write(itap95) iblock,inxt
c     write(itap95) iflag1
c     backspace itap95
      fword=1
      fword2=iword4
      nword=nwkmx2
      do 110 i=1,npass
         if(i.eq.npass) nword=nwordl
cpws     call wread(itap96,sigma,nword,fword,fword)
cpws     call wwrit(itap94,sigma,nword,fword2,fword2)
  110 continue
c     write(itap95) iflag2
c     backspace itap95
      return
      end
