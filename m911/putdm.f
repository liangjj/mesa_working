*deck @(#)putdm.f	5.1  11/6/94
      subroutine putdm(acrcy,ilbl,ilbm)
      implicit real*8 (a-h,o-z)
      integer refwlk,symorb,bmax,orbfrm
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /diag/ rep,fzcore,eguess,eci,cnverg,sqcdif,czero
     *,             refwlk,mxiter,icnvg,iter,nroot
c
      dimension ilbl(26),ilbm(26)
c
cbl   intsrt=1+2*i2sec(128)
c
c..bhl
c      intsrt=2+2*i2sec(128)
c      ninr=nbf*norbs*2
c      intsrt=intsrt+i2sec(ninr)
c      nijr=nbf*(nbf+1)
c      intsrt=intsrt+2*i2sec(nijr)
c      nbfsqr= nbf * nbf *2
c      intsrt=intsrt+i2sec(nbfsqr)
c      call rsetsa(itap20,intsrt)
c      intsrt=intsrt+i2sec(112)
c      call ncdlbl(itap20,ilbl,ilbm,ngroup,nmax,nsym,acrcy,fzcore,rep)
c..bhl
c
      write(itape6,14)ilbl
  14  format(/,' label from integrals ...',26a3)
      write(itape6,15)nsym
  15  format(' the number of symmetry types is',i3)
      write(itape6,16)nmax
  16  format(' the integral group size is',i8)
      write(itape6,17)ngroup
  17  format(' the number of these groups is',i5)
      write(itape6,18)
  18  format(/)
      write(itape6,19)fzcore
      write(itape6,20)rep
      write(itape6,21)acrcy
  19  format(' frozen core energy is',f14.8)
  20  format(' nuclear repulsion  is',f14.8)
  21  format(' loop cutoff value  is',e14.4)
c..bhl
      nkl=(norbs*(norbs+1))/2
c      intsrt=intsrt+2*i2sec(nkl)
      nin=nsym*norbs
c      intsrt=intsrt+2*i2sec(nin)
c      call rsetsa(itap20,intsrt)
c..bhl
      return
c
c     call rgetsa(itap20,int1)
c     if(int1.eq.intsrt) go to 35
c     write(itape6,30)int1,intsrt
c  30 format(1h0,'error in start of integral buffer ',2i10)
c     call lnkerr(' m911: error ')
c900  write(itape6,910) ngrps,ngroup
c910  format(' in getint groups do not match ngrps=',i6,' ngroup=',i6)
c     call lnkerr(' m911: error ')
c     return
c  35 continue
c     return
      end
