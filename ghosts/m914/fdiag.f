*deck @(#)fdiag.f	2.1  10/10/91
      subroutine fdiag(ints,v,d,dvdmat,root,dvdvec,ndvdmx,thresh,nguess,
     #                 prtflg,nattim)
c
c
c     14 august 1987   pws at lanl
c     adding nroots-at-a-time option (nattim)
c
      implicit real*8 (a-h,o-z)
c
cvax  extended dummy ints,v,d,dvdmat
c
      real*8 d(nwksmx,mxvc),v(nwksmx,mxvc),ints(nmax),dvdmat(lowtri)
      real*8 root(*),dvdvec(*),thresh
      character*8 prtflg,status
c
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf,ikl34
      common /symq/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8),ibl(8)
      common /d901/ rep,fzcore,eguess,eci,refwlk,mxiter,cnverg,icnvg
     *,             iter,sqcdif,czero,nroot
      common /io/     itape5,itape6
      common /all/  val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
     *,             lvfrm1,nlwki,nlwkj,imax,imin
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common /nvectr/ nvc,mxvc
c
      logical pagein
      integer arr,refwlk,symorb,bmax,fword,fword2,orbfrm
c
      sqcdif=1.0d+00
      icnvg=-1
c   prepare for unit guess on lowest diagonal matrix element
      jj=0
      nwkmx2=intowp(nwksmx)
      nword=nwkmx2
      fword=1
ctemp
c     call iosys('write real diagonals on bliu',nwks,d,0,' ')
c
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if
c
      call bliu('initialize',status,d,thresh,nwks,mxiter,nroots,itape6,
     #         nattim,rep+fzcore,cnverg,dvdvec)
cend
      ndvdit=3
      nvc=1
      call iosys('write real tempd on bliu',nwks,d,0,' ')
      do 61 guess=1,nguess
         call iosys('read real tempd from bliu',-1,d,0,' ')
         eguess=0.0d+00
         do 50 j=1,nwks
            if(d(j,1).gt.eguess) go to 40
            refwlk=j
            eguess=d(j,1)
   40       continue
  50     continue
  60  continue
         d(refwlk,1)=0.0
         call iosys('write real tempd to bliu',nwks,d,0,' ')
c
c
         call rzero(v,nwks)
         v(refwlk,1)=1.0d+00
      if (prtflg.ne.'minimum') then
         if(guess.eq.1) then
            write(itape6,63) eguess+rep+fzcore,refwlk
         else
            write(itape6,62) eguess+rep+fzcore,refwlk
   62    format(27x,g20.12,'(',i8,')')
   63    format(5x,'reference energies:   ',g20.12,'(',i8,')')
         endif
      end if
         call rzero(d,nwks)
         call loopy
         call bliu('with vectors',0,v,d,nwks,mxiter,0,0,0,0,0,0)
   61 continue
c
c
      ndvdit=3
c
c
c      iword4=nwks2+1
c      if(irstrt.eq.0) iword6=iword4
c      eguess=eguess+rep+fzcore
c      write(itape6,70)refwlk
c  70  format(' reference configuration=',i8)
c      write(itape6,80)eguess
c  80  format(' reference energy=',g16.10)
c      if(irstrt.gt.0) go to 300
c      iter=0
c      if(iguess.ne.0) go to 200
cc   use a unit vector as a starting guess
c      eci=eguess
c      nroot=1
c      ndvdit=1
c      write(itape6,100)
c 100  format(1h ,'unit vector used as a starting guess')
cc   write out vector for davidson algorithm
c      fword=intowp(refwlk-1)+1
c      one=1.0d+00
ccpws  call wwrit(itap93,one,intowp(1),fword,junk)
c      if(nwks.le.nwksmx) v(refwlk,1)=1.0d+00
cc   prepare final vector file
ccsel  nsect=nroots*i2sec(nwks2)+i2sec(112)
ccsel  call wfile(itap12,nsect)
ccpws  call srew(itap12)
cc   do partial iteration for reference interactions
c      if(ngroup.eq.1) go to 110
ccpws  call srew(itape2)
ccvax  call rsetsa(itape2,intsrt)
ccpws  call setwa(itape2,intsrt)
ccpws  call sread(itape2,ints,nmax2)
c  110 continue
c      call loopy
ccvax  call davidson
c      call davidson(v,d,root,dvdvec,dvdmat,ndvdmx)
      return
c   read in old vector from itap12 as a starting guess
  200 continue
      if(iguess.gt.0) go to 205
      nvec=-iguess
      write(itape6,201) nvec
  201 format(1x,'reading ',i2,' approximate vector(s) from input.')
      go to 210
  205 continue
      nvec=iguess
      write(itape6,206) iguess
  206 format(1x,'reading ',i2,' old vector(s) from tape 12.')
      write(itape6,207) refwlk,eciold
  207 format(1x,'old reference=',i8,/,' old ci result=',g17.11)
  210 continue
      nroot=irooti
      if(nroot.gt.nvec) nroot=nvec
      ndvdit=nvec
cpws  fword=iround(112)+1
      fword2=1
      do 230 i=1,nvec
         if(iguess.gt.0) go to 215
         write(itape6,211) i
  211    format(1x,'vector ',i2)
         renorm=1.0d+00
         fword2=fword2+nwks2
         go to 225
  215    continue
         nword=nwkmx2
         renorm=0.0d+00
         do 220 j=1,npass
            if(j.eq.npass) nword=nwordl
            nloop=nword/intowp(1)
            do 216 k=1,nloop
  216       renorm=renorm+v(k,1)*v(k,1)
cpws        call wread(itap12,v,nword,fword,fword)
cpws0    call wwrit(itap93,v,nword,fword2,fword2)
  220    continue
  225    if(i.eq.1) go to 230
cps         call orthog(d,v,int,i,renorm)
  230 continue
      if(nvec.eq.1) return
      do 235 j=1,nwksmx
  235 d(j,1)=0.0d+00
      fword2=1
cpws  if(npass.eq.1) call wread(itap93,v,nwks2,fword2,junk)
      do 280 i=1,nvec-1
cpws     call srew(itape2)
cvax     call rsetsa(itape2,intsrt)
cpws     call setwa(itape2,intsrt)
cpws     call sread(itape2,ints,nmax2)
         call loopy
         offset=1
         nword=nwkmx2
         do 270 j=1,npass
            if(j.eq.npass) nword=nwordl
            nloop=nword/intowp(1)
            fword=offset
            ksize=(i-1)*i/2
cpws        call wread(itap94,d,nword,fword,junk)
cpws        call wread(itap93,v,nword,fword2,junk)
            do 240 k=1,nloop
               v(k,1)=v(k,1)*d(k,1)
  240       continue
            fword2=fword2+nwks2
cpws        call wread(itap94,d,nword,fword2,junk)
            do 245 k=1,nloop
               v(k,1)=v(k,1)+d(k,1)
  245       continue
            do 260 l=1,i
               z=0.0d+00
               ksize=ksize+1
cpws           call wread(itap93,d,nword,fword,junk)
               do 250 k=1,nloop
                  z=z+v(k,1)*d(k,1)
  250          continue
               dvdmat(ksize)=dvdmat(ksize)+z
               fword=fword+nwks2
  260       continue
            offset=offset+nword
  270    fword2=fword2+nword-nwks2
         iword4=iword4+nwks2
         iword3=iword3+nwks2
cpws     if(nwks.eq.nwksmx) call wread(itap93,v,nwks2,iword3,junk)
         if(irstrt.eq.0) iword6=iword4
         do 280 k=1,nwksmx
            d(k,1)=0.0d+00
  280 continue
      return
c
c   read in restart information and started mid interation
 300  continue
      irstrt=-1
      write(itape6,310) iter
 310  format(1h0,80('*')/'0calculation restarted from iteration',i4)
cvax  call davidson
cpws      call davidson(v,d,root,dvdvec,dvdmat,ndvdmx)
      delta=eci-eguess
      write(itape6,320)nroot,iter,eci,delta,sqcdif
 320  format(1h0,'root',i3,' iter',i3,' energy ',g17.11,' delta ',e9.3
     *,      ' sqcdif ',e9.3)
      return
      end
