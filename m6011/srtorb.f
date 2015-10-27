*deck  @(#)srtorb.f	1.1 9/7/91
c***begin prologue     m6011
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6011, link 6011, orbital sort
c***author             schneider, barry (lanl)
c***source             m6011
c***purpose            sort initial numerical orbital file
c***                   into a subset kept for m6005
c***description        
c***                   
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6011
      program srtorb
      implicit integer (a-z)
      common a(1)
      common/ memory / ioff
      common /io/ inp, iout
      real*8 z, norm
      character *8 cpass, chrkey
      character *3 ans
      character *800 card
      character *128 orbfl, grdfl
      dimension z(1), lstin(200), list(200), norm(200)
      equivalence (z,a)
      call drum
      write (iout,100)
      write (iout,10)
      write (iout,100)
      call posinp('$srtorb',cpass)
      call cardin(card)
      call iosys ('read character "orbital filename" from rwf',
     1             -1,0,0,orbfl)
      call iosys ('read character "grid filename" from rwf',-1,0,0,
     1             grdfl)
      call iosys ('open orbs as old',0,0,0,orbfl)
      call iosys ('read integer "no. regions" from orbs',1,nreg,0,' ')
      call iosys ('read integer "no. grid pts" from orbs',1,npnts,0,' ')
      call iosys ('read integer "point buffer" from orbs',1,pntbuf,
     1            0,' ')
      call iosys ('read integer "no. cont" from orbs',1,ncon,0,' ')
      nkept=intkey(card,'no-functions-kept',ncon,' ')
      call intarr(card,'function-list',lstin,nkept,' ')
      call izero(list,200)
      if (nkept.eq.ncon) then
          do 50 i=1,ncon
             list(i)=i
   50     continue
      else
          do 60 i=1,nkept
             ii=lstin(i)
             list(ii)=i
   60     continue
      endif
      call iosys ('read integer "final pts" from orbs',1,nolst,0,' ')
      call iosys ('write integer "no. kept" to orbs',1,nkept,0,' ')
      call iosys ('write integer "function list" to orbs',200,list,
     1            0,' ')
      call iosys ('open grid as old',0,0,0,grdfl)
      memmax=intkey(card,'maximum-memory',1000000,' ')
      size=(nreg+1)*pntbuf
      call getscm(0,a,canget,'how much core',0)
      canget=canget-4*pntbuf
      if (memmax.lt.0) memmax=canget
      noarr=memmax/size
      noarr=min(noarr,ncon+nkept)
      numin=noarr/2
      numin=min(numin,ncon)
      numout=noarr-numin 
      if (numin.lt.1) then
          call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
          call iosys ('close orbs',0,0,0,' ')
          call lnkerr ('cannot get space for input array')
      endif
      words=wpadti(noarr*pntbuf+4*pntbuf)
      call iosys ('read integer maxsiz from rwf',1,maxsiz,0,' ')
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,a,ngot,'m6011',0)
      write (iout,20) noarr
      sizin=npnts*ncon
      sizout=nkept*npnts
      call iosys ('create real "sorted orbs" on orbs',sizout,0,0,' ')
      grid=ioff
      fin=grid+4*pntbuf
      fout=fin+numin*pntbuf
      nwdin=0
      nwdout=0
      nwrin=0
      nwrout=0
      call rzero(norm,200)
      do 200 ireg=1,nreg
         noptrg=pntbuf
         if (ireg.eq.nreg) then
             noptrg=nolst
         endif
         call srt (z(fin),z(fout),z(grid),ncon,numin,numout,nkept,
     1             lstin,noptrg,nwrin,nwdin,nwrout,nwdout,norm,ireg)
  200 continue
      if (nwdin.ne.npnts*ncon) then
          call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
          call iosys ('close orbs',0,0,0,' ')
          call iosys ('rewind all on grid read-and-write',0,0,0,' ')
          call iosys ('close grid',0,0,0,' ')
          call lnkerr('error in input word count to orbs')
      endif
      if (nwdout.ne.npnts*nkept) then
          call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
          call iosys ('close orbs',0,0,0,' ')
          call iosys ('rewind all on grid read-and-write',0,0,0,' ')
          call iosys ('close grid',0,0,0,' ')
          call lnkerr('error in output word count to orbs')
      endif
      write (iout,30) nwrin, nwrout
      write (iout,40) nwdin, nwdout
      write (iout,300) (norm(i),i=1,nkept)
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxsiz,0,' ')
      call chainx(0)
      stop
   10 format(/,20x,'link:m6011:numerical orbital sort')
   20 format(/,5x,'store',1x,i4,1x,'arrays in memory per pass')
   30 format (/,5x,'no. input reads',1x,i5,1x,'no. output writes',1x,i5)
   40 format (/,5x,'no. input words',1x,i8,1x,'no. output words',1x,i8)
  100 format(/,20x,'***********************************************')
  300 format(/,5x,'normalization integrals',(/,5x,5e15.8)) 
      end
