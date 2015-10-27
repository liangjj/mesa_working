      subroutine denom(hdenom,nbig,hpvhp,lmtop,nstate,hpvbtrn,nbfmax
     1 ,nchnl,hbbtrn,nchan,nscat,nlm,kchan,eground,obbtrn,iprint,ntot,
     2 vopt,iopt,nsch,hpptemp,nmotot,hbbtrnd,istat,itell,ihow,isym,jsym,
     $ hbbtrne,nopen)
       implicit real*8 (a-h,o-z)
c
c construct denominator (h-e) from various pieces
c
      complex*16 hdenom(nbig,nbig),hpvhp(lmtop,lmtop,nstate)
      complex*16 hpvbtrn(lmtop,nbfmax,nchnl**2)
      real*8 hbbtrn(nbfmax,nbfmax,nstate)
      real*8 hbbtrne(nbfmax,nbfmax,nstate)
      real*8 hpptemp(nbig,nbig),  vopt(nbig,nbig)
      integer nsch(nbfmax,nchan), nscat(nchnl),nlm(nchnl)
      real*8 kchan(nchnl), obbtrn(nbfmax,nbfmax,nstate)
      dimension hbbtrnd(nbfmax,nbfmax,nstate)
      character*8 istatic,istat,itell,ihow,isym,jsym
      data istatic/8hstatic  /
c
c note hpptemp is a real array which is EQUIVALENCED  to hdenom
c  It is used only during the binary read of hpp from the mesa file
c
c
c  read size of Hpp and Hopt together with number of MO's in
c  target space from mesa file
c
      if(istat.eq.istatic) go to 603
c      if(ihow.eq.itell)go to 608
c      if(istat.ne.istatic) go to 603
      write(6,194) 
 194  format(' this is not a static calculation')
      read(9) npvec,nsmall
      write(6,877)npvec,nsmall
 877  format(/' npvec,nsmall from mesa:',2i5)
c  read hpp-E in with mesa indexing
      call rdbinsqr(hpptemp,nbig,npvec,9)
c      write(6,117)
117   format(//' hpptemp matrix')
c      do 61 i=1,npvec
c 61   write(6,118) i,(hpptemp(j,i),j=1,npvec)
 118  format(1x,i3,3(f8.5,3x,f8.5,3x),/,
     &     (4x,3(f8.5,3x,f8.5,3x)))
c
c  read hopt in with mesa indexing
c
      call rdbinsqr(vopt,nbig,npvec,9)
c      write(6,127)
127   format(//' vopt matrix')
c      do 62 i=1,npvec
c 62      write(6,118) i,(vopt(j,i),j=1,npvec)
c
c  add hpp-E and hopt together, still in mesa indexing
c
      if(iopt.eq.1) then
      do 600 ip = 1,npvec
      do 600 jp=1,npvec
      vopt(ip,jp)=vopt(ip,jp)+hpptemp(ip,jp)
 600  continue
      else
c
c  if iopt syas to turn hopt off, load just hpp-E
c.
      do 602 ip=1,npvec
      do 602 jp=1,npvec
  602 vopt(ip,jp) = hpptemp(ip,jp)
      endif
c      write(6,128)
 128  format(//' vopt matrix')
c      do 63 i=1,npvec
c 63      write(6,118) i,(vopt(j,i),j=1,npvec)
      go to 605
603   continue

      nsmall=0
      ist=0
      iprev=0
      do 606 ic=1,nchan
      nsic=nscat(ic)
      jprev=0
      do 607 jc=1,ic
      en=0.0
      if(ic.eq.jc) en=0.5*kchan(ic)**2
c ***********tnr*********closed-channel fix-up
      if(ic.gt.nopen)en=-en
      ist=ist+1
      nsjc=nscat(jc)
      do 604 isc=1,nsic
      do 604 jsc=1,nsjc
      iref=isc+iprev
      jref=jsc+jprev
      hdenom(iref,jref)=hbbtrnd(isc,jsc,ist)-en*obbtrn(isc,jsc,ist)
c     $ +hbbtrne(isc,jsc,ist)
604   hdenom(jref,iref)=hdenom(iref,jref)
607   jprev=jprev+nsjc
606   iprev=iprev+nsic
      if(ihow.eq.itell)then
      do 888 i=1,iprev
      do 888 j=1,iprev
 888     hpptemp(i,j)=hdenom(i,j)
      call wrbinsqr(hpptemp,nbig,iprev,99)
      endif
      go to 608
605   continue
      do 601 i=1,nbig
      do 601 j=1,nbig
  601 hdenom(i,j) = 0.0
c
c symmetry case
      if(isym.eq.jsym)then
      do 705 i=1,npvec
      do 705 j=1,npvec
 705     hdenom(i,j)=vopt(i,j)
      iprev=npvec
      jprev=npvec
      else
c
c  move heff=hpp+hopt - E into hdenom, reindexing to conform
c  with assignment of MO's to each channel.
c  mesa indexing will have all MO's associated with every channel
c.
      iprev=0
      do 11 ic=1,nchan
      jprev=0
      nsic=nscat(ic)
      do 10 jc=1,ic
      nsjc=nscat(jc)
      ist=ic*(ic-1)/2+jc
      do 5 isc=1,nsic
      isub=isc+iprev
      do 5 jsc=1,nsjc
      jsub=jsc+jprev
c
c  This is the crucial computation of indices in mesa convention.
c  nsmall = number of MO's contributing to target configurations
c  nmotot = total number of MO's in transformation matrix
c  so (nmotot-nsmall) = number of MO's associated with every channel
c     by mesa, which does not build symmetry adapted configurations.
c  nsch(isc,ic) = index in transformation matrix of scattering orbital
c    (MO) isc in channel ic, but nsmall target MO's precede the first scattering
c    orbital in that matrix.  So we substract nsmall to get the mesa index
c    relative to the beginning of the entries for that channel.
c
      iref = nsch(isc,ic) - nsmall + (ic-1)*(nmotot-nsmall)
      jref = nsch(jsc,jc) - nsmall + (jc-1)*(nmotot-nsmall)
      hdenom(isub,jsub) = vopt(iref,jref)
5     hdenom(jsub,isub) = hdenom(isub,jsub)
10    jprev= jprev + nsjc
11    iprev=iprev + nsic
      endif
608   continue
      istart=0
      do 21 ic=1,nchan
21    istart=istart+nscat(ic)
      if(ihow.eq.itell)go to 609
      if(istart.ne.iprev.or.istart.ne.jprev) then
      write(6,101)
101   format(//' indexing error in kohnopt routine hdenom')
      stop
      endif
 609  continue
c
c free free part
c
      iprev=0
c.....tnr...9-94.......
c      do 31 ic=1,nchan
      do 31 ic=1,nopen
      nlmic=nlm(ic)
      jprev=0
      do 30 jc=1,ic
      nlmjc=nlm(jc)
      ist=ic*(ic-1)/2 + jc
      do 29 ilm=1,nlmic
      isub=ilm+iprev+istart
      do 29 jlm=1,nlmjc
      jsub=jlm+jprev+istart
      hdenom(isub,jsub)=hpvhp(ilm,jlm,ist)
29    hdenom(jsub,isub)=hdenom(isub,jsub)
30    jprev=jprev+nlmjc
31    iprev=iprev+nlmic
c
c bound-free part
c
      iprev=istart
c.....tnr...9-94.......
c      do 41 ic=1,nchan
      do 41 ic=1,nopen
      nlmic=nlm(ic)
      jprev=0
      do 40 jc=1,nchan
      nsjc=nscat(jc)
      icc=nchan*(ic-1) + jc
      do 39 ilm=1,nlmic
      isub=ilm+iprev
      do 39 jsc=1,nsjc
      jsub=jsc+jprev
      hdenom(isub,jsub) = hpvbtrn(ilm,jsc,icc)
39    hdenom(jsub,isub) =hdenom(isub,jsub)
40    jprev=jprev+nsjc
41    iprev=iprev+nlmic
      ntot=iprev
      if(ihow.eq.itell)then
      iw=istart+1
      if(isym.ne.jsym)then
      iw=1
      endif
      write(6,*)'iw,ntot',iw,ntot
      write(7)((hdenom(i,j),i=1,ntot),j=iw,ntot)
      endif
c
      if(iprint.ne.0) then
      write(6,107)
107   format(//' denominator matrix of (h-e)')
      do 60 i=1,ntot
60    write(6,108) i,(hdenom(j,i),j=1,ntot)
108   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      endif
      return
      end
