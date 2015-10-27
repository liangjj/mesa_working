      program phot
      implicit real*8(a-h,o-z)
      parameter (nbig=700, msmall=100,nchnl=4)
c      parameter (mtop=2,ltop=2)
c      parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))
      parameter (lmtop=15)
      parameter (lm2=2*lmtop)
      parameter (nbfmax=500)
      parameter (maxene=200)
c
      integer nscat(nchnl),nhole(nchnl)
      integer nsch(nbfmax,nchnl)
      real*8 echan(nchnl), energy(maxene)
      real*8 kchan(nchnl),crz(maxene,nchnl)
      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
c
      dimension denp(nbfmax,nchnl),rdenq(nbfmax*nbfmax)
      dimension cdenq(2,nbfmax,nbfmax)
      complex*16 denq(nbfmax,nbfmax),tmat(3,lmtop)
        dimension hbx(nbfmax,nbfmax),aip(nchnl)
        dimension hby(nbfmax,nbfmax),ome(maxene,nchnl)
        dimension hbz(nbfmax,nbfmax),cry(maxene,nchnl)
        dimension cs(nchnl,3),phase(nchnl),crx(maxene,nchnl)
        complex*16 htop(nbig,msmall),c(nbig),csum
        complex*16 hpvbx(lmtop,nbfmax,nchnl),amp(nchnl)
        complex*16 hpvby(lmtop,nbfmax,nchnl)
        complex*16 hpvbz(lmtop,nbfmax,nchnl)
        complex*16 ampx(nchnl),ampy(nchnl),ampz(nchnl),temp
        complex*16 ai,aa,bb,cp,cgamma
        real*8 carg
                integer iww(3)
        character*6 cart(3),ylm(3)
        data cart/'   X  ','   Y  ','   Z  '/
        data  ylm/' 1 -1 ',' 1  1 ',' 1  0 '/
        data iww/1,2,3/
c
       equivalence(cdenq,denq)
c
       data  pi/3.14159265358/
c
      ai=(0.,1.) 
c
      open(5,file='inphot')
      open(6,file='outphot')
      open(7,file='pltphot')
      open(8,file='dipbf',form='unformatted')
      open(9,file='bcoef',form='unformatted')
      open(10,file='denphot',form='unformatted')
      open(11,file='punphot',form='unformatted')
c

c
       read(8) nener,nchan,(nlm(ic),ic=1,nchan)
     $ ,(nscat(ic),ic=1,nchan)
      read(8) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
      read(8)((nsch(j,ic),j=1,nscat(ic)),ic=1,nchan)
      read(8) eground
      read(8) nmotot
      read(8)((hbx(ig,jg),ig=1,nmotot),jg=1,nmotot)
      read(8)((hby(ig,jg),ig=1,nmotot),jg=1,nmotot)
      read(8)((hbz(ig,jg),ig=1,nmotot),jg=1,nmotot)
      read(9)nchan,nall,nfree
c
      read(10) nsmall,nbf,nen,nchn
c
      write(6,77001) nsmall,nbf,nen,nchn
77001 format(//,' MESA Input ',/,
     $' nsmall = ',i5,/,
     $' nbf    = ',i5,/,
     $' nen    = ',i5,/,
     $' nchn   = ',i5,//)
c
      if(nchn.ne.nchan .or. nener.ne.nen 
     $ .or. nbf.ne.nmotot) then
       write(6,*)' Consistency Error in Tape Input'
       write(6,*)' nchan nener nmotot ',nchan,nener,nmotot
       write(6,*)' nchn nen nbf ',nchn,nen,nbf
       stop
      end if
c
c read in hole states
c
c     read in IP in eV
c
      read(5,*) iprt,nww
      if(nww.ne.0) then 
         read(5,*) (iww(i),i=1,nww) 
      else
         nww=3
      end if
      read(5,*)(aip(i),i=1,nchan)
      write(6,*)' IP is ',(aip(i),i=1,nchan)
      write(6,109)(cart(iww(i)),i=1,nww)
      write(6,110)(ylm(iww(i)),i=1,nww)
 109  format(/,' The Following Transition Amplitudes Will Be Stored ',/
     $       ' Cartesian ',3(2x,a6))
 110  format(' Y(lm)     ',3(2x,a6)) 
c
c      read(5,*)(nhole(i),i=1,nchan)
c      read(5,*)(phase(i),i=1,nchan)
c      write(6,101) (nhole(i),i=1,nchan)
c      write(6,301) (phase(i),i=1,nchan)
cc
c 101  format(' holes in orbital number ',10i3//)
c 301  format(' state phases: ',10f5.1//)
c
c  Read in P-Space Density
c
      do 101 ii=1,nchan
c
         read(10)(denp(iq,ii),iq=1,nsmall)
c     
      if(iprt.gt.0) then 
      write(6,*)' Reading P-Density ',ii
       if(iprt.gt.1) then
          write(6,1044)(denp(iq,ii),iq=1,nsmall)
 1044     format(5(2x,f12.8))
       end if
      end if
c
 101  continue
c
      icx=0
c
c open a loop on energies
c
      do 1000 iene=1,nener
c
      if(iprt.gt.0) write(6,*)' Energy ',iene
c
      tx=0.
      ty=0.
      tz=0.
c
      read(9) (kchan(ic),ic=1,nchan)
      read(9)((htop(i,j),i=1,nall),j=1,nfree)
      read(8) (kchan(ic),ic=1,nchan)
      do 222 ic=1,nchan
      nlmic=nlm(ic)
      ngjc=nmotot
      read(8)((hpvbx(il,ig,ic),il=1,nlmic),ig=1,ngjc)
      read(8)((hpvby(il,ig,ic),il=1,nlmic),ig=1,ngjc)
      read(8)((hpvbz(il,ig,ic),il=1,nlmic),ig=1,ngjc)
c     
      if(iprt.gt.2) then
         write(6,2221) ic
 2221    format(//,' hpvbz:  channel  ',i5)
        do 2222 iu=1,nsmall,4
           iue=min(iu+3,nsmall)
           write(6,2225)(io,io=iu,iue)
 2225      format(/,6x,4(8x,i5,9x))
 2223      format(i4,2x,4(2x,2(2x,f8.4)))
        do 2224 il=1,nlmic
           write(6,2223)il,(hpvbz(il,ig,ic),ig=iu,iue)
 2224   continue
 2222 continue
      end if
c     
 222  continue
c
c loop over incident channels
c
      jstart=1
      do 504 ic=1,nchan
c
      if(iprt.gt.0) write(6,*)' ion channel ',ic
c
      omega=aip(ic)/27.211+kchan(ic)**2/2.
c
      xnorm=sqrt(8.*pi*omega/3./137.*100.*.529*.529)
c
      nlmic=nlm(ic)
      jstop=jstart+nlmic-1
      ngjc=nmotot
c
c
c     loop over physical l-channels
c
      ilm=0
      totx=0.
      toty=0.
      totz=0.
      do 111 iph=jstart,jstop
         ilm=ilm+1
      if(iprt.gt.0) then
      write(6,*)' physical l,m:  ',lch(ilm,ic),mch(ilm,ic)
      end if 
c     
c     load solution vector
c
      do 1 i=1,nall
      c(i)=htop(i,iph)
 1    continue      
      if(iprt.gt.2) then
      write(6,*)' kohn coefficient vector:',ic,ilm
      write(6,112)(c(i),i=1,nall)
      end if
 112   format(2e12.4)
      ampx(ic)=0.
      ampy(ic)=0.
      ampz(ic)=0.
      ii=0
      do 2 jc=1,nchan
      njc=nscat(jc)
cc
c      nh=nhole(jc)
c      ph=phase(jc)
cc
      do 3 kc=1,njc
      ii=ii+1
      do 33 nh=1,nsmall
      ampx(ic)=ampx(ic)-c(ii)*hbx(nh,nsch(kc,jc))*denp(nh,jc)
      ampy(ic)=ampy(ic)-c(ii)*hby(nh,nsch(kc,jc))*denp(nh,jc)
      ampz(ic)=ampz(ic)-c(ii)*hbz(nh,nsch(kc,jc))*denp(nh,jc)
 33   continue
c
c      write(6,986)ii,c(ii),hbz(nh,nsch(kc,jc)),ampz(ic)
 986  format(i5,6e12.4)
c
 3     continue
 2      continue
c
        if(iprt.gt.0) then
        write(6,*)' bound-bound part'
        write(6,102)ampx(ic),ampy(ic),ampz(ic)
        end if
c
        do 4 jc=1,nchan
           njc=nlm(jc)
c           nh=nhole(jc)
c           ph=phase(jc)

           do 5 kc=1,njc
              ii=ii+1
           do 55 nh=1,nsmall
              ampx(ic)=ampx(ic)-c(ii)*hpvbx(kc,nh,jc)*denp(nh,jc)
              ampy(ic)=ampy(ic)-c(ii)*hpvby(kc,nh,jc)*denp(nh,jc)
              ampz(ic)=ampz(ic)-c(ii)*hpvbz(kc,nh,jc)*denp(nh,jc)
 55           continue
c         write(6,986)ii,c(ii),hpvbz(kc,nh,jc),ampz(ic)
 5            continue
 4            continue
c
c              nh=nhole(ic)
c              ph=phase(ic)
c
c
           do 66 nh=1,nsmall
              ampx(ic)=ampx(ic)+imag(hpvbx(ilm,nh,ic))*denp(nh,ic)
              ampy(ic)=ampy(ic)+imag(hpvby(ilm,nh,ic))*denp(nh,ic)
              ampz(ic)=ampz(ic)+imag(hpvbz(ilm,nh,ic))*denp(nh,ic)
 66        continue
c   
              if(iprt.gt.0) then 
              write(6,*)' full scattered wave part'
              write(6,102)ampx(ic),ampy(ic),ampz(ic)
              end if
c
cc
c Read the Q-Space Density
cc
              if(iprt.gt.0)   write(6,*)' Reading Q-Density '
c
      do 67 k=1,2
c      icx=icx+1
c      write(6,*)' Reading Q-Density ',icx 
      read(10)(rdenq(ii),ii=1,nsmall*nsmall) 
      ix=0
       do 68 ii=1,nsmall
        do 69 jj=1,nsmall
         ix=ix+1
         cdenq(k,ii,jj)=rdenq(ix)
 69     continue 
 68    continue
 67   continue
ccc
c Q-Space Contributions
cc
      if(iprt.gt.0) write(6,*)' Q-Space Part '
c
      do 71 ii=1,nsmall
      do 72 jj=1,nsmall
      ampx(ic)=ampx(ic)-hbx(ii,jj)*denq(ii,jj)
      ampy(ic)=ampy(ic)-hby(ii,jj)*denq(ii,jj)
      ampz(ic)=ampz(ic)-hbz(ii,jj)*denq(ii,jj)
 72   continue
 71   continue
c
c store transition amplitudes for Differenital Cross Sections
c 
      aa=float(lch(ilm,ic))+1.-ai/kchan(ic)
      bb=cgamma(aa)
      sl=carg(bb)-float(lch(ilm,ic))*pi/2.
      cp=exp(+ai*sl)
c
c      tmat(1,ilm)=ampx(ic)*xnorm*ai**lch(ilm,ic)*cp
c      tmat(2,ilm)=ampy(ic)*xnorm*ai**lch(ilm,ic)*cp
c      tmat(3,ilm)=ampz(ic)*xnorm*ai**lch(ilm,ic)*cp
c
      tmat(1,ilm)=ampx(ic)*xnorm*cp
      tmat(2,ilm)=ampy(ic)*xnorm*cp
      tmat(3,ilm)=ampz(ic)*xnorm*cp
c
      jstart=jstop+1
      cs(ic,1)=abs(ampx(ic))**2
      cs(ic,2)=abs(ampy(ic))**2
      cs(ic,3)=abs(ampz(ic))**2
c
      tx=tx+abs(tmat(1,ilm))**2
      ty=ty+abs(tmat(2,ilm))**2
      tz=tz+abs(tmat(3,ilm))**2
c
      if(iprt.gt.0) then
      write(6,*)' final solution'
      write(6,102)ampx(ic),ampy(ic),ampz(ic)
      write(6,102)tmat(1,ilm),tmat(2,ilm),tmat(3,ilm)
      write(6,103)kchan(ic),cs(ic,1),cs(ic,2),cs(ic,3)
      end if
c
 102  format(6e14.6)
 103  format(4e14.6)
      totx=totx+cs(ic,1)
      toty=toty+cs(ic,2)
      totz=totz+cs(ic,3)
c
 111  continue
c
c
c write transition amplitudes to a punch file
c
      write(11) 1,ic+1,nww,ilm,1.,kchan(ic)
      write(11)((tmat(iww(it),jt),it=1,nww),jt=1,ilm)
c
      if(iprt.gt.0) then
      write(6,*)' TMAT on file ic ilm ',ic,ilm
      write(6,88004)((it,iww(it),tmat(iww(it),jt),it=1,nww),jt=1,ilm)
88004 format(3(3x,2i5,4x,f12.8,2x,f12.8))
      end if
c
c
c      totx=16*pi*omega/3./137.*totx*100.*.529*.529
c      toty=16*pi*omega/3./137.*toty*100.*.529*.529
c      totz=16*pi*omega/3./137.*totz*100.*.529*.529
c
      totx=8*pi*omega/3./137.*totx*100.*.529*.529
      toty=8*pi*omega/3./137.*toty*100.*.529*.529
      totz=8*pi*omega/3./137.*totz*100.*.529*.529
      write(6,88001)iene,ic
88001 format(/,' ** total cross sections ** ',/,
     $         ' energy ',i5,' ion channel ',i5)
      write(6,88002)totx,toty,totz
      write(6,88003)tx,ty,tz
88002 format(3(2x,e15.7),//)
88003 format(' Trace TMAT ',/,3(2x,e15.7),//)
      ome(iene,ic)=omega
      crx(iene,ic)=totx
      cry(iene,ic)=toty
      crz(iene,ic)=totz
504   continue
c
c
c close big loop on incident energies
c
 1000 continue
      do 1001 i=1,nchan
      write(7,1003) i
 1001 write(7,1002)(ome(j,i),27.211*ome(j,i),
     $ crx(j,i),cry(j,i),crz(j,i),
     $ j=1,nener)
 1002 format(5f12.5)
 1003 format(/,' Channel  ',i5)
      stop
      end
      function carg(z)
      complex*16 z
      real*8 carg
      carg=atan2(imag(z),dble(z))
      return
      end
      double complex function cgamma(x)
      implicit real*8 (a - h, o - z)
      complex*16 x
      parameter (
     &    pi = 3.14159265358979324d+00, 
     &    pv = 7.31790632447016203d+00, 
     &    pu = 3.48064577727581257d+00, 
     &    pr = 3.27673720261526849d-02, 
     &    p1 = 1.05400280458730808d+01, 
     &    p2 = 4.73821439163096063d+01, 
     &    p3 = 9.11395751189899762d+01, 
     &    p4 = 6.62756400966213521d+01, 
     &    p5 = 1.32280130755055088d+01, 
     &    p6 = 2.93729529320536228d-01)
      parameter (
     &    q1 = 9.99999999999975753d-01, 
     &    q2 = 2.00000000000603851d+00, 
     &    q3 = 2.99999999944915534d+00, 
     &    q4 = 4.00000003016801681d+00, 
     &    q5 = 4.99999857982434025d+00, 
     &    q6 = 6.00009857740312429d+00)
      xr = dreal(x)
      xi = dimag(x)
      if (xr .lt. 0) then
          wr = 1 - xr
          wi = -xi
      else
          wr = xr
          wi = xi
      end if
      ur = wr + q6
      vr = ur * (wr + q5) - wi * wi
      vi = wi * (wr + q5) + ur * wi
      yr = p6 + (p5 * ur + p4 * vr)
      yi = p5 * wi + p4 * vi
      ur = vr * (wr + q4) - vi * wi
      ui = vi * (wr + q4) + vr * wi
      vr = ur * (wr + q3) - ui * wi
      vi = ui * (wr + q3) + ur * wi
      yr = yr + (p3 * ur + p2 * vr)
      yi = yi + (p3 * ui + p2 * vi)
      ur = vr * (wr + q2) - vi * wi
      ui = vi * (wr + q2) + vr * wi
      vr = ur * (wr + q1) - ui * wi
      vi = ui * (wr + q1) + ur * wi
      yr = yr + (p1 * ur + vr)
      yi = yi + (p1 * ui + vi)
      ur = vr * wr - vi * wi
      ui = vi * wr + vr * wi
      t = ur * ur + ui * ui
      vr = (yr * ur + yi * ui) + pr * t
      vi = yi * ur - yr * ui
      yr = wr + pv
      ur = 0.5d0 * log(yr * yr + wi * wi) - 1
      ui = atan2(wi, yr)
      yr = exp(ur * (wr - 0.5d0) - ui * wi - pu) / t
      yi = ui * (wr - 0.5d0) + ur * wi
      ur = yr * cos(yi)
      ui = yr * sin(yi)
      yr = ur * vr - ui * vi
      yi = ui * vr + ur * vi
      if (xr .lt. 0) then
          wr = pi * xr
          wi = exp(pi * xi)
          vi = 1 / wi
          ur = (vi + wi) * sin(wr)
          ui = (vi - wi) * cos(wr)
          vr = ur * yr + ui * yi
          vi = ui * yr - ur * yi
          ur = 2 * pi / (vr * vr + vi * vi)
          yr = ur * vr
          yi = ur * vi
      end if
      cgamma = cmplx(yr, yi)
      return
      end

