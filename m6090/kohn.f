*deck kohn
c***begin prologue     kohn
c***date written       920417   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, link 6080
c***author             schneider, barry (lanl)
c***source             m6090
c***purpose            one channel, s wave kohn code
c***description        kohn calculation for radial potential and s waves
c***                   where all integrals are performed numerically.
c***                   primarily for a test of useful variants of the kohn principle
c*** 
c
c***references       
c
c***routines called
c***end prologue       kohn
      program kohn
      implicit integer(a-z)
      common z(1)
      dimension a(1)
      equivalence (a,z)
      common /io/ inp, iout
      common /memory/ioff
      real *8 z, eryd, fpkey, rgs, wtt, charge
      real *8 rfs, rmatch, rk, rstr, rfnr
      character *4096 ops
      character *3 itoc
      character *80 cpass
      character *8 chrkey, mshtyp
      character *1600 card
      character *80 title
      character *16 vtyp
      logical logkey, typcal, prntbs, prntmt, prntbb, prntfb
      logical prntff, posinp
      dimension rgs(1000), wtt(1000)
      dimension eryd(100)
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      prntbs=logkey(ops,'print=m6090=basis-information',.false.,' ')
      prntmt=logkey(ops,'print=m6090=matrix-information',.false.,' ')
      prntbb=logkey(ops,'print=m6090=bound-bound',.false.,' ')
      prntfb=logkey(ops,'print=m6090=free-bound',.false.,' ')
      prntff=logkey(ops,'print=m6090=free-free',.false.,' ')
      typcal=logkey(ops,'m6090=fitted-forms',.false.,' ')
      write (iout,100)
c
      if ( posinp('$kohn',cpass) ) then
           call cardin (card)
           nbf=intkey(card,'no-internal-basis-functions',1,' ') 
           nbfn=intkey(card,'no-fitting-functions',2,' ')
           nobas=max(nbfn,nbf)
           charge=fpkey(card,'charge',0.d+00,' ')
           rmatch=fpkey(card,'matching-radius',10.d0,' ')
           nrg=intkey(card,'no-integration-regions',1,' ')
      endif
      if ( posinp ('$energy',cpass) ) then
           call cardin (card)
           nen=intkey(card,'no-energies',1,' ')
           call fparr (card,'energies',eryd,nen,' ')
      endif
c
      write (iout,200)
      write (iout,300) nbf, charge, nrg, nen
      if (typcal) then
          write(iout,310) nbfn
      else
          write(iout,320)
      endif
*
      call rzero(rgs,1000)
      call rzero(wtt,1000)
c     ----- nrg=no. regions of integration , nps=no. points in a region
c     ----- rstart=starting value for a region , rfin=final value -----
c     ----- mshtyp=0 is gauss-legendre quadrature =1 trapezoidal rule --
c     ----- rgs and wtt hold the point and integration weights -----
c
      ic=0
      do 10 iz=1,nrg
         if ( posinp('$region-'//itoc(iz),cpass) ) then
              call cardin(card)
              mshtyp=chrkey(card,'type-quad','legendre',' ')
              npnr=intkey (card,'no-points',3,' ')
              rstr=fpkey (card,'starting-r',0.d0,' ')
              rfnr=fpkey (card,'ending-r',10.d0,' ')
              write(iout,400) iz,npnr,rstr,rfnr
              rfs=rfnr-rstr
              nps=npnr
              do 20 jz=1,nps
                 ic=ic+1
                 call lgndrx (nps,jz,wtt(ic),rgs(ic))
                 wtt(ic)=wtt(ic)*rfs
                 rgs(ic)=rgs(ic)*rfs+rstr
   20         continue
         endif
   10 continue
      npnts=ic
      if (typcal) then
          do 40 i=1,npnts
             if (rgs(i).ge.rmatch) go to 50
   40    continue
         lasti=npnts
         go to 60 
   50    lasti=i-1      
   60    continue
         write(iout,500) lasti, rmatch
      endif  
c----------------------------------------------------------------------c
c                      memory allocation                               c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,maxsiz,0,' ')
      words=5*nbfn+2*npnts*nbfn+8*npnts+10000
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')      
      call getscm(words,z(1),ngot,'m6080',0)
      alpha=ioff
      power=wpadti(alpha+nobas)
      ifns=iadtwp(power+nobas)
      iddfns=ifns+npnts*nbf
      vmat=iddfns+npnts*nbf
      if(typcal) then
         fns=vmat+npnts
         ddfns=iddfns+npnts*nbfn
         mtch=ddfns+nbfn*npnts
         mtrx=mtch+3*nbfn
         rhs=mtrx+18
         diff=rhs+6
         cfn=diff+2*npnts
      else
         cfn=vmat+npntsn
      endif
      ddcfn=cfn+4*npnts
      hbbp=ddcfn+4*npnts
      hfbp=hbbp+nbf*nbf
      hffp=hfbp+4*nbf
      sbbp=hffp+2
      sfbp=sbbp+nbf*nbf
      hbbt=sfbp+4*nbf
      hfbt=hbbt+nbf*nbf
      eig=hfbt+4*nbf
      work=eig+nbf
      scrb=work+nbf
      scrc=scrb
      temp=scrb
      tempc=scrb
      mtrx=scrb
      rhs=mtrx+2*(nbf+1)*(nbf+1)
      ipvt=wpadti(rhs+2*(nbf+1))
      scrwds=max(nbf*npnts,4*npnts,nbf*nbf,4*nbf)
      addwds=2*(nbf+1)*(nbf+1)+3*(nbf+1)
      words=scrb+max(scrwds,addwds)
      if(  posinp('$basis',cpass) ) then
           call basis(z(ifns),z(iddfns),z(alpha),a(power),rgs,npnts,nbf)
           title='internal basis functions'
           call prntrm(title,z(ifns),npnts,nbf,npnts,nbf,iout)
           title='internal basis functions second derivatives'
           call prntrm(title,z(iddfns),npnts,nbf,npnts,nbf,iout)
      endif
      if (typcal) then
          if ( posinp('$fitfn',cpass) ) then
               call fitbas(z(fns),z(ddfns),z(alpha),a(power),rgs,
     1                     z(mtch),rmatch,npnts,nbfn)      
               if (prntbs) then
                   title='fitting functions'
                   call prntrm(title,z(fns),npnts,nbfn,npnts,nbfn,iout)
                   title='fitting function second derivatives'
                   call prntrm(title,z(ddfns),npnts,nbfn,npnts,nbfn,iout)
               endif
          endif
      endif
c
      call modpot(z(vmat),rgs,npnts,charge,vtyp)
      call hbound(z(ifns),z(iddfns),z(vmat),wtt,z(hbbp),z(sbbp),z(hbbt),
     1            z(eig),z(work),z(temp),z(scrb),nbf,npts,count,prntbb)
      do 70 ien=1,nen
         rk=eryd(ien)*.5d0
         write(iout,600) ien, eryd(ien)
         if (typcal) then
             call match(z(fns),z(ddfns),z(cfn),z(ddcfn),z(diff),
     1                  rgs,z(mtch),z(mtrx),z(rhs),rk,rmatch,
     2                  lasti,a(ipvt),npnts,nbfn)
         else
             call freefn(z(cfn),z(ddcfn),rgs,rk,npnts)
         endif
         if (prntbs) then
             title='complex fitting function'
             call prntcmn(title,z(cfn),npnts,2,npnts,2,iout,'e')
             title='complex fitting function second derivatives'
             call prntcmn(title,z(ddcfn),npnts,2,npnts,2,iout,'e')
             if (typcal) then
                 title='difference'
                 call prntcmn(title,z(diff),npnts,1,npnts,1,iout,'e')
             endif
         endif
         call hsmat(z(ifns),z(iddfns),z(cfn),z(ddcfn),z(vmat),rk,wtt,
     1              z(hbbp),z(hfbp),z(hffp),z(sbbp),z(sfbp),z(hbbt),
     2              z(hfbt),z(scrb),z(scrc),z(temp),z(tempc),nbf,
     3              npts,count,prntfb,prntff)
         call tophse(z(hbbt),z(hfbt),z(hffp),rk,z(mtrx),z(rhs),
     1               a(ipvt),nbf)
   70 continue
      call iosys ('write integer maxsiz to rwf',1,maxsiz,0,' ') 

  100 format (//,25x,'***** one channel kohn calculation *****')
  200 format (///,25x,'***** input data *****')
  300 format(//,5x,'no. basis functions',1x,i3,1x,'charge',1x,e15.8,/,
     1          5x,'no. integration regions',1x,i3,1x,'no. energies',
     2          1x,i3)
  310 format(//,25x,'using cutoff functions',/,10x,'will fit to',
     1               1x,i2,1x,'functions') 
  320 format(//,25x,'using standard regular and irregular forms') 
  400 format (//,25x,'region',1x,i3,/,5x,'no. points',1x,i4,1x,
     1               'r left',1x,e15.8,1x,'r right',1x,e15.8)
  500 format(//,15x,'using box matched free functions',//,
     1          5x,'matching point',1x,i4,1x,'value',1x,e15.8)
  600 format(/,10x,'the',1x,i3,1x,'energy = ',e15.8)
      call chainx (0)
      stop
      end










