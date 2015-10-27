*deck @(#)meopt.f
c***begin prologue     meopt
c***date written       920409   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6060, link 6060
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            overlap matrix elements between vlamdas of complex
c***                   optical potential, the kohn free functions and
c***                   the gaussians.
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6060
      program meopt
      implicit integer (a-z)
      parameter (dimc=30 , dime=100 , dimlm=100 , dimbf=200, dimmo=200)
      real*8 z, echan, energy, kchan, rmin, rmax, rdel, alpha, gamma
      real *8 rd26, ec
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /io / inp, iout
      common / memory / ioff
c----------------------------------------------------------------------c
c        the small arrays below are the only explicitly dimensioned    c
c             arrays in the code (see parameter statement)             c
c----------------------------------------------------------------------c
      dimension echan(dimc), energy(dime), kchan(dimc), nlm(dimc)
      dimension lch(dimlm,dimc), mch(dimlm,dimc), ngauss(dimc)
      dimension ngch(dimbf,dimc), nmoc(dimc), nmoch(dimbf,dimc)
      dimension list(dimbf), ntrgt(dimc), nbscat(dimc)
      dimension orblst(dimmo,dimc), finlst(dimmo,dimc), nbtot(dimc)
      logical logkey, ptest, bsym, typ, posinp, north, rtest, zroset
      character *4096 ops
      character *1600 card
      character *24 filnm
      character *32 xform
      character *13 grdtyp
      character *128 filgrd, filorb, filylm, filbes
      character *128  filint, filkne
      character *13 chrkey, srtfl
      character *10 cpass
      character*8  bcondx
      character *16 fptoc
      character *3 itoc
      data pi /3.14159265358979323846d+00/
c----------------------------------------------------------------------c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      ptest=logkey(ops,'print=m6060',.false.,' ')
      rtest=logkey(ops,'real-check',.false.,' ')
      zroset=logkey(ops,'zero-potential',.false.,' ')
      north=logkey(ops,'no-orthogonality',.false.,' ')
      srtfl='"sorted orbs"'
      typ=.true.
      bcondx='t-matrix'
      filgrd='grid'
      filorb='orbs'
      filylm='ylm'
      filbes='bessfn'
      filint='optint'
      if( posinp('$optint',cpass) ) then
          call cardin(card)
	  nchan=intkey(card,'no-channels',1,' ')
          srtfl=chrkey(card,'orbital-file-type','con array',' ')
          if (srtfl.eq.'con array') then
              srtfl='"con array"'
              typ=.false.
          else
              srtfl='"sorted orbs"'
              typ=.true.
          endif
c----------------------------------------------------------------------c
c        read boundary condition specification , filenames etc         c
c----------------------------------------------------------------------c
          bcondx=chrkey(card,'boundary-condition','t-matrix',' ')
      endif
      call iosys ('read character "grid filename" from rwf',-1,
     1             0,0,filgrd)
      call iosys ('read character "orbital filename" from rwf',-1,
     1             0,0,filorb)
      call iosys ('read character "spherical harmonic filename" '//
     1            'from rwf',-1,0,0,filylm)
      call iosys ('read character "bessel function filename" '//
     1            'from rwf',-1,0,0,filbes)
      call iosys ('read character "optical potential filename" '//
     1            'from rwf',-1,0,0,filint)         
      call iosys ('read character "kohn data filename" from rwf',0,0,0,
     1             filkne)
      write (iout,1)
      call iosys ('open kohndt as old',0,0,0,filkne)
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('open optint as old',0,0,0,filint)
      call iosys ('open bessel as old',0,0,0,filbes)
      call iosys ('open orbs as old',0,0,0,filorb)
      call iosys ('open ylms as old',0,0,0,filylm)
      call iosys ('write character "boundary cond" to optint',0,0,
     1             0,bcondx)
c----------------------------------------------------------------------c
c            read in spline information for free functions             c
c----------------------------------------------------------------------c
      call iosys ('read integer maximum-l from bessel',1,lmax,0,' ')
      call iosys ('read integer "total pts" from bessel',1,nr,0,' ')
      call iosys ('read real rmin from bessel',1,rmin,0,' ')
      call iosys ('read real rmax from bessel',1,rmax,0,' ')
      call iosys ('read real "r spacing" from bessel',1,rdel,0,' ')
      rd26= rdel* rdel/6.
c----------------------------------------------------------------------c
c             read grid parameters from basis set file                 c
c----------------------------------------------------------------------c
      call iosys ('read character "grid type" from orbs',0,0,0,grdtyp)
      call iosys ('read integer "no. grid pts" from orbs',1,ngrid,0,
     1            ' ')
      call iosys ('read integer "point buffer" from orbs',1,pntbuf,
     1            0,' ')
      call iosys ('read integer "no. cont" from orbs',1,ncon,0,' ')
      call iosys ('read integer "no. regions" from orbs',1,nreg,0,' ')
      call iosys ('read integer "final pts" from orbs',1,nolst,0,' ')
      nkept=ncon
      do 100 i=1,nkept
         list(i)=i
  100 continue
      if (typ) then
          write (iout,1500)
          call iosys ('read integer "no. kept" from orbs',1,nkept,
     1                0,' ')
          call iosys ('read integer "function list" from orbs',dimbf,
     1                list,0,' ')
      endif
c----------------------------------------------------------------------c
c              read parameters from ylm file and check                 c
c----------------------------------------------------------------------c
      call iosys ('read integer "no. grid pts" from ylms',1,ngchk,0,
     1            ' ')
      call iosys ('read integer "point buffer" from ylms',1,pntchk,0,
     1            ' ')
      call iosys ('read integer "max l in ylm" from ylms',1,lmxylm,0,
     1            ' ')
      call iosys ('read integer "max m in ylm" from ylms',1,mumax,0,
     1            ' ')
      call iosys ('read integer "no. vlamdas" from optint',1,nolam,
     1             0,' ')
c----------------------------------------------------------------------c
c         read in information about energies and channels              c
c         the scattering energies are incident electron energies       c
c----------------------------------------------------------------------c
      call iosys ('read integer "no. open chan" from kohndt',1,nchan,
     1             0,' ')
      ntri=nchan*(nchan+1)/2
      call iosys ('read integer "no. energies" from kohndt',1,
     1             nener,0,' ')
      call iosys ('read real "chan energies" from kohndt',nchan,
     1             echan,0,' ')
      call iosys ('read real "scatt energies" from kohndt',nener,
     1             energy,0,' ')
      call iosys ('read integer "total no. mos" from kohndt',1,nmotot,
     1             0,' ')
      if( posinp('$energy',cpass) ) then
          call cardin(card)
          call fparr(card,'channel-energies',echan,nchan,' ')
          nener=intkey(card,'no-energies',1,' ')
          call fparr(card,'scattering-energies',energy,nener,' ')
      endif
c     ignd=ismin(nchan,echan,1)
      write (iout,3) nchan,(echan(i),i=1,nchan)
      write (iout,5) rmin, rmax, rdel, lmax, alpha, gamma, ncut
      write (iout,6) bcondx
c----------------------------------------------------------------------c
c             read channel quantum numbers and assignment of           c
c             l**2 basis functions to channels                         c
c----------------------------------------------------------------------c
      if( posinp('$target',cpass) ) then
          call cardin(card)
          bsym=logkey(card,'bound-symmetry=off',.false.,' ')
          if (bsym) then
              call iosys ('write character symmetry to optint',0,0,
     1                     0,'on')
          else
              call iosys ('write character symmetry to optint',0,0,
     1                     0,'off')
          endif
          maxlm=0
          matbf=0
          matbv=0
          ntchn=0
          do 60 ch1=1,nchan
             if( posinp('$chan-'//itoc(ch1),cpass) ) then
                 call cardin(card)
                 nlm(ch1)=intkey(card,'no-lm',1,' ')
                 call intarr(card,'l-values',lch(1,ch1),nlm(ch1),' ')
                 call intarr(card,'m-values',mch(1,ch1),nlm(ch1),' ')
                 maxlm=max(maxlm,nlm(ch1))
                 ntchn=ntchn+nlm(ch1)
                 if (.not.bsym) then
                     ngauss(ch1)=intkey(card,'no-aos',1,' ')
                     call intarr(card,'aos',ngch(1,ch1),ngauss(ch1),' ')
                     nmoc(ch1)=intkey(card,'no-mos',1,' ')
                     call intarr(card,'mos',nmoch(1,ch1),nmoc(ch1),' ')
                     nbscat(ch1)=intkey(card,'no-expansion-mos',1,' ')
                     call intarr(card,'expansion-mos',orblst(1,ch1),
     1                           nbscat(ch1),' ')
                     ntrgt(ch1)=intkey(card,'no-target-mos',1,' ')
                     ii=nbscat(ch1)+1
                     call intarr(card,'target-mos',orblst(ii,ch1),
     1                           ntrgt(ch1),' ')
                 else
                     ngauss(ch1)=ncon
                     nbscat(ch1)=nmotot
                     nmoc(ch1)=ncon
                     ntrgt(ch1)=0
                     do 70 orb=1,nmotot
                        orblst(orb,ch1)=orb
   70                continue
                     do 80 orb=1,ncon
                        ngch(orb,ch1)=orb
                        nmoch(orb,ch1)=orb
   80                continue
                 endif
                 nbtot(ch1)=nbscat(ch1)+ntrgt(ch1)
                 matbv=matbv+nbscat(ch1)
                 matbf=matbf+nbtot(ch1)
                 write(iout,7) ch1
                 write (iout,8) (lch(i,ch1),i=1,nlm(ch1))
                 write (iout,9) (mch(i,ch1),i=1,nlm(ch1))
                 write (iout,10) (ngch(i,ch1),i=1,ngauss(ch1))
                 write (iout,25) (nmoch(i,ch1),i=1,nmoc(ch1))
                 write (iout,40) (orblst(i,ch1),i=1,nbscat(ch1))
                 if (ntrgt(ch1).ne.0) then
                     write (iout,50) (orblst(i,ch1),i=ii,nbscat(ch1)+
     1                                ntrgt(ch1))
                 endif
                 do 90 j=1,nlm(ch1)
                    mm=2*iabs(mch(j,ch1))
                    if (mch(j,ch1).ge.0) then
                        mch(j,ch1)=mm
                    else
                        mch(j,ch1)=mm-1
                    endif
   90            continue
               endif
   60     continue
      endif
      call newlst (nbtot,nbscat,ntrgt,finlst,nchan,dimc,dimmo)
      call iosys ('write integer "no. energies" to optint',1,
     1             nener,0,' ')
      call iosys ('write real "chan energies" to optint',nchan,
     1             echan,0,' ')
      call iosys ('write real "scatt energies" to optint',nener,
     1             energy,0,' ')
      call iosys ('write integer "no. channels" to optint',1,
     1             nchan,0,' ')
      call iosys ('write real eground to optint',1,echan(1),0,' ')
      call iosys ('write integer "lm-vals-chan" to optint',nchan,
     1             nlm,0,' ')
      call iosys ('write integer "no.-gauss-chan" to optint',nchan,
     1             ngauss,0,' ')
      call iosys ('write integer "no.-mos-chan" to optint',nchan,
     1             nmoc,0,' ')
      do 200 ch1=1,nchan
         filnm='"l val ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to optint',nlm(ch1),
     1                lch(1,ch1),0,' ')
         filnm='"m val ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to optint',nlm(ch1),
     1                mch(1,ch1),0,' ')
         filnm='"aos ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to optint',ngauss(ch1),
     1                ngch(1,ch1),0,' ')
         filnm='"mos ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to optint',nmoc(ch1),
     1                nmoch(1,ch1),0,' ')
  200 continue
      call iosys ('write integer "total no. channels" to optint',
     1             1,ntchn,0,' ')
      call iosys ('write integer "no. variational orbitals" to optint',
     1             1,matbv,0,' ')
c----------------------------------------------------------------------c
c               get memory and begin calculation one energy            c
c                           at a time                                  c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c   dimensions of bessel and spherical functions are determined by     c
c   the largest l and m values needed. the matrix elements are         c
c   dimensioned by channels and maximum number of l or basis functions c
c                               needed                                 c
c----------------------------------------------------------------------c
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      wdplm=2*maxlm*nolam*nchan
      wdmlm=wdplm
      wdlmb=2*nolam*nkept*nchan
      wdlmbm=2*nolam*nmotot*nchan
      wdpb=2*maxlm*nkept*nchan
      wdmb=2*maxlm*nkept*nchan
      wdpbm=2*maxlm*nmotot*nchan
      wdmbm=wdpbm
      beswdr=nr*(lmax+1)
      beswds=2*beswdr
      newwds=2*pntbuf*maxlm*nchan
      ylmwds=pntbuf*(lmax+1)*(mumax+mumax+1)
      baswds=nkept*pntbuf
      newwds=2*pntbuf*maxlm*nchan
      cmtwds=4*maxlm*pntbuf
      bmtwds=2*nkept*pntbuf
      scrwds=max(maxlm*nolam,nkept*nolam,maxlm*nkept,
     1           nolam*ncon,maxlm*ncon)
      scrwds=2*scrwds
      wdff=2*maxlm*maxlm
      wdfb=2*maxlm*nmotot
      wdout=6*ntchn*ntchn+4*ntchn*matbf+2*matbf*matbf
      wd1=iadtwp(lmax+1)+2*nolam*pntbuf+nr+2*beswdr+newwds+6*pntbuf+
     1      ylmwds+baswds+cmtwds+bmtwds+scrwds
      wd2=wdff*(2*ntri+nchan*nchan)+2*nmotot*nmotot*ntri+
     1          2*wdfb*nchan*nchan+wdout
      words=2*nolam+ncon*nmotot+wdplm+wdmlm+wdlmb+wdlmbm+wdpb+wdmb+
     1      wdpbm+wdmbm+max(wd1,wd2)
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6060',0)
      write (iout,300) words
      eigval=ioff
      vec=eigval+2*nolam
      ovplm=vec+ncon*nmotot
      ovmlm=ovplm+wdplm
      ovlmb=ovmlm+wdmlm
      ovpb=ovlmb+wdlmb
      ovmb=ovpb+wdpb
      ovlmbm=ovmb+wdmb
      ovpbm=ovlmbm+wdlmbm
      ovmbm=ovpbm+wdpbm
      vlamda=ovmbm+wdmbm
      ipow=wpadti(vlamda+2*nolam*pntbuf)      
      x=iadtwp(ipow+lmax+1)      
      hs=x+nr
      cj=hs+beswdr
      hp=cj+beswdr
      rvec=hp+newwds
      krvec=rvec+pntbuf
      grid=krvec+pntbuf
      ylm=grid+4*pntbuf
      basis=ylm+ylmwds
      cmat=basis+baswds
      bmat=cmat+cmtwds
      scrc=bmat+bmtwds
c----------------------------------------------------------------------c
c                   core reallocated here                              c
c----------------------------------------------------------------------c
      vpp=vlamda
      vpm=vpp+wdff*ntri
      vmm=vpm+wdff*nchan*nchan
      vbb=vmm+wdff*ntri
      vbp=vbb+2*nmotot*nmotot*ntri
      vbm=vbp+wdfb*nchan*nchan      
      vppo=vbm+wdfb*nchan*nchan
      vpmo=vppo+2*ntchn*ntchn
      vmmo=vpmo+2*ntchn*ntchn
      vpbo=vmmo+2*ntchn*ntchn
      vmbo=vpbo+2*ntchn*matbf
      vbbo=vmbo+2*ntchn*matbf
      call iosys ('read real "complex eigenvalues" from optint',2*nolam,
     1             z(eigval),0,' ')
      if (rtest) then
          call zimag(z(eigval),nolam,1)
      endif
      call mkpow(a(ipow),lmax)
      call iosys('read character "transformation vector" from kohndt',
     $           -1,0,0,xform)
      call iosys ('read real '//xform//' from kohndt',ncon*nmotot,
     1             z(vec),0,' ')
      call iosys ('read real points from bessel',nr,z(x),0,' ')
      call iosys('read real "complex kohn function" from bessel',
     1            beswds,z(hs),0,' ')            
      call iosys('read real "function spline coefficients" from '//
     1           'bessel',beswds,z(cj),0,' ')            
      if (zroset) then
          write(iout,*) '   zeroing potential'
      endif
c----------------------------------------------------------------------c
c                      loop over energies                              c
c----------------------------------------------------------------------c
      do 400 ene=1,nener
c----------------------------------------------------------------------c
c                 construct channel momenta                            c
c----------------------------------------------------------------------c
         do 410 ch1=1,nchan
            ec = .5d+00*energy(ene) - (echan(ch1)-echan(1))
            if(ec.le.0.0) then
               call lnkerr('no closed channels allowed at present:quit')
            endif
            kchan(ch1) = sqrt(2.0*ec)
  410    continue
         write(iout,420) energy(ene), (kchan(i),i=1,nchan)
c----------------------------------------------------------------------c
c              rewind all energy-independent files                     c
c----------------------------------------------------------------------c
         call iosys ('rewind all on grid read-and-write',0,0,0,' ')
         call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
         call iosys ('rewind all on optint read-and-write',0,0,0,' ')
         call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
c----------------------------------------------------------------------c
c                   zero matrix elements                               c
c----------------------------------------------------------------------c
         call rzero(z(ovplm),wdplm)
         call rzero(z(ovmlm),wdmlm)
         call rzero(z(ovlmb),wdlmb)
         call rzero(z(ovlmbm),wdlmbm)
         call rzero(z(ovpb),wdpb)
         call rzero(z(ovmb),wdmb)
         call rzero(z(ovpbm),wdpbm)
         call rzero(z(ovmbm),wdmbm)
c----------------------------------------------------------------------c
c                   loop over grid                                     c
c----------------------------------------------------------------------c
         npnts=pntbuf
         wrdylm=pntbuf*(lmax+1)
         do 430 ireg=1,nreg
            if (ireg.eq.nreg) then
                npnts=nolst
                wrdylm=npnts*(lmax+1)
            endif
            call iosys ('read real '//grdtyp//' from grid without '//
     1                  'rewinding',4*npnts,z(grid),0,' ')
c----------------------------------------------------------------------c
c                read in a block of gaussians                          c
c----------------------------------------------------------------------c
            call iosys ('read real '//srtfl//' from orbs without '//
     1                  'rewinding',npnts*nkept,z(basis),0,' ')
c----------------------------------------------------------------------c
c                read in a block of ylm's                              c
c----------------------------------------------------------------------c
            ylmc=ylm
            do 440 m=0,mumax
               if(m.eq.0) then
                  call iosys ('read real ylm from ylms without '//
     1                        'rewinding',wrdylm,z(ylmc),0,' ')
                  ylmc=ylmc+wrdylm
               else
                  call iosys ('read real ylm from ylms without '//
     1                        'rewinding',wrdylm,z(ylmc),0,' ')
                  ylmc=ylmc+wrdylm
                  call iosys ('read real ylm from ylms without '//
     1                        'rewinding',wrdylm,z(ylmc),0,' ')
                  ylmc=ylmc+wrdylm
               endif
  440       continue
c----------------------------------------------------------------------c
c                compute the bessel functions from the spline          c
c----------------------------------------------------------------------c
            call mbes (z(hs),z(cj),z(hp),z(grid),z(rvec),z(x),
     1                 z(krvec),kchan,rmin,rdel,rd26,a(ipow),nlm,lch,
     2                 npnts,nr,lmax,maxlm,nchan,dimlm,dimc)
c----------------------------------------------------------------------c
c              read vlamdas for this block                             c
c----------------------------------------------------------------------c
            call iosys ('read real "complex vlamdas" from optint '//
     1                  'without rewinding',2*npnts*nolam,z(vlamda),
     2                   0,' ')
            if (rtest) then
                call zimag(z(vlamda),npnts,nolam)
            endif
c----------------------------------------------------------------------c
c             multiply integration weight into vlamdas                 c
c----------------------------------------------------------------------c
            call sclvlm(z(vlamda),z(grid),npnts,nolam)
c----------------------------------------------------------------------c
c               calculate overlaps with vlamdas                        c
c----------------------------------------------------------------------c
            call ovrlp(z(ovplm),z(ovmlm),z(ovlmb),z(ovpb),z(ovmb),
     1                 z(basis),z(hp),z(ylm),z(vlamda),z(grid),z(cmat),
     2                 z(bmat),z(scrc),nlm,lch,mch,ngauss,ngch,npnts,
     3                 nchan,nkept,list,lmax,dimlm,dimc,dimbf,maxlm,
     4                 nolam,bcondx)
  430    continue

         filnm='ch-en-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',nchan,kchan,
     1               0,' ')
c-----------------------------------------------------------------------c
c              transform to molecular orbital basis                     c
c-----------------------------------------------------------------------c
         if(ptest) then
            call aoout(z(ovplm),z(ovmlm),z(ovlmb),z(ovpb),z(ovmb),
     1                 energy(ene),nchan,maxlm,nolam,nkept)
         endif
         call tomo(z(vec),z(ovlmb),z(ovlmbm),z(ovpb),z(ovpbm),z(ovmb),
     1             z(ovmbm),z(scrc),z(scrc),maxlm,nolam,ncon,nmotot,
     2             nchan,ngauss,ngch,nkept,list,dimbf,dimc)
         filnm='plam-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',wdplm,z(ovplm),
     1                0,' ')
         filnm='mlam-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',wdmlm,z(ovmlm),
     1                0,' ')
         filnm='blam-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',wdlmbm,
     1                z(ovlmbm),0,' ')
         filnm='pbm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',wdpbm,
     1                z(ovpbm),0,' ')
         filnm='mbm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',wdmbm,
     1                z(ovmbm),0,' ')
c----------------------------------------------------------------------c
c                    output section                                    c
c----------------------------------------------------------------------c
         if(ptest) then
            call moout(z(ovlmbm),z(ovpbm),z(ovmbm),energy(ene),nchan,
     1                 maxlm,nolam,nmotot)
         endif
c----------------------------------------------------------------------c
c                construct optical potential matrix elements           c
c----------------------------------------------------------------------c
         call optpot(z(ovplm),z(ovmlm),z(ovlmbm),z(ovpbm),z(ovmbm),
     1               z(vpp),z(vpm),z(vmm),z(vbb),z(vbp),z(vbm),
     2               energy(ene),z(eigval),nchan,maxlm,nolam,
     3               nmotot,ntri,north)
c----------------------------------------------------------------------c
c               put in channel representation                          c
c----------------------------------------------------------------------c
         call optmat(z(vpp),z(vpm),z(vmm),z(vbp),z(vbm),z(vbb),z(vppo),
     1              z(vpmo),z(vmmo),z(vpbo),z(vmbo),z(vbbo),nbscat,
     2              orblst,finlst,nlm,nchan,ntchn,maxlm,nmotot,matbv,
     3              itri,dimmo,dimc,ptest,zroset)
         filnm='vppo-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',2*ntchn*ntchn,
     1                z(vppo),0,' ')
         filnm='vpmo-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',2*ntchn*ntchn,
     1                z(vpmo),0,' ')
         filnm='vmmo-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',2*ntchn*ntchn,
     1                z(vmmo),0,' ')
         filnm='vpbo-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',2*ntchn*matbv,
     1                z(vpbo),0,' ')
         filnm='vmbo-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',2*ntchn*matbv,
     1                z(vmbo),0,' ')
         filnm='vbbo-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to optint',2*matbv*matbv,
     1                z(vbbo),0,' ')
c close big loop on incident energies
c
  400 continue
c----------------------------------------------------------------------c
c                   rewind and close all files                         c
c----------------------------------------------------------------------c
      call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('close kohndt',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('rewind all on optint read-and-write',0,0,0,' ')
      call iosys ('close optint',0,0,0,' ')
      call iosys ('rewind all on bessel read-and-write',0,0,0,' ')
      call iosys ('close bessel',0,0,0,' ')
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
      call iosys ('close ylms',0,0,0,' ')
      call iosys ('rewind all on optint read-and-write',0,0,0,' ')
      call iosys ('close optint',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
    1 format(//,20x,'*****  m6060:optical potential matrix elements ***
     1**')
    3 format(/,5x,'target energies for ',i3,' channels:',/,(2x,5e15.8))
    4 format(/,5x,'incident energies:',/,(2x,5e15.8))
    5 format(/,5x,'spline data',/,5x,'rmin',1x,f10.5,1x,'rmax',f10.5,
     11x,'step size',1x,f10.5,/,5x,'max. l',1x,i3,1x,'inner cutoff',1x,
     2 f6.3,1x,'outer exponential cutoff',1x,f6.3,1x,'outer n cutoff',
     3 1x,i3)
    6 format(/,5x,'boundary conditions are',1x,a10)
    7 format(/,5x,'information for channel',1x,i3)
    8 format(/,5x,'l-values',(/,15x,10(i2,1x)))
    9 format(/,5x,'m-values',(/,15x,10(i2,1x)))
   10 format(/,5x,'bound-atomic-orbitals',(/,17x,10(i3,1x)))
   25 format(/,5x,'bound-molecular-orbitals',(/,17x,10(i3,1x)))
   40 format(/,5x,'expansion molecular orbitals',(/,17x,10(i3,1x)))
   50 format(/,5x,'target molecular orbitals',(/,17x,10(i3,1x)))
  300 format(/,5x,'need',1x,i8,1x,'words for calculation')
  420 format(/,5x,'incident energy(ry)',1x,f10.6,/,5x,'channel momenta',
     1       (/,10x,5(e12.5,1x)))
 1000 format(/,5x,'atomic orbital representation')
 1500 format(/,5x, 'will use short orbital list')
 2000 format(/,5x,'molecular orbital representation')
      call chainx(0)
      stop
      end




