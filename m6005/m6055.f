*Deck M6005
c***begin prologue     M6005
c***date written       XXXXXX   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6005, link 6005, Orbital Decomposition
c***author             Schneider, Barry (lanl)
c***source             m6005
c***purpose            Free-Free And Bound-Free Matrix Elements
c***description        Calculation Of Free-Free And Bound-Free Integrals
c***                   Of Unperturbed Hamiltonian Needed In Kohn
c***                   Calculations.
c***
c***                   For Free-Free Integrals:Computes M = (H - 0.5*kchan**2)
c***                                                      +     +
c***                   Matrix Elements Of Two Classes: < h  M  h > And
c***                      +     -                      -
c***                   < h  M  h >. The Definition Of h Changes Depending
c***                   On whether T-Matrix Or S-Matrix Boundary Conditions
c**                    Are Specified. For Bound-Free Integrals:Computes
c***                      +
c***                   < h (H-E) b >
c***                   For T-Matrix Boundary Conditions The Imaginary Part
c***                                                     -
c***                   Of The Free Function Is Used For h  Which Is
c***                    Asymptotically A Regular Bessel Function.
c***references
c
c***routines called    IOsys, Util and MDutil
c***end prologue       M6005
      program ffbf
      implicit integer (a-z)
      parameter (dimc=30 , dime=100 , dimlm=100 , dimbf=200)
      real*8 z, echan, energy, kchan, rmin, rmax, rdel, alpha, gamma
      real *8 rd26, ec
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /io / inp, iout
      common / memory / ioff
c----------------------------------------------------------------------c
c        The Small Arrays Below Are The Only Explicitly Dimensioned    c
c             Arrays In The Code (See Parameter Statement)             c
c----------------------------------------------------------------------c
      dimension echan(dimc), energy(dime), kchan(dimc), nlm(dimc)
      dimension lch(dimlm,dimc), mch(dimlm,dimc), ngauss(dimc)
      dimension ngch(dimbf,dimc), nmoc(dimc), nmoch(dimbf,dimc)
      dimension list(dimbf)
      logical logkey, ptest, bsym, typ, posinp
      logical prntbf, prntff, prntv, prnty, prntb, prntbs, prntgr
      character *4096 ops
      character *1600 card
      character *16 filnm
      character *32 xform
      character *3 ans
      character *13 grdtyp, grdchk
      character *8 filgrd, filorb, filylm, filpot, filkne
      character *13 chrkey, srtfl
      character *10 cpass
      character*8 filbes, filint, bcondx
      character *10 fptoc
      character *3 itoc
c----------------------------------------------------------------------c
c               Ancient History                                        c
c               Code Now Dynamically Dimensioned                       c
c----------------------------------------------------------------------c
c     parameter (nblok=1000,mtop=5,ltop=5,nchnl=2)                     c
c     parameter (nsblok=1100)                                          c
c     parameter (lmtop=(2*mtop+1)*(ltop+1)-mtop*(mtop+1))              c
c     parameter (nstate=nchnl**2)                                      c
c     parameter (nsymst = nchnl*(nchnl+1)/2)                           c
c     parameter (nbfmax=50)                                            c
c     parameter (maxene=20)                                            c
c  nblok = number of points in block of grid etc. values               c
c  mtop = largest m value                                              c
c  ltop = largest l value                                              c
c  nchnl = largest number of channels allowed                          c
c  nbfmax  =  largest number of l**2 basis functions allowed           c
c      (transformation to MO's occurs in scattering code)              c
c  nsblok = max number of spline points for bessel fcns                c
c  maxene = maximum # of incident energies                             c
c----------------------------------------------------------------------c
c      complex hp(nblok,0:ltop,nchnl), hd(nblok,0:ltop,nchnl)
c      complex cj(nsblok,0:ltop), cy(nsblok,0:ltop)
c      complex hs(nsblok,0:ltop), hsder(nsblok,0:ltop)
c      real cjr(2*nsblok,0:ltop), cyr(2*nsblok,0:ltop)
c      real hsr(2*nsblok,0:ltop), hsdr(2*nsblok,0:ltop)
c      real rvec(nblok),krvec(nblok), x(nsblok)
c      integer ipow(0:ltop)
c      real grid(nblok,3)
c      real vbuf(nblok*nchnl*(nchnl+1)/2)
c      real vpot(nblok,nsymst)
c      real wt(nblok)
c      real echan(nchnl), energy(maxene)
c      real kchan(nchnl)
c      real ylm(nblok,0:ltop,0:2*ltop)
c      real buff(nblok*4),biff(2*(nblok*(ltop+1)))
c      real basis(nbfmax*nblok)
c      dimension ibuff(10000b,3)
c      dimension nlm(nchnl),lch(lmtop,nchnl),mch(lmtop,nchnl)
c      dimension ngauss(nchnl),ngch(nbfmax,nchnl)
c      complex cdotu,cvec(nblok),hpvb(lmtop,nbfmax,nstate)
c      complex ovbf(lmtop,nbfmax,nchnl)
c      complex hpvhp(lmtop,lmtop,nsymst),hpvhm(lmtop,lmtop,nchnl**2)
c      complex evec(nblok), dvec(nblok)
c      common/ivparms/ istate,ngrid3
c      common/parms/mgrid,npt,nbf
c      common/parmbes/ngrid2,nchunk,nener,nchan,lmax2
c      common/parm/nbig,mpt,lmax,mumax
       data pi /3.14159265358979323846d+00/
c      equivalence (cy(1),cyr(1)),(cj(1),cjr(1)),(hs(1),hsr(1))
c      equivalence (hsder(1),hsdr(1))
c----------------------------------------------------------------------c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      ptest=logkey(ops,'PRINT=M6005=ALL',.false.,' ')
      if (ptest) then
          prntbf=.true.
          prntff=.true.
          prntv=.true.
          prnty=.true.
          prntb=.true.
          prntbs=.true.
          prntgr=.true.
      else
          prntbf=logkey(ops,'PRINT=M6005=BOUND-FREE',.false.,' ')
          prntff=logkey(ops,'PRINT=M6005=FREE-FREE',.false.,' ')
          prntv=logkey(ops,'PRINT=M6005=POTENTIAL',.false.,' ')
          prnty=logkey(ops,'PRINT=M6005=YLM',.false.,' ')
          prntb=logkey(ops,'PRINT=M6005=BOUND',.false.,' ')
          prntbs=logkey(ops,'PRINT=M6005=BESSEL',.false.,' ')
          prntgr=logkey(ops,'PRINT=M6005=GRID',.false.,' ')
      endif
      srtfl='"sorted orbs"'
      typ=.true.
      bcondx='T-matrix'
      filgrd='GRID'
      filorb='ORBS'
      filpot='POT'
      filylm='YLM'
      filbes='BESSFN'
      filint='KNINTS'
      if( posinp('$KOHNINT',cpass) ) then
          call cardin(card)
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
          bcondx=chrkey(card,'boundary-condition','T-matrix',' ')
          call iosys ('does "grid filename" exist on rwf',-1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "grid filename" from rwf',-1,
     1                     0,0,filgrd)
          else
              filgrd=chrkey(card,'grid-file-name','GRID',' ')
          endif
          call iosys ('does "orb file" exist on rwf',-1,0,0,filorb)
          if (ans.ne.'no') then
              call iosys ('read character "orb file" from rwf',-1,0,0,
     1                     filorb)
          else
              filorb=chrkey(card,'numerical-orbital-file-name',
     1                      'ORBS',' ')
          endif
              call iosys ('does "pot file" exist on rwf',-1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "pot file" from rwf',-1,0,0,
     1                    filpot)
          else
              filpot=chrkey(card,'numerical-potential-file-name',
     1                      'VSTAT',' ')
          endif
          call iosys ('does "ylm file" exist on rwf',-1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "ylm file" from rwf',-1,0,0,
     1                    filylm)
          else
              filylm=chrkey(card,'ylm-file-name','YLMS',' ')
          endif
          call iosys ('does "bessel file" exist on rwf',-1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "bessel file" from rwf',-1,
     1                     0,0,filbes)
          else
              filbes=chrkey(card,'bessel-file-name','BESSEL',' ')
          endif
          filint=chrkey(card,'kohn-integral-file-name','KNINTS',' ')
          call captlz(filint,filint)
      endif
      call iosys ('read character "kohn data file" from rwf',0,0,0,
     1             filkne)
      write (iout,1)
      call iosys ('write character "kohn intg file" to rwf',0,0,0,
     1            filint)
      call iosys ('open kohndt as old',0,0,0,filkne)
      call iosys ('open grid as old',0,0,0,filgrd)
      call iosys ('open kohnint as new on ssd',262144,0,0,filint)
      call iosys ('write character "boundary cond" to kohnint',0,0,
     1             0,bcondx)
c----------------------------------------------------------------------c
c            read in spline information for free functions             c
c----------------------------------------------------------------------c
      call iosys ('open bessel as old',0,0,0,filbes)
      call iosys ('read integer maximum-l from bessel',1,lmax,0,' ')
      call iosys ('read integer "total pts" from bessel',1,nr,0,' ')
      call iosys ('read real rmin from bessel',1,rmin,0,' ')
      call iosys ('read real rmax from bessel',1,rmax,0,' ')
      call iosys ('read real "r spacing" from bessel',1,rdel,0,' ')
      rd26= rdel* rdel/6.
c----------------------------------------------------------------------c
c             read grid parameters from basis set file                 c
c----------------------------------------------------------------------c
      call iosys ('open orbs as old',0,0,0,filorb)
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
      call iosys ('open ylms as old',0,0,0,filylm)
      call iosys ('read character "grid type" from ylms',-1,0,0,grdchk)
      call iosys ('read integer "no. grid pts" from ylms',1,ngchk,0,
     1            ' ')
      call iosys ('read integer "point buffer" from ylms',1,pntchk,0,
     1            ' ')
      call iosys ('read integer "max l in ylm" from ylms',1,lmxylm,0,
     1            ' ')
      call iosys ('read integer "max m in ylm" from ylms',1,mumax,0,
     1            ' ')
      if (grdtyp.ne.grdchk) then
          call lnkerr('Mismatch In Type Grids On Orbs And Ylm')
      endif
      if (ngrid.ne.ngchk) then
          call lnkerr('Mismatch In No. Points On Orbs And Ylm')
      endif
      if (pntbuf.ne.pntchk) then
          call lnkerr('Mismatch In Point Buffer On Orbs And Ylm')
      endif
      call iosys ('open vstat as old',0,0,0,filpot)
      call iosys ('read character "grid type" from vstat',-1,0,0,grdchk)
      call iosys ('read integer "no. states" from vstat',1,nchan,0,' ')
      nstri=nchan*(nchan+1)/2
      if (grdchk.ne.grdtyp) then
          call lnkerr('Mismatch In Type Grids On Orbs And Vstat')
      endif
c----------------------------------------------------------------------c
c         read in information about energies and channels              c
c         the scattering energies are incident electron energies       c
c----------------------------------------------------------------------c
      call iosys ('read integer "no. energies" from kohndt',1,
     1             nener,0,' ')
      call iosys ('read real "chan energies" from kohndt',nchan,
     1             echan,0,' ')
      call iosys ('read real "scatt energies" from kohndt',nener,
     1             energy,0,' ')
      if( posinp('$ENERGY',cpass) ) then
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
c             L**2 basis functions to channels                         c
c----------------------------------------------------------------------c
      if( posinp('$TARGET',cpass) ) then
          call cardin(card)
          bsym=logkey(card,'bound-symmetry=off',.false.,' ')
          if (bsym) then
              call iosys ('write character symmetry to kohnint',0,0,
     1                     0,'on')
          else
              call iosys ('write character symmetry to kohnint',0,0,
     1                     0,'off')
          endif
          maxlm=0
          maxbf=0
          do 60 ch1=1,nchan
             if( posinp('$CHAN-'//itoc(ch1),cpass) ) then
                 call cardin(card)
                 nlm(ch1)=intkey(card,'no-lm',1,' ')
                 maxlm=max(maxlm,nlm(ch1))
                 call intarr(card,'l-values',lch(1,ch1),nlm(ch1),' ')
                 call intarr(card,'m-values',mch(1,ch1),nlm(ch1),' ')
                 if (.not.bsym) then
                     ngauss(ch1)=intkey(card,'no-aos',1,' ')
                     maxbf=max(maxbf,ngauss(ch1))
                     call intarr(card,'aos',ngch(1,ch1),ngauss(ch1),' ')
                     nmoc(ch1)=intkey(card,'no-mos',1,' ')
                     call intarr(card,'mos',nmoch(1,ch1),nmoc(ch1),' ')
                 else
                     ngauss(ch1)=ncon
                     maxbf=max(maxbf,ngauss(ch1))
                     nmoc(ch1)=ncon
                     do 65 orb=1,ncon
                        ngch(orb,ch1)=orb
                        nmoch(orb,ch1)=orb
   65                continue
                 endif
                 write(iout,7) ch1
                 write (iout,8) (lch(i,ch1),i=1,nlm(ch1))
                 write (iout,9) (mch(i,ch1),i=1,nlm(ch1))
                 write (iout,10) (ngch(i,ch1),i=1,ngauss(ch1))
                 write (iout,25) (nmoch(i,ch1),i=1,nmoc(ch1))
                 do 70 j=1,nlm(ch1)
                    mm=2*iabs(mch(j,ch1))
                    if (mch(j,ch1).ge.0) then
                        mch(j,ch1)=mm
                    else
                        mch(j,ch1)=mm-1
                    endif
   70            continue
               endif
   60     continue
      endif
      call iosys ('write integer "no. energies" to kohnint',1,
     1             nener,0,' ')
      call iosys ('write real "chan energies" to kohnint',nchan,
     1             echan,0,' ')
      call iosys ('write real "scatt energies" to kohnint',nener,
     1             energy,0,' ')
      call iosys ('write integer "no. channels" to kohnint',1,
     1             nchan,0,' ')
      call iosys ('write real eground to kohnint',1,echan(1),0,' ')
      call iosys ('write integer "lm-vals-chan" to kohnint',nchan,
     1             nlm,0,' ')
      call iosys ('write integer "no.-gauss-chan" to kohnint',nchan,
     1             ngauss,0,' ')
      call iosys ('write integer "no.-mos-chan" to kohnint',nchan,
     1             nmoc,0,' ')
      do 200 ch1=1,nchan
         filnm='"l val ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to kohnint',nlm(ch1),
     1                lch(1,ch1),0,' ')
         filnm='"m val ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to kohnint',nlm(ch1),
     1                mch(1,ch1),0,' ')
         filnm='"aos ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to kohnint',ngauss(ch1),
     1                ngch(1,ch1),0,' ')
         filnm='"mos ch-'//itoc(ch1)//'"'
         call iosys ('write integer '//filnm//' to kohnint',nmoc(ch1),
     1                nmoch(1,ch1),0,' ')
  200 continue
c----------------------------------------------------------------------c
c               Get Memory And Begin Calculation One Energy            c
c                           At A Time                                  c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c   Dimensions Of Bessel And Spherical Functions Are Determined By     c
c   The Largest L And M Values Needed. The Matrix Elements Are         c
c   Dimensioned By Channels And Maximum Number Of L Or Basis Functions c
c                               Needed                                 c
c----------------------------------------------------------------------c
      call iosys ('read maxsiz from rwf',1,maxcor,0,' ')
      call iosys ('read integer "total no. mos" from kohndt',1,nmotot,
     1             0,' ')
      hpvbwd=2*maxlm*nkept*nchan*nchan
      ovbfwd=2*maxlm*nkept*nchan
      ovbfww=2*maxlm*nmotot*nchan
      hpvbww=2*maxlm*nmotot*nchan*nchan
      hpvhpw=2*maxlm*maxlm*nstri
      hpvhmw=2*maxlm*maxlm*nchan*nchan
      beswdr=nr*(lmax+1)
      beswds=2*beswdr
      newwds=2*pntbuf*maxlm*nchan
      vwords=pntbuf*nstri
      ylmwds=pntbuf*(lmax+1)*(mumax+mumax+1)
      baswds=nkept*pntbuf
      words=ncon*nmotot+hpvbwd+ovbfwd+hpvhpw+hpvhmw+
     1      ovbfww+hpvbww+2*maxlm*max(maxlm,ncon)*nchan*nchan+
     2      maxlm*max(maxlm,nkept)+iadtwp(lmax+1)+nr+
     3      2*beswdr+2*beswds+2*newwds+6*pntbuf+vwords+ylmwds+
     4      baswds+4*pntbuf*maxlm+max(4*pntbuf*maxlm,2*pntbuf*nkept)+
     5      maxlm*max(2,pntbuf)
      if (words.gt.maxcor) then
          call lnkerr('Not Enough Memory For Calculation:Quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'M6005',0)
      write (iout,300) words
      vec=ioff
      hpvb=vec+ncon*nmotot
      ovbf=hpvb+hpvbwd
      hpvhp=ovbf+ovbfwd
      hpvhm=hpvhp+hpvhpw
      ovbm=hpvhm+hpvhmw
      hpvm=ovbm+ovbfww
      scrc=hpvm+hpvbww
      scrr=scrc+2*maxlm*max(maxlm,ncon)*nchan*nchan
      ipow=wpadti(scrr+maxlm*max(maxlm,nkept))
      x=iadtwp(ipow+lmax+1)
      hs=x+nr
      hsder=hs+beswds
      cj=hsder+beswdr
      cy=cj+beswds
      hp=cy+beswdr
      newwds=2*pntbuf*maxlm*nchan
      hd=hp+newwds
      rvec=hd+pntbuf*maxlm*nchan
      krvec=rvec+pntbuf
      grid=krvec+pntbuf
      vpot=grid+4*pntbuf
      ylm=vpot+vwords
      basis=ylm+ylmwds
      cmat=basis+baswds
      cmatt=cmat
      rtrick=cmatt
      dmat=cmat+4*pntbuf*maxlm
      rmat=dmat+max(4*pntbuf*maxlm,2*pntbuf*nkept)
      call mkpow(a(ipow),lmax)
      call iosys('read character "Transformation Vector" from kohndt',
     $           -1,0,0,xform)
      call iosys ('read real '//xform//' from kohndt',ncon*nmotot,
     1             z(vec),0,' ')
      call iosys ('read real points from bessel',nr,z(x),0,' ')
      call iosys('read real "complex kohn function" from bessel',
     1            beswds,z(hs),0,' ')            
      call iosys('read real "second derivative of complex kohn '//
     1           'function" from bessel',beswdr,z(hsder),0,' ')            
      call iosys('read real "function spline coefficients" from '//
     1           'bessel',beswds,z(cj),0,' ')            
      call iosys('read real "second derivative spline coefficients" '//
     1            'from  bessel',beswdr,z(cy),0,' ')            
c----------------------------------------------------------------------c
c                      Loop Over Energies                              c
c----------------------------------------------------------------------c
      do 400 ene=1,nener
c----------------------------------------------------------------------c
c                 construct channel momenta                            c
c----------------------------------------------------------------------c
         do 410 ch1=1,nchan
            ec = .5d+00*energy(ene) - (echan(ch1)-echan(1))
            if(ec.le.0.0) then
               call lnkerr('No Closed Channels Allowed At Present:Quit')
            endif
            kchan(ch1) = sqrt(2.0*ec)
  410    continue
         write(iout,420) energy(ene), (kchan(i),i=1,nchan)
c----------------------------------------------------------------------c
c              Rewind All Energy-Independent Files                     c
c----------------------------------------------------------------------c
         call iosys ('rewind all on grid read-and-write',0,0,0,' ')
         call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
         call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
         call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
c----------------------------------------------------------------------c
c                   Zero Matrix Elements                               c
c----------------------------------------------------------------------c
         call rzero(z(hpvhp),hpvhpw)
         call rzero(z(hpvhm),hpvhmw)
         call rzero(z(hpvb),hpvbwd)
         call rzero(z(ovbf),ovbfwd)
c----------------------------------------------------------------------c
c                   Loop Over Grid                                     c
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
c                Compute The Bessel Functions From The Spline          c
c----------------------------------------------------------------------c
            call mkbes (z(hs),z(hsder),z(cj),z(cy),z(hp),z(hd),z(grid),
     1                  z(rvec),z(x),z(krvec),kchan,rmin,rdel,rd26,
     2                  a(ipow),nlm,lch,npnts,nr,lmax,maxlm,nchan,
     3                  dimlm,dimc)
c----------------------------------------------------------------------c
c           read static potentials for this block                      c
c           weights are already multiplied into potentials             c
c----------------------------------------------------------------------c
            call iosys ('read real "static potential" from vstat '//
     1                  'without rewinding',npnts*nstri,z(vpot),0,' ')
            if (prntgr) then
                call grdprn(z(grid),npnts,ireg)
            endif
            if (prntv) then
                call vprnt(z(vpot),npnts,nstri,ireg)
            endif
            call scalev(z(vpot),z(grid),npnts,nstri)
            if (prnty) then
                call yprnt(z(ylm),npnts,lmax,mumax,ireg)
            endif
            if (prntb) then
                call bndprn(z(basis),npnts,nkept,ireg)
            endif
            if (prntbs) then
                call bsprnt(z(hp),z(hd),nchan,npnts,nlm,lch,maxlm,
     1                      dimlm,dimc,ireg)
            endif
c----------------------------------------------------------------------c
c     At this point a block of everything on the grid has been         c
c     read in.  Now call routines which accumulate the bound-free      c
c     and free-free matrix elements from contributions at              c
c                        each grid point                               c
c----------------------------------------------------------------------c
            call ffints(z(hpvhp),z(hpvhm),z(grid),z(hp),z(hd),z(ylm),
     1                  z(vpot),z(cmatt),z(cmat),z(dmat),z(rmat),
     2                  z(rtrick),z(scrc),z(scrr),nlm,lch,mch,npnts,
     3                  nchan,lmax,nstri,dimlm,dimc,maxlm,bcondx)
            call bfints(z(hpvb),z(ovbf),z(grid),z(basis),z(hp),z(hd),
     1                  z(ylm),z(vpot),z(cmat),z(rmat),z(dmat),z(scrc),
     2                  z(scrr),nlm,lch,mch,ngauss,ngch,npnts,nchan,
     3                  nkept,list,lmax,nstri,dimlm,dimc,dimbf,maxlm)
  430    continue
         filnm='ch-en-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',nchan,kchan,
     1               0,' ')
c-----------------------------------------------------------------------c
c              Transform To Molecular Orbital Basis                     c
c-----------------------------------------------------------------------c
         call tomobs(z(vec),z(ovbf),z(ovbm),z(hpvb),z(hpvm),z(scrc),
     1               maxlm,ncon,nmotot,nchan,ngauss,ngch,nkept,
     2               list,dimbf,dimc)
         filnm='ovlp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',ovbfww,z(ovbm),
     1               0,' ')
         filnm='bf-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',hpvbww,z(hpvm),
     1               0,' ')
         filnm='ffp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',hpvhpw,
     1               z(hpvhp),0,' ')
         filnm='ffm-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',hpvhmw,
     1               z(hpvhm),0,' ')
c----------------------------------------------------------------------c
c                    output section                                    c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c      write bound-free matrices to printer if requested               c
c----------------------------------------------------------------------c
         if( prntbf) then
             write (iout,1000)
             call bfprnt(z(ovbf),z(hpvb),nlm,nchan,dimc,maxlm,nkept)
             write (iout,2000)
             call bfprnt(z(ovbm),z(hpvm),nlm,nchan,dimc,maxlm,nmotot)
         endif
         if (prntff) then
             call ffprnt(z(hpvhp),z(hpvhm),nlm,nchan,nstri,dimc,maxlm)
         endif
 
c close big loop on incident energies
c
  400 continue
c----------------------------------------------------------------------c
c                   Rewind And Close All Files                         c
c----------------------------------------------------------------------c
      call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys ('close kohndt',0,0,0,' ')
      call iosys ('rewind all on grid read-and-write',0,0,0,' ')
      call iosys ('close grid',0,0,0,' ')
      call iosys ('rewind all on kohnint read-and-write',0,0,0,' ')
      call iosys ('close kohnint',0,0,0,' ')
      call iosys ('rewind all on bessel read-and-write',0,0,0,' ')
      call iosys ('close bessel',0,0,0,' ')
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
      call iosys ('close ylms',0,0,0,' ')
      call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
      call iosys ('close vstat',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
    1 format(//,20x,'*****  M6005:Bound-Free And Free-Free Integrals ***
     1**')
    3 format(/,5x,'Target Energies For ',i3,' Channels:',/,(2x,5e15.8))
    4 format(/,5x,'Incident Energies:',/,(2x,5e15.8))
    5 format(/,5x,'Spline Data',/,5x,'Rmin',1x,f10.5,1x,'Rmax',f10.5,
     11x,'Step Size',1x,f10.5,/,5x,'Max. L',1x,i3,1x,'Inner Cutoff',1x,
     2 f6.3,1x,'Outer Exponential Cutoff',1x,f6.3,1x,'Outer N Cutoff',
     3 1x,i3)
    6 format(/,5x,'Boundary Conditions Are',1x,a10)
    7 format(/,5x,'Information For Channel',1x,i3)
    8 format(/,5x,'L-Values',(/,15x,10(i2,1x)))
    9 format(/,5x,'M-Values',(/,15x,10(i2,1x)))
   10 format(/,5x,'Bound-Atomic-Orbitals',(/,17x,10(i3,1x)))
   25 format(/,5x,'Bound-Molecular-Orbitals',(/,17x,10(i3,1x)))
  300 format(/,5x,'Need',1x,i8,1x,'Words For Calculation')
  420 format(/,5x,'Incident Energy(Ry)',1x,f10.6,/,5x,'Channel Momenta',
     1       (/,10x,5(e12.5,1x)))
 1000 format(/,5x,'Atomic Orbital Representation')
 1500 format(/,5x, 'Will Use Short Orbital List')
 2000 format(/,5x,'Molecular Orbital Representation')
      call chainx(0)
      stop
      end
*Deck Mkbes
c***begin prologue     Mkbes
c***date written       XXXXXX   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6005, link 6005, Bessel Functions
c***author             Rescigno, Tom (LLNL)
c***source             m6005
c***purpose            Bessel Functions From Spline Fit
c***description        Calculation Of Bessel Functions From Spline
c***                   Fits In Link M6004
c***routines called    IOsys, Util and MDutil
c***end prologue       Mkbes
      subroutine mkbes(hs,hsder,cj,cy,hp,hd,grid,rvec,x,krvec,kchan,
     1                         rmin,rdel,rd26,ipow,nlm,lch,npt,nr,
     2                         lmax,maxlm,nchan,dimlm,dimc)
      implicit integer (a-z)
      real *8 grid, rvec, x, krvec, kchan, rmin, rdel, rd26, ksq, k52
      real *8 a, b, hsder, cy, hd
      complex *16 hs, cj, hp
      dimension grid(4,npt), rvec(npt), krvec(npt), kchan(nchan)
      dimension hs(nr,0:lmax), hsder(nr,0:lmax), cj(nr,0:lmax)
      dimension cy(nr,0:lmax), hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension ipow(0:lmax), x(nr), nlm(dimc), lch(dimlm,dimc)
      do 10 i = 1,npt
         rvec(i)=sqrt(grid(1,i)**2 + grid(2,i)**2 +grid(3,i)**2)
   10 continue
      do 20 ch1=1,nchan
         do 30 i=1,npt
            krvec(i)=kchan(ch1)*rvec(i)
   30    continue
         ksq = sqrt(kchan(ch1))
         k52 = ksq*kchan(ch1)**2
         do 40 lv=1,nlm(ch1)
               l=lch(lv,ch1)
            do 50 i=1,npt
c********************
c the expression for klo is split into two fortran
c statements so that klo is never less than 1
c*****************
               klo=(krvec(i)-rmin)/rdel
               klo=klo+1
               a=(x(klo+1)-krvec(i))/rdel
               b=(krvec(i)-x(klo))/rdel
               hp(i,lv,ch1)=(a*hs(klo,l)+b*hs(klo+1,l)+(a*(a*a-1.)*
     1                    cj(klo,l)+b*(b*b-1.)*cj(klo+1,l))*rd26)*ksq
               hd(i,lv,ch1)=a*hsder(klo,l)+b*hsder(klo+1,l)+(a*(a*a-1.)*
     1                    cy(klo,l)+b*(b*b-1.)*cy(klo+1,l))*rd26
               hd(i,lv,ch1)=hd(i,lv,ch1)*k52/krvec(i)**(1-ipow(l))
   50       continue
   40    continue
   20 continue
      return
      end
*Deck FFints
c***begin prologue     FFints
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           Kohn Integrals
c***author             Schneider, Barry (LANL)
c***source             m6005
c***purpose            Free-Free Matrix Elements
c***description        Calculation Of Free-Free Integrals
c***references
c
c***routines called    IOsys, Util and MDutil
c***end prologue       FFints
      subroutine ffints(hpvhp,hpvhm,grid,hp,hd,ylm,vpot,cmat,cmatt,dmat,
     1                  rmat,rtrick,scrc,scrr,nlm,lch,mch,npt,nchan,
     2                  lmax,nstri,dimlm,dimc,maxlm,bcond)
      implicit integer (a-z)
      character *(*) bcond
      real *8 ylm, vpot, grid, scrr, hd, sdot, rmat, rtrick
      complex *16 cmat, cmatt, dmat, hpvhp, scrc
      complex *16 hpvhm, hp, ai
      dimension ylm(npt,0:lmax,0:2*lmax), vpot(npt,nstri)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt,2), dmat(maxlm,npt,2), grid(4,npt)
      dimension hpvhp(1:maxlm,1:maxlm,nstri), scrc(*), scrr(*)
      dimension hpvhm(1:maxlm,1:maxlm,nchan,nchan), cmatt(npt,maxlm)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension rmat(maxlm,2), rtrick(2*npt,maxlm)
      ai=cmplx(0.d+00,1.d+00)
c
c compute potential matrices
c
c----------------------------------------------------------------------c
c         First Two Loops Over Channel 1 And Associated L, M           c
c         To Load Some Temporary Matrices For Vectorization            c
      do 10 ch1=1,nchan
         call czero(cmat(1,1,1),maxlm*npt)
         call czero(dmat(1,1,1),maxlm*npt)
         call filcf1(cmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm)
         call fildf1(dmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm,bcond)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c          Next Two Loops Over Channel 2 And Associated L, M           c
c          Same Strategy On Temporary Matrices                         c
         do 20 ch2=1,ch1
            call czero(cmat(1,1,2),maxlm*npt)
            call czero(dmat(1,1,2),maxlm*npt)
            ist = ch1*(ch1-1)/2 + ch2
            call filcf2(cmat(1,1,2),ylm,hp(1,1,ch2),vpot(1,ist),
     1                  nlm(ch2),lch(1,ch2),mch(1,ch2),maxlm,npt,
     2                  lmax,dimlm)
            call fildf2(dmat(1,1,2),ylm,hp(1,1,ch2),vpot(1,ist),
     1                  nlm(ch2),lch(1,ch2),mch(1,ch2),maxlm,npt,
     2                  lmax,dimlm,bcond)
c----------------------------------------------------------------------c
c                  Accumulate Into Matrix Elements                     c
c----------------------------------------------------------------------c
            call vcebct(scrc,cmat(1,1,1),cmat(1,1,2),scrr,nlm(ch1),
     1                  npt,nlm(ch2))
c           call cebc(scrc,cmat(1,1,1),cmat(1,1,2),nlm(ch1),npt,
c    1                nlm(ch2))
            call cvadd(hpvhp(1,1,ist),maxlm,scrc,nlm(ch1),nlm(ch1),
     1                 nlm(ch2))
c           call capbctx(hpvhp(1,1,ist),maxlm,cmat(1,1,1),nlm(ch1),
c    1                   cmat(1,1,2),nlm(ch2),nlm(ch1),npt,nlm(ch2))
            if (ch1.ne.ch2) then
c               call cebc(scrc,cmat(1,1,1),dmat(1,1,2),nlm(ch1),npt,
c    1                    nlm(ch2))
                call vcebct(scrc,cmat(1,1,1),dmat(1,1,2),scrr,nlm(ch1),
     1                      npt,nlm(ch2))
                call cvadd(hpvhm(1,1,ch1,ch2),maxlm,scrc,nlm(ch1),
     1                     nlm(ch1),nlm(ch2))
c               call capbctx(hpvhm(1,1,ch1,ch2),maxlm,cmat(1,1,1),
c    1                       nlm(ch1),dmat(1,1,2),nlm(ch2),nlm(ch1),
c    2                       npt,nlm(ch2))
            endif
c           call cebc(scrc,cmat(1,1,2),dmat(1,1,1),nlm(ch2),npt,
c    1                nlm(ch1))
            call vcebct(scrc,cmat(1,1,2),dmat(1,1,1),scrr,nlm(ch2),npt,
     1                  nlm(ch1))
            call cvadd(hpvhm(1,1,ch2,ch1),maxlm,scrc,nlm(ch2),nlm(ch2),
     1                 nlm(ch1))
c           call capbctx(hpvhm(1,1,ch2,ch1),maxlm,cmat(1,1,2),
c    1                   nlm(ch2),dmat(1,1,1),nlm(ch1),nlm(ch2),
c    2                   npt,nlm(ch1))
   20    continue
   10 continue
c----------------------------------------------------------------------c
c            compute matrix of (kchan**2/2 - T)                        c
c            diagonal in channel indices so its simplier than above    c
c       rtrick (real) and cmatt (complex) must be equivalenced         c
c       in the call to this routine for it to function properly        c
c----------------------------------------------------------------------c
      do 30 ch1=1,nchan
         clst=ch1*(ch1+1)/2
         call czero(cmatt,maxlm*npt)
         do 40 nolm1=1,nlm(ch1)
            l1=lch(nolm1,ch1)
            m1=mch(nolm1,ch1)
            do 50 grpt=1,npt
               cmatt(grpt,nolm1) = grid(4,grpt)*ylm(grpt,l1,m1)
     1                             *ylm(grpt,l1,m1)
     2                             *hp(grpt,nolm1,ch1)
   50       continue
   40    continue
         do 60 nolm1=1,nlm(ch1)
            rmat(nolm1,1)=sdot(npt,rtrick(1,nolm1),2,
     1                                hd(1,nolm1,ch1),1)
            rmat(nolm1,2)=sdot(npt,rtrick(2,nolm1),2,
     1                             hd(1,nolm1,ch1),1)
            hpvhp(nolm1,nolm1,clst)=hpvhp(nolm1,nolm1,clst) +
     1                               rmat(nolm1,1) + ai*rmat(nolm1,2)
   60    continue
         if(bcond.eq.'S-matrix') then
            do 70 nolm1=1,nlm(ch1)
               hpvhm(nolm1,nolm1,ch1,ch1)=hpvhm(nolm1,nolm1,ch1,ch1) +
     1                                    rmat(nolm1,1) +
     2                                    ai*rmat(nolm1,2)
   70       continue
         endif
   30 continue
      return
      end
*Deck BFints
c***begin prologue     BFints
c***date written       XXXXXX   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           Kohn Integrals
c***author             Schneider, Barry (LANL)
c***source             m6005
c***purpose            Bound-Free Matrix Elements
c***description        Calculation Of Bound-Free Integrals
c***references
c
c***routines called    IOsys, Util and MDutil
c***end prologue       BFints
      subroutine bfints(hpvb,ovbf,grid,basis,hp,hd,ylm,vpot,cmat,rmat,
     1                  bmat,scrc,scrr,nlm,lch,mch,ngauss,ngch,npt,
     2                  nchan,nkept,list,lmax,nstri,dimlm,dimc,dimbf,
     3                  maxlm)
      implicit integer (a-z)
      real *8 hd, ylm, vpot, grid, basis, scrr, rmat
      complex *16  cmat, bmat, hpvb, scrc
      complex *16 ovbf, hp
      dimension ylm(npt,0:lmax,0:2*lmax), vpot(npt,nstri)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt), bmat(npt,nkept), grid(4,npt)
      dimension rmat(maxlm,npt)
      dimension hpvb(1:maxlm,nkept,nchan,nchan), scrc(*)
      dimension ovbf(1:maxlm,nkept,nchan), scrr(*)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension basis(npt,nkept), list(dimbf)
      dimension ngauss(dimc), ngch(dimbf,dimc)
c----------------------------------------------------------------------c
c          First Two Loops Over Channel 1 And Associated L, M          c
 
      do 10 ch1=1,nchan
         call czero(cmat,maxlm*npt)
         call filcf1(cmat,ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             Next Two Loops Over Channel 2 And Associated             c
c                        Bound Orbitals                                c
         do 20 ch2=1,nchan
            call czero(bmat,npt*nkept)
            chlrg=max(ch1,ch2)
            chsml=ch1+ch2-chlrg
            itri=chlrg*(chlrg-1)/2+chsml
            do 30 bfn=1,ngauss(ch2)
               bfnch2 = ngch(bfn,ch2)
               blocte=list(bfnch2)
               if (blocte.eq.0) then
                   call lnkerr('Bound Orbital Out Of Range')
               endif
c----------------------------------------------------------------------c
               do 40 grpt=1,npt
                  bmat(grpt,blocte)=vpot(grpt,itri)*basis(grpt,blocte)
   40          continue
   30       continue
            call vcebc(scrc,cmat,bmat,scrr,nlm(ch1),npt,nkept)
            call cvadd(hpvb(1,1,ch1,ch2),maxlm,scrc,nlm(ch1),nlm(ch1),
     1                 nkept)
   20    continue
   10 continue
c----------------------------------------------------------------------c
c            compute matrix of (E - T)                                 c
c----------------------------------------------------------------------c
      do 70 ch1=1,nchan
         call czero(cmat,maxlm*npt)
         call rzero(rmat,maxlm*npt)
         call filc12(cmat,rmat,ylm,grid,hp(1,1,ch1),hd(1,1,ch1),
     1               nlm(ch1),lch(1,ch1),mch(1,ch1),maxlm,npt,lmax,
     2               dimlm)
         call czero(bmat,npt*nkept)
         do 80 bfn=1,ngauss(ch1)
            bfnch1=ngch(bfn,ch1)
            blocte=list(bfnch1)
            if (blocte.eq.0) then
                call lnkerr('Bound Orbital Out Of Range')
            endif
            do 90 grpt=1,npt
               bmat(grpt,blocte)=basis(grpt,blocte)
   90       continue
   80    continue
         call vcebc(scrc,cmat,bmat,scrr,nlm(ch1),npt,nkept)
         call cvadd(ovbf(1,1,ch1),maxlm,scrc,nlm(ch1),nlm(ch1),nkept)
         call ebc(scrr,rmat,basis,nlm(ch1),npt,nkept)
         call crvadd(hpvb(1,1,ch1,ch1),maxlm,scrr,nlm(ch1),nlm(ch1),
     1              nkept)
   70 continue
      return
      end
*Deck CRvadd
c***Begin Prologue     CRvadd
c***Date Written       890602   (yymmdd)
c***Revision Date               (yymmdd)
c***Keywords           Complex Matrix Add
c***Author             Schneider, Barry (Lanl)
c***Source             Mylib
c***Purpose            Add Complex And Real Matrix
c***          
c***Description        
c***                   
c***                   
c
c***References       
c
c***Routines called   
c***End Prologue      
      subroutine crvadd(a,n1,b,n2,n3,n4)
      implicit integer (a-z)
      real *8 b
      complex*16 a
      dimension a(n1,n4), b(n2,n4)
      do 10 i=1,n3
         do 20 j=1,n4
            a(i,j)=a(i,j)+b(i,j)
   20    continue
   10 continue
      return
      end
*Deck Tomobs
c***begin prologue     Tomobs
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           Transformation, Kohn
c***author             Schneider, Barry (LANL)
c***source             M6005
c***purpose            Bound-Free Ao To Mo Transformation
c***
c***references         None
c
c***routines called    Cebc
c***end prologue       Tomobs
 
      subroutine tomobs(trans,ovbf,ovbm,hpvb,hpvm,scrc,maxlm,ncon,nmo,
     1                  nchan,ngauss,ngch,nkept,list,dimbf,dimc)
      implicit integer (a-z)
      real *8 trans
      complex *16 ovbf, ovbm, hpvb, hpvm, scrc
      dimension ovbf(1:maxlm,nkept,nchan), ovbm(1:maxlm,nmo,nchan)
      dimension hpvb(1:maxlm,nkept,nchan,nchan), trans(ncon,nmo)
      dimension hpvm(1:maxlm,nmo,nchan,nchan), list(dimbf)
      dimension scrc(1:maxlm,ncon,nchan,nchan)
      dimension ngauss(dimc), ngch(dimbf,dimc)
      nzero=maxlm*ncon*nchan*nchan
      call czero(scrc,nzero)
      do 10 ch1=1,nchan
         do 20 bfn=1,ngauss(ch1)
            bfnch = ngch(bfn,ch1)
            blocte=list(bfnch)
            do 30 lm=1,maxlm
               scrc(lm,bfnch,ch1,1)=ovbf(lm,blocte,ch1)
   30       continue
   20    continue
   10 continue
      do 40 ch1=1,nchan
         call ecbcx(ovbm(1,1,ch1),maxlm,scrc(1,1,ch1,1),maxlm,trans,
     1             ncon,maxlm,ncon,nmo)
   40 continue
      call czero(scrc,nzero)
      do 50 ch1=1,nchan
         do 60 ch2=1,nchan
            do 70 bfn=1,ngauss(ch2)
               bfnch = ngch(bfn,ch2)
               blocte=list(bfnch)
               do 80 lm=1,maxlm
                  scrc(lm,bfnch,ch1,ch2)=hpvb(lm,blocte,ch1,ch2)
   80          continue
   70       continue
   60    continue
   50 continue
      do 90 ch1=1,nchan
         do 100 ch2=1,nchan
            call ecbcx(hpvm(1,1,ch1,ch2),maxlm,scrc(1,1,ch1,ch2),maxlm,
     1                 trans,ncon,maxlm,ncon,nmo)
  100    continue
   90 continue
      return
      end
*Deck Mkpow
c***begin prologue     Mkpow
c***date written       XXXXXX   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           Kohn, Integrals
c***author             Rescigno, Tom (LLNL)
c***source             m6005
c***purpose            Fill Array Ipow
c***references
c
c***routines called    None
c***end prologue       Mkpow
      subroutine mkpow(ipow,lmax)
      implicit integer (a-z)
      dimension ipow(0:lmax)
      data zero, one/ 0, 1/
      ipow(0)=zero
      do 10 i=1,lmax
         ipow(i)=one
   10 continue
      return
      end
*Deck BFprnt
c***begin prologue     BFprnt
c***date written       XXXXXX   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           Kohn Integrals
c***author             Schneider, Barry (LANL)
c***source             m6005
c***purpose            Print Bound-Free Matrix Elements
c***references
c
c***routines called
      subroutine bfprnt(ovbf,hpvb,nlm,nchan,dimc,maxlm,ncon)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 rowv
      character *80 title
      character *3 itoc
      character *8 rowt, colt
      complex *16 ovbf, hpvb
      dimension ovbf(maxlm,ncon,nchan), hpvb(maxlm,ncon,nchan,nchan)
      dimension nlm(dimc)
      colv=-99
      rowv=-99.
      do 10 ch1=1,nchan
         do 20 ch2=1,nchan
            title='HPvB Matrix Ch1-'//itoc(ch1)//' Ch2-'//itoc(ch2)
            call cmprir(hpvb(1,1,ch1,ch2),rowv,colv,nlm(ch1),
     1                  ncon,maxlm,ncon,title,rowt,colt,iout)
   20    continue
   10 continue
      do 50 ch1=1,nchan
         title='Bound Free Overlap Matrix Ch1-'//itoc(ch1)
         call cmprir(ovbf(1,1,ch1),rowv,colv,nlm(ch1),ncon,
     1               maxlm,ncon,title,rowt,colt,iout)
   50 continue
      return
      end
*Deck FFprnt
c***begin prologue     FFprnt
c***date written       XXXXXX   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           Kohn Integrals
c***author             Schneider, Barry (LANL)
c***source             m6005
c***purpose            Print Free-Free Matrix Elements
c***references
c
c***routines called
      subroutine ffprnt(hpvhp,hpvhm,nlm,nchan,nstri,dimc,maxlm)
      implicit integer (a-z)
      character *80 title
      character *8 colt, rowt
      character *3 itoc
      complex *16 hpvhp, hpvhm
      real *8 rowv
      common /io/ inp, iout
      dimension hpvhp(maxlm,maxlm,nstri)
      dimension hpvhm(maxlm,maxlm,nchan,nchan)
      dimension nlm(dimc)
      ist=0
      do 10 ch1=1,nchan
         do 20 ch2=1,ch1
            ist=ist+1
            title='HPvHP Matrix Ch1-'//itoc(ch1)//' Ch2-'//itoc(ch2)
            rowt='LM Index'
            colt=rowt
            rowv=-99.
            colv=-99
            call cmprir(hpvhp(1,1,ist),rowv,colv,nlm(ch1),nlm(ch2),
     1               maxlm,maxlm,title,rowt,colt,iout)
   20    continue
   10 continue
      do 70 ch1=1,nchan
         do 80 ch2=1,nchan
            title='HPvHM Matrix Ch1-'//itoc(ch1)//' Ch2-'//itoc(ch2)
            rowt='LM Index'
            colt=rowt
            rowv=-99.
            colv=-99
            call cmprir(hpvhm(1,1,ch1,ch2),rowv,colv,nlm(ch1),nlm(ch2),
     1               maxlm,maxlm,title,rowt,colt,iout)
   80    continue
   70 continue
      return
      end
*Deck Grdprn
c***begin prologue     Grdprn
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, Print Grid
c***author             Schneider, Barry (lanl)
c***source             m6005
c***purpose            Print Co-Ordinate Grid
c***references         None
c
c***routines called    None
c***end prologue       Grdprn
      subroutine grdprn(grid,npnts,reg)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 grid
      dimension grid(4,npnts)
      write (iout,10) reg
      write (iout,20)
      do 30 grpt=1,npnts
         write (iout,40) grid(1,grpt), grid(2,grpt), grid(3,grpt)
   30 continue
   40 format(5x,3f10.6)
      return
   10 format (/,5x,'Grid Points For Region',1x,i4)
   20 format(/,5x,'     X     ','     Y     ','     Z     ')
      end
*Deck Vprnt
c***begin prologue     Vprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, Print Potential
c***author             Schneider, Barry (lanl)
c***source             m6005
c***purpose            Print Static Potentials On Grid
c***references         None
c
c***routines called    None
c***end prologue       Vprnt
      subroutine vprnt(v,npnts,ntri,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character *8 colt, rowt
      common /io/ inp, iout
      real *8 v, rowv
      dimension v(npnts,ntri)
      title='Potential Matrices For Region-'//itoc(reg)
      rowv=-99.
      colv=-99
      colt='State'
      rowt='Point'
      call mprir(v,rowv,colv,npnts,ntri,npnts,ntri,title,rowt,colt,
     1            iout)
      return
      end
*Deck Yprnt
c***begin prologue     Yprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, Spherical Harmonic Print
c***author             Schneider, Barry (lanl)
c***source             m6005
c***purpose            Print Spherical Harmonics
c***references         None
c
c***routines called    None
c***end prologue       Yprnt
      subroutine yprnt(ylm,npnts,lmax,mumax,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character *8 colt, rowt
      common /io/ inp, iout
      real *8 ylm, rowv
      dimension ylm(npnts,0:lmax,0:2*mumax)
      dimension lval(0:100)
      do 1 l=0,lmax
         lval(l)=l
    1 continue
      rowt='Point'
      colt='L'
      rowv=-99.
      write (iout,10) reg
      mcnt=0
      do 20 m=0,mumax
         iup=2
         if (m.eq.0) then
             iup=1
         endif
         do 40 mm=1,iup
            title='Mu-'//itoc(m)//' Comp-'//itoc(mm)
            call mprir(ylm(1,m,mcnt),rowv,lval(m),npnts,lmax+1-m,npnts,
     1                 lmax+1,title,rowt,colt,iout)
            mcnt=mcnt+1
   40    continue
   20 continue
   10 format(/,5x,'Ylm For Region',1x,i4)
      return
      end
*Deck Bndprn
c***begin prologue     Bndprn
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, Print
c***author             Schneider, Barry (lanl)
c***source             m6005
c***purpose            Print Bound Orbitals
c***references         None
c
c***routines called    None
c***end prologue       Bndprn
      subroutine bndprn(orbs,npnts,ncon,reg)
      implicit integer (a-z)
      common /io/ inp, iout
      character *80 title
      character *8 rowt, colt
      character *3 itoc
      real *8 orbs, rowv
      dimension orbs(npnts,ncon)
      title='Bound Orbitals Region-'//itoc(reg)
      rowv=-99.
      colv=-99
      colt='Con Ao'
      rowt='Point'
      call mprir(orbs,rowv,colv,npnts,ncon,npnts,ncon,title,rowt,colt,
     1            iout)
      return
      end
*Deck BSprnt
c***begin prologue     BSprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, Bessel Print
c***author             Schneider, Barry (lanl)
c***source             m6005
c***purpose            Print Channel "Bessel" Functions
c***routines called    None
c***end prologue       Bsprnt
      subroutine bsprnt(hp,hd,nchan,npnts,nlm,lch,maxlm,dimlm,
     1                  dimc,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character*8 rowt, colt
      common /io/ inp, iout
      complex *16 hp, hd
      real *8 rowv
      dimension hp(npnts,maxlm,nchan), hd(npnts,maxlm,nchan)
      dimension nlm(dimc), lch(dimlm,dimc)
      rowv=-99.
      rowt='Point'
      colt='L'
      write (iout,10) reg
      do 20 ch1=1,nchan
         title='Hp Functions Channel-'//itoc(ch1)
         call cmprir(hp(1,1,ch1),rowv,lch(1,ch1),npnts,nlm(ch1),npnts,
     1               maxlm,title,rowt,colt,iout)
         title='Hd Functions Channel-'//itoc(ch1)
         call cmprir(hd(1,1,ch1),rowv,lch(1,ch1),npnts,nlm(ch1),npnts,
     1               maxlm,title,rowt,colt,iout)
   20 continue
      return
   10 format(/,5x,'Free Functions For Region',1x,i4)
      end
*Deck Scalev
c***begin prologue     Scalev
c***date written       890512   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             Schneider, Barry (Lanl)
c***source             M6005
c***purpose            Multiply Potential By Integration Weights
c
c***references
c
c***routines called   None
c***end prologue      Scalev
      subroutine scalev(v,grid,npnts,nstri)
      implicit integer (a-z)
      real *8 grid, v
      dimension v(npnts,nstri), grid(4,npnts)
      do 10 i=1,nstri
         do 20 grpt=1,npnts
            v(grpt,i)=v(grpt,i)*grid(4,grpt)
   20    continue
   10 continue
      return
      end
*Deck Filcf1
c***begin prologue     Filcf1
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           Utility, Matrix
c***author             Schneider, Barry (lanl)
c***source             M6005
c***purpose            Fill Matrix In Packed Form
c***description
c***
c***
c***references         None
c
c***routines called
c***end prologue       Filcf1
      subroutine filcf1(c,a,b,nc,l,m,maxl,npt,lmax,dim)
      implicit integer (a-z)
      complex *16 b, c
      real *8 a
      dimension c(nc,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim)
      do 10 nolm1=1,nc
         l1=l(nolm1)
         m1=m(nolm1)
         do 20 grpt=1,npt
            c(nolm1,grpt)=a(grpt,l1,m1)*b(grpt,nolm1)
   20    continue
   10 continue
      return
      end
*Deck Filcf2
c***begin prologue     Filcf2
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           Utility, Matrix
c***author             Schneider, Barry (lanl)
c***source             M6005
c***purpose            Fill Matrix In Packed Form
c***description
c***
c***
c***references         None
c
c***routines called
c***end prologue       Filcf2
      subroutine filcf2(c,a,b,v,nc,l,m,maxl,npt,lmax,dim)
      implicit integer (a-z)
      complex *16 b, c
      real *8 a, v
      dimension c(nc,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim), v(npt)
      do 10 nolm1=1,nc
         l1=l(nolm1)
         m1=m(nolm1)
         do 20 grpt=1,npt
            c(nolm1,grpt)=a(grpt,l1,m1)*v(grpt)*b(grpt,nolm1)
   20    continue
   10 continue
      return
      end
*Deck Fildf1
c***begin prologue     Fildf1
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           Utility, Matrix
c***author             Schneider, Barry (lanl)
c***source             M6005
c***purpose            Fill Matrix In Packed Form
c***description
c***
c***
c***references         None
c
c***routines called
c***end prologue       Fildf1
      subroutine fildf1(d,a,b,nd,l,m,maxl,npt,lmax,dim,bcond)
      implicit integer (a-z)
      complex *16  b, d
      real *8 a
      character *(*) bcond
      dimension d(nd,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim)
c----------------------------------------------------------------------c
c                 Test Boundary Conditions                             c
c----------------------------------------------------------------------c
      if (bcond.eq.'S-matrix') then
          do 10 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 20 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*conjg(b(grpt,nolm1))
   20        continue
   10     continue
      elseif (bcond.eq.'T-matrix') then
           do 30 nolm1=1,nd
              l1=l(nolm1)
              m1=m(nolm1)
              do 40 grpt=1,npt
                 d(nolm1,grpt)=a(grpt,l1,m1)*aimag(b(grpt,nolm1))
   40        continue
   30     continue
      else
          call lnkerr('Boundary Condition Error In FFints')
      endif
      return
      end
*Deck Fildf2
c***begin prologue     Fildf2
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           Utility, Matrix
c***author             Schneider, Barry (lanl)
c***source             M6005
c***purpose            Fill Matrix In Packed Form
c***description
c***
c***
c***references         None
c
c***routines called
c***end prologue       Fildf2
      subroutine fildf2(d,a,b,v,nd,l,m,maxl,npt,lmax,dim,bcond)
      implicit integer (a-z)
      complex *16 b, d
      real *8 a, v
      character *(*) bcond
      dimension d(nd,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim), v(npt)
c----------------------------------------------------------------------c
c                 Test Boundary Conditions                             c
c----------------------------------------------------------------------c
      if (bcond.eq.'S-matrix') then
          do 10 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 20 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*v(grpt)*conjg(b(grpt,nolm1))
   20        continue
   10     continue
      elseif (bcond.eq.'T-matrix') then
          do 30 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 40 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*v(grpt)*aimag(b(grpt,nolm1))
   40       continue
   30    continue
      else
         call lnkerr('Boundary Condition Error In FFints')
      endif
      return
      end
*Deck Filc12
c***begin prologue     Filc12
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           Utility, Matrix
c***author             Schneider, Barry (lanl)
c***source             M6005
c***purpose            Fill Matrix In Packed Form
c***description
c***
c***
c***references         None
c
c***routines called
c***end prologue       Filc12
      subroutine filc12(cc,cr,a,b,h1,h2,nc,l,m,maxl,npt,lmax,dim)
      implicit integer (a-z)
      complex*16 cc, h1
      real *8 a, b, h2, cr
      dimension cc(nc,npt), a(npt,0:lmax,0:2*lmax), b(4,npt)
      dimension h1(npt,maxl), h2(npt,maxl), l(dim), m(dim)
      dimension cr(nc,npt)
      do 10 nolm1=1,nc
         l1=l(nolm1)
         m1=m(nolm1)
         do 20 grpt=1,npt
            cc(nolm1,grpt)=a(grpt,l1,m1)*b(4,grpt)*h1(grpt,nolm1)
            cr(nolm1,grpt)=a(grpt,l1,m1)*b(4,grpt)*h2(grpt,nolm1)
   20    continue
   10 continue
      return
      end
