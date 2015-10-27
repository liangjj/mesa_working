*deck @(#)ffbf.f	1.1 9/8/91
c***begin prologue     m6005
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6005, link 6005, orbital decomposition
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            free-free and bound-free matrix elements
c***description        calculation of free-free and bound-free integrals
c***                   of unperturbed hamiltonian needed in kohn
c***                   calculations.
c***
c***                   for free-free integrals:computes m = (h - 0.5*kchan**2)
c***                                                      +     +
c***                   matrix elements of two classes: < h  m  h > and
c***                      +     -                      -
c***                   < h  m  h >. the definition of h changes depending
c***                   on whether t-matrix or s-matrix boundary conditions
c**                    are specified. for bound-free integrals:computes
c***                      +
c***                   < h (h-e) b >
c***                   for t-matrix boundary conditions the imaginary part
c***                                                     -
c***                   of the free function is used for h  which is
c***                    asymptotically a regular bessel function.
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6005
      program ffbf
      implicit integer (a-z)
      parameter (dimc=30 , dime=100 , dimlm=100 , dimbf=200)
      real*8 z, echan, energy, kchan, rmin, rmax, rdel, alpha, gamma
      real*8 rd26, ec
      character*8 drctv
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
      dimension list(dimbf)
      logical logkey, ptest, bsym, typ, posinp, vdrctv, zerov
      logical prntbf, prntff, prntv, prnty, prntb, prntbs, prntgr
      character *4096 ops
      character *1600 card
      character *24 filnm
      character *32 xform
      character *13 grdtyp, grdchk
      character *128 filgrd, filorb, filpot, filylm, filbes
      character *128  filint, filkne
      character *13 chrkey, srtfl
      character *10 cpass
      character*8  bcondx
      character *16 fptoc
      character *3 itoc
c----------------------------------------------------------------------c
c               ancient history                                        c
c               code now dynamically dimensioned                       c
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
c      (transformation to mo's occurs in scattering code)              c
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
      ptest=logkey(ops,'print=m6005=all',.false.,' ')
      if (ptest) then
          prntbf=.true.
          prntff=.true.
          prntv=.true.
          prnty=.true.
          prntb=.true.
          prntbs=.true.
          prntgr=.true.
      else
          prntbf=logkey(ops,'print=m6005=bound-free',.false.,' ')
          prntff=logkey(ops,'print=m6005=free-free',.false.,' ')
          prntv=logkey(ops,'print=m6005=potential',.false.,' ')
          prnty=logkey(ops,'print=m6005=ylm',.false.,' ')
          prntb=logkey(ops,'print=m6005=bound',.false.,' ')
          prntbs=logkey(ops,'print=m6005=bessel',.false.,' ')
          prntgr=logkey(ops,'print=m6005=grid',.false.,' ')
      endif
      zerov=logkey(ops,'zero-potential',.false.,' ')
      srtfl='"sorted orbs"'
      typ=.true.
      bcondx='t-matrix'
      filgrd='grid'
      filorb='orbs'
      filpot='pot'
      filylm='ylm'
      filbes='bessfn'
      filint='knints'
      if( posinp('$kohnint',cpass) ) then
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
c----------------------------------------------------------------------c
c         calculate potential integrals or only one electron types     c
c----------------------------------------------------------------------c
      call iosys ('read character "grid filename" from rwf',-1,
     1             0,0,filgrd)
      call iosys ('read character "orbital filename" from rwf',-1,
     1             0,0,filorb)
      call iosys ('read character "potential filename" from rwf',
     1             -1,0,0,filpot)
      call iosys ('read character "spherical harmonic filename" '//
     1            'from rwf',-1,0,0,filylm)
      call iosys ('read character "bessel function filename" '//
     1            'from rwf',-1,0,0,filbes)
      call iosys ('read character "kohn integral filename" from rwf',-1,
     1             0,0,filint)         
      call iosys ('read character "kohn data filename" from rwf',0,0,0,
     1             filkne)
      write (iout,1)
      if (zerov) then
          write(iout,*) '     zeroing potential'
      endif
      call iosys ('open kohndt as old',0,0,0,filkne)
      call iosys ('read character "integral directive" from kohndt',-1,
     1             0,0,drctv)
      vdrctv=.false.
      if (drctv.eq.'one') then
          write (iout,3000)
          vdrctv=.true.
      endif
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
          call lnkerr('mismatch in type grids on orbs and ylm')
      endif
      if (ngrid.ne.ngchk) then
          call lnkerr('mismatch in no. points on orbs and ylm')
      endif
      if (pntbuf.ne.pntchk) then
          call lnkerr('mismatch in point buffer on orbs and ylm')
      endif
      if (.not.vdrctv) then
           call iosys ('open vstat as old',0,0,0,filpot)
      endif
      call iosys ('read integer "no. states" from kohndt',1,
     1             nchan,0,' ')
      nstri=nchan*(nchan+1)/2
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
              call iosys ('write character symmetry to kohnint',0,0,
     1                     0,'on')
          else
              call iosys ('write character symmetry to kohnint',0,0,
     1                     0,'off')
          endif
          maxlm=0
          maxbf=0
          do 60 ch1=1,nchan
             if( posinp('$chan-'//itoc(ch1),cpass) ) then
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
      call iosys ('read integer "total no. mos" from kohndt',1,nmotot,
     1             0,' ')
      exdim=nchan
      if (vdrctv) then
          exdim=1
      endif
      hpvbwd=2*maxlm*nkept*nchan*exdim
      ovbfwd=2*maxlm*nkept*nchan
      ovbfww=2*maxlm*nmotot*nchan
      hpvbww=2*maxlm*nmotot*nchan*exdim
      ovppw=0
      ovpmw=0
      hpvhpw=2*maxlm*maxlm*nstri
      hpvhmw=2*maxlm*maxlm*nchan*exdim
      beswdr=nr*(lmax+1)
      beswds=2*beswdr
      newwds=2*pntbuf*maxlm*nchan
      vwords=pntbuf*nstri
      ylmwds=pntbuf*(lmax+1)*(mumax+mumax+1)
      baswds=nkept*pntbuf
      if (vdrctv) then
          hpvhpw=2*maxlm*nchan
          hpvhmw=0
          ovppw=hpvhpw
          ovpmw =ovppw
          vwords=0
      endif
      words=ncon*nmotot+hpvbwd+ovbfwd+hpvhpw+hpvhmw+
     1      ovbfww+hpvbww+2*maxlm*max(maxlm,ncon)*nchan*exdim+
     2      maxlm*max(maxlm,nkept)+iadtwp(lmax+1)+nr+
     3      2*beswdr+2*beswds+2*newwds+6*pntbuf+vwords+ylmwds+
     4      baswds+4*pntbuf*maxlm+max(4*pntbuf*maxlm,2*pntbuf*nkept)+
     5      maxlm*max(2,pntbuf)+ovppw+ovpmw
      if (words.gt.maxcor) then
          call lnkerr('not enough memory for calculation:quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6005',0)
      write (iout,300) words
      vec=ioff
      hpvb=vec+ncon*nmotot
      ovbf=hpvb+hpvbwd
      hpvhp=ovbf+ovbfwd
      hpvhm=hpvhp+hpvhpw
      ovpp=hpvhm+hpvhmw
      ovpm=ovpp+ovppw
      ovbm=ovpm+ovpmw
      hpvm=ovbm+ovbfww
      scrc=hpvm+hpvbww
      scrr=scrc+2*maxlm*max(maxlm,ncon)*nchan*exdim
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
      call iosys('read character "transformation vector" from kohndt',
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
c                   zero matrix elements used in all cases             c
c----------------------------------------------------------------------c
         call rzero(z(hpvhp),hpvhpw)
         call rzero(z(hpvb),hpvbwd)
         call rzero(z(ovbf),ovbfwd)
c----------------------------------------------------------------------c
c              rewind all energy-independent files                     c
c----------------------------------------------------------------------c
         call iosys ('rewind all on grid read-and-write',0,0,0,' ')
         call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
         call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
         if (.not.vdrctv) then
             call iosys ('rewind all on vstat read-and-write',0,0,
     1                    0,' ')
c----------------------------------------------------------------------c
c             hpvhm is only needed when the potential is needed        c
c----------------------------------------------------------------------c  
             call rzero(z(hpvhm),hpvhmw)
         else
c----------------------------------------------------------------------c
c             free-free overlaps only needed for special cases         c
c----------------------------------------------------------------------c
             call rzero(z(ovpp),ovppw)
             call rzero(z(ovpm),ovpmw)
         endif
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
            call mkbes (z(hs),z(hsder),z(cj),z(cy),z(hp),z(hd),z(grid),
     1                  z(rvec),z(x),z(krvec),kchan,rmin,rdel,rd26,
     2                  a(ipow),nlm,lch,npnts,nr,lmax,maxlm,nchan,
     3                  dimlm,dimc)
c----------------------------------------------------------------------c
c           read static potentials for this block                      c
c----------------------------------------------------------------------c
            if (.not.vdrctv) then
                call iosys ('read real "static potential" from'//
     1                      ' vstat without rewinding',npnts*nstri,
     2                        z(vpot),0,' ')
                if (zerov) then
                    call rzero(z(vpot),npnts*nstri)
                endif
                if (prntv.and..not.zerov) then
                    call vprnt(z(vpot),npnts,nstri,ireg)
                endif
c----------------------------------------------------------------------c
c             multiply integration weight into potentials              c
c----------------------------------------------------------------------c
                call scalev(z(vpot),z(grid),npnts,nstri)
            endif
            if (prntgr) then
                call grdprn(z(grid),npnts,ireg)
            endif
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
c     at this point a block of everything on the grid has been         c
c     read in.  now call routines which accumulate the bound-free      c
c     and free-free matrix elements from contributions at              c
c                        each grid point                               c
c----------------------------------------------------------------------c
            if (vdrctv) then
                call ffint1(z(ovpp),z(ovpm),z(hpvhp),z(grid),z(hp),
     1                      z(hd),z(ylm),z(cmatt),z(rmat),z(rtrick),
     2                      z(dmat),nlm,lch,mch,npnts,nchan,lmax,
     3                      dimlm,dimc,maxlm)
                call bfint1(z(hpvb),z(ovbf),z(grid),z(basis),z(hp),
     1                      z(hd),z(ylm),z(cmat),z(rmat),z(dmat),
     2                      z(scrc),z(scrr),nlm,lch,mch,ngauss,ngch,
     3                      npnts,nchan,nkept,list,lmax,dimlm,dimc,
     4                      dimbf,maxlm)
            else
                call ffints(z(hpvhp),z(hpvhm),z(grid),z(hp),z(hd),
     1                      z(ylm),z(vpot),z(cmatt),z(cmat),z(dmat),
     2                      z(rmat),z(rtrick),z(scrc),nlm,lch,mch,npnts,
     3                      nchan,lmax,nstri,dimlm,dimc,maxlm,bcondx)
                call bfints(z(hpvb),z(ovbf),z(grid),z(basis),z(hp),
     1                      z(hd),z(ylm),z(vpot),z(cmat),z(rmat),
     2                      z(dmat),z(scrc),z(scrr),nlm,lch,mch,ngauss,
     3                      ngch,npnts,nchan,nkept,list,lmax,nstri,
     4                      dimlm,dimc,dimbf,maxlm)
            endif
  430    continue
         filnm='ch-en-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',nchan,kchan,
     1               0,' ')
c-----------------------------------------------------------------------c
c              transform to molecular orbital basis                     c
c-----------------------------------------------------------------------c
         call tomobs(z(vec),z(ovbf),z(ovbm),z(hpvb),z(hpvm),z(scrc),
     1               maxlm,ncon,nmotot,nchan,ngauss,ngch,nkept,
     2               list,dimbf,dimc,exdim)
         filnm='ovlp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',ovbfww,z(ovbm),
     1               0,' ')
         filnm='bf-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',hpvbww,z(hpvm),
     1               0,' ')
         filnm='ffp-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',hpvhpw,
     1               z(hpvhp),0,' ')
         if (vdrctv) then
             filnm='ppo-'//fptoc(energy(ene))
             call iosys ('write real '//filnm//' to kohnint',ovppw,
     1                    z(ovpp),0,' ')
             filnm='pmo-'//fptoc(energy(ene))
             call iosys ('write real '//filnm//' to kohnint',ovpmw,
     1                    z(ovpm),0,' ')
         else
             filnm='ffm-'//fptoc(energy(ene))
             call iosys ('write real '//filnm//' to kohnint',hpvhmw,
     1                   z(hpvhm),0,' ')         
         endif
c----------------------------------------------------------------------c
c                    output section                                    c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c      write bound-free matrices to printer if requested               c
c----------------------------------------------------------------------c
         if( prntbf) then
             write (iout,1000)
             call bfprnt(z(ovbf),z(hpvb),nlm,nchan,dimc,maxlm,
     1                   nkept,exdim)
             write (iout,2000)
             call bfprnt(z(ovbm),z(hpvm),nlm,nchan,dimc,maxlm,
     1                   nmotot,exdim)
         endif
         if (prntff) then
             call ffprnt(z(hpvhp),z(hpvhp),z(hpvhm),z(hpvhm),z(ovpp),
     1                   z(ovpm),nlm,nchan,nstri,dimc,maxlm,vdrctv)
         endif
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
      call iosys ('rewind all on kohnint read-and-write',0,0,0,' ')
      call iosys ('close kohnint',0,0,0,' ')
      call iosys ('rewind all on bessel read-and-write',0,0,0,' ')
      call iosys ('close bessel',0,0,0,' ')
      call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
      call iosys ('close orbs',0,0,0,' ')
      call iosys ('rewind all on ylms read-and-write',0,0,0,' ')
      call iosys ('close ylms',0,0,0,' ')
      if (.not.vdrctv) then
          call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
          call iosys ('close vstat',0,0,0,' ')
      endif
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
    1 format(//,20x,'*****  m6005:bound-free and free-free integrals ***
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
  300 format(/,5x,'need',1x,i8,1x,'words for calculation')
  420 format(/,5x,'incident energy(ry)',1x,f10.6,/,5x,'channel momenta',
     1       (/,10x,5(e12.5,1x)))
 1000 format(/,5x,'atomic orbital representation')
 1500 format(/,5x, 'will use short orbital list')
 2000 format(/,5x,'molecular orbital representation')
 3000 format (/,5x,'only one-electron integrals will be computed')
      call chainx(0)
      stop
      end
