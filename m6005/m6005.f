*deck m6005
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
          call iosys ('does "grid filename" exist on rwf',-1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "grid filename" from rwf',-1,
     1                     0,0,filgrd)
          else
              filgrd=chrkey(card,'grid-file-name','grid',' ')
          endif
          call iosys ('does "orbital filename" exist on rwf',-1,0,0,
     1                 filorb)
          if (ans.ne.'no') then
              call iosys ('read character "orbital filename" from rwf',
     1                     -1,0,0,filorb)
          else
              filorb=chrkey(card,'numerical-orbital-file-name',
     1                      'orbs',' ')
          endif
              call iosys ('does "potential filename" exist on rwf',-1,
     1                     0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "potential filename" from rwf',
     1                     -1,0,0,filpot)
          else
              filpot=chrkey(card,'numerical-potential-file-name',
     1                      'vstat',' ')
          endif
          call iosys ('does "spherical harmonic filename" exist on rwf',
     1                -1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "spherical harmonic '//
     1                    'filename" from rwf',-1,0,0,filylm)
          else
              filylm=chrkey(card,'ylm-file-name','ylms',' ')
          endif
          call iosys ('does "bessel function filename" exist on rwf',
     1                 -1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "bessel function filename" '//
     1                     'from rwf',-1,0,0,filbes)
          else
              filbes=chrkey(card,'bessel-file-name','bessel',' ')
          endif
          call iosys ('does "kohn integral filename" exist on rwf',
     1                -1,0,0,ans)
          if (ans.ne.'no') then
              call iosys ('read character "kohn integral filename" '//
     1                     'from rwf',-1,0,0,filint)
          else
              filint=chrkey(card,'kohn-integral-file-name','knints',' ')
          endif
          call locase(filint,filint)
      endif
      call iosys ('read character "kohn data filename" from rwf',0,0,0,
     1             filkne)
      write (iout,1)
      call iosys ('write character "kohn integral filename" to rwf',
     1             0,0,0,filint)
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
      call iosys ('read real "inner cutoff" from bessel',1,alpha,0,' ')
      call iosys ('read real "outer exp cutoff" from bessel',1,
     1            gamma,0,' ')
      call iosys ('read integer "outer n cutoff" from bessel',1,
     1            ncut,0,' ')
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
      call iosys ('open vstat as old',0,0,0,filpot)
      call iosys ('read character "grid type" from vstat',-1,0,0,grdchk)
      call iosys ('read integer "no. states" from vstat',1,nchan,0,' ')
      nstri=nchan*(nchan+1)/2
      if (grdchk.ne.grdtyp) then
          call lnkerr('mismatch in type grids on orbs and vstat')
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
      hpvbwd=2*maxlm*nkept*nchan*nchan
      ovbfwd=2*maxlm*nkept*nchan
      ovbfww=2*maxlm*nmotot*nchan
      hpvbww=2*maxlm*nmotot*nchan*nchan
      hpvhpw=2*maxlm*maxlm*nstri
      hpvhmw=2*maxlm*maxlm*nchan*nchan
      beswds=2*nr*(lmax+1)
      newwds=2*pntbuf*maxlm*nchan
      vwords=pntbuf*nstri
      ylmwds=pntbuf*(lmax+1)*(mumax+mumax+1)
      baswds=nkept*pntbuf
      words=ncon*nmotot+hpvbwd+ovbfwd+hpvhpw+hpvhmw+
     1      ovbfww+hpvbww+2*maxlm*max(maxlm,ncon)*nchan*nchan+
     2      maxlm*max(maxlm,nkept)+iadtwp(lmax+1)+nr+
     3      4*beswds+2*newwds+6*pntbuf+vwords+ylmwds+
     4      baswds+4*pntbuf*maxlm+max(4*pntbuf*maxlm,2*pntbuf*nkept)
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
      ovbm=hpvhm+hpvhmw
      hpvm=ovbm+ovbfww
      scrc=hpvm+hpvbww
      scrr=scrc+2*maxlm*max(maxlm,ncon)*nchan*nchan
      ipow=wpadti(scrr+maxlm*max(maxlm,nkept))
      x=iadtwp(ipow+lmax+1)
      hs=x+nr
      hsder=hs+beswds
      cj=hsder+beswds
      cy=cj+beswds
      hp=cy+beswds
      newwds=2*pntbuf*maxlm*nchan
      hd=hp+newwds
      rvec=hd+newwds
      krvec=rvec+pntbuf
      grid=krvec+pntbuf
      vpot=grid+4*pntbuf
      ylm=vpot+vwords
      basis=ylm+ylmwds
      cmat=basis+baswds
      dmat=cmat+4*pntbuf*maxlm
      call mkpow(a(ipow),lmax)
      call iosys('read character "transformation vector" from kohndt',
     $           -1,0,0,xform)
      call iosys ('read real '//xform//' from kohndt',ncon*nmotot,
     1             z(vec),0,' ')
      call iosys ('read real points from bessel',nr,z(x),0,' ')
      call iosys('read real hs from bessel',beswds,z(hs),0,' ')
      call iosys('read real hsder from bessel',beswds,z(hsder),0,' ')
      call iosys('read real cj from bessel',beswds,z(cj),0,' ')
      call iosys('read real cy from bessel',beswds,z(cy),0,' ')
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
         call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
         call iosys ('rewind all on orbs read-and-write',0,0,0,' ')
c----------------------------------------------------------------------c
c                   zero matrix elements                               c
c----------------------------------------------------------------------c
         call rzero(z(hpvhp),hpvhpw)
         call rzero(z(hpvhm),hpvhmw)
         call rzero(z(hpvb),hpvbwd)
         call rzero(z(ovbf),ovbfwd)
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
c     at this point a block of everything on the grid has been         c
c     read in.  now call routines which accumulate the bound-free      c
c     and free-free matrix elements from contributions at              c
c                        each grid point                               c
c----------------------------------------------------------------------c
            call ffints(z(hpvhp),z(hpvhm),z(grid),z(hp),z(hd),z(ylm),
     1                  z(vpot),z(cmat),z(cmat),z(dmat),z(scrc),z(scrr),
     2                  nlm,lch,mch,npnts,nchan,lmax,nstri,dimlm,dimc,
     3                  maxlm,bcondx)
            call bfints(z(hpvb),z(ovbf),z(grid),z(basis),z(hp),z(hd),
     1                  z(ylm),z(vpot),z(cmat),z(dmat),z(scrc),z(scrr),
     2                  nlm,lch,mch,ngauss,ngch,npnts,nchan,nkept,
     3                  list,lmax,nstri,dimlm,dimc,dimbf,maxlm)
  430    continue
         filnm='ch-en-'//fptoc(energy(ene))
         call iosys ('write real '//filnm//' to kohnint',nchan,kchan,
     1               0,' ')
c-----------------------------------------------------------------------c
c              transform to molecular orbital basis                     c
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
      call iosys ('rewind all on vstat read-and-write',0,0,0,' ')
      call iosys ('close vstat',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
    1 format(//,20x,'*****  m6005:bound-free and free-free integrals ***
     1**')
    3 format(/,5x,'target energies for ',i3,' channels:',/,(2x,5d15.8))
    4 format(/,5x,'incident energies:',/,(2x,5d15.8))
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
      call chainx(0)
      stop
      end
*deck mkbes
c***begin prologue     mkbes
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           m6005, link 6005, bessel functions
c***author             rescigno, tom (llnl)
c***source             m6005
c***purpose            bessel functions from spline fit
c***description        calculation of bessel functions from spline
c***                   fits in link m6004
c***routines called    iosys, util and mdutil
c***end prologue       mkbes
      subroutine mkbes(hs,hsder,cj,cy,hp,hd,grid,rvec,x,krvec,kchan,
     1                         rmin,rdel,rd26,ipow,nlm,lch,npt,nr,
     2                         lmax,maxlm,nchan,dimlm,dimc)
      implicit integer (a-z)
      real *8 grid, rvec, x, krvec, kchan, rmin, rdel, rd26, ksq, k52
      real *8 a, b
      complex *16 hs, hsder, cj, cy, hp, hd
      dimension grid(4,npt), rvec(npt), krvec(npt), kchan(nchan)
      dimension hs(nr,0:lmax), hsder(nr,0:lmax), cj(nr,0:lmax)
      dimension cy(nr,0:lmax), hp(npt,maxlm,nchan), hd(npt,maxlm,npt)
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
*deck ffints
c***begin prologue     ffints
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            free-free matrix elements
c***description        calculation of free-free integrals
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       ffints
      subroutine ffints(hpvhp,hpvhm,grid,hp,hd,ylm,vpot,cmat,cmatt,dmat,
     1                  scrc,scrr,nlm,lch,mch,npt,nchan,lmax,nstri,
     2                  dimlm,dimc,maxlm,bcond)
      implicit integer (a-z)
      character *(*) bcond
      real *8 ylm, vpot, grid, scrr
      complex *16 cmat, cmatt, dmat, hpvhp, scrc
      complex *16 hpvhm, hp, hd, cdotu
      dimension ylm(npt,0:lmax,0:2*lmax), vpot(npt,nstri)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt,2), dmat(maxlm,npt,2), grid(4,npt)
      dimension hpvhp(1:maxlm,1:maxlm,nstri), scrc(*), scrr(*)
      dimension hpvhm(1:maxlm,1:maxlm,nchan,nchan), cmatt(npt,maxlm,2)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
c
c compute potential matrices
c
c----------------------------------------------------------------------c
c         first two loops over channel 1 and associated l, m           c
c         to load some temporary matrices for vectorization            c
      do 10 ch1=1,nchan
         call czero(cmat(1,1,1),maxlm*npt)
         call czero(dmat(1,1,1),maxlm*npt)
         call filcf1(cmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm)
         call fildf1(dmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm,bcond)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c          next two loops over channel 2 and associated l, m           c
c          same strategy on temporary matrices                         c
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
c                  accumulate into matrix elements                     c
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
c            compute matrix of (kchan**2/2 - t)                        c
c            diagonal in channel indices so its simplier than above    c
c----------------------------------------------------------------------c
      do 30 ch1=1,nchan
         clst=ch1*(ch1+1)/2
         call czero(cmatt,maxlm*npt*2)
         do 40 nolm1=1,nlm(ch1)
            l1=lch(nolm1,ch1)
            m1=mch(nolm1,ch1)
            if(bcond.eq.'s-matrix') then
               do 50 grpt=1,npt
1                 cmatt(grpt,nolm1,1)=grid(4,grpt)*ylm(grpt,l1,m1)
     1                                *ylm(grpt,l1,m1)
     2                                *hp(grpt,nolm1,ch1)
                  cmatt(grpt,nolm1,2)=conjg(hd(grpt,nolm1,ch1))
   50          continue
            elseif (bcond.eq.'t-matrix') then
               do 60 grpt=1,npt
                   cmatt(grpt,nolm1,1)=grid(4,grpt)*ylm(grpt,l1,m1)
     1                                 *ylm(grpt,l1,m1)
     2                                 *hp(grpt,nolm1,ch1)
                   cmatt(grpt,nolm1,2)=aimag(hd(grpt,nolm1,ch1))
   60          continue
            endif
   40    continue
         do 70 nolm1=1,nlm(ch1)
            hpvhp(nolm1,nolm1,clst)=hpvhp(nolm1,nolm1,clst)+
     1                              cdotu(npt,cmatt(1,nolm1,1),1,
     2                                    hd(1,nolm1,ch1),1)
            hpvhm(nolm1,nolm1,ch1,ch1)=hpvhm(nolm1,nolm1,ch1,ch1)+
     1                                 cdotu(npt,cmatt(1,nolm1,2),1,
     2                                       cmatt(1,nolm1,1),1)
   70    continue
   30 continue
      return
      end
*deck bfints
c***begin prologue     bfints
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            bound-free matrix elements
c***description        calculation of bound-free integrals
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       bfints
      subroutine bfints(hpvb,ovbf,grid,basis,hp,hd,ylm,vpot,cmat,bmat,
     1                  scrc,scrr,nlm,lch,mch,ngauss,ngch,npt,nchan,
     2                  nkept,list,lmax,nstri,dimlm,dimc,dimbf,maxlm)
      implicit integer (a-z)
      real *8 ylm, vpot, grid, basis, scrr
      complex *16  cmat, bmat, hpvb, scrc
      complex *16 ovbf, hp, hd
      dimension ylm(npt,0:lmax,0:2*lmax), vpot(npt,nstri)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt,2), bmat(npt,nkept), grid(4,npt)
      dimension hpvb(1:maxlm,nkept,nchan,nchan), scrc(*)
      dimension ovbf(1:maxlm,nkept,nchan), scrr(*)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension basis(npt,nkept), list(dimbf)
      dimension ngauss(dimc), ngch(dimbf,dimc)
c----------------------------------------------------------------------c
c          first two loops over channel 1 and associated l, m          c
 
      do 10 ch1=1,nchan
         call czero(cmat(1,1,1),maxlm*npt)
         call filcf1(cmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             next two loops over channel 2 and associated             c
c                        bound orbitals                                c
         do 20 ch2=1,nchan
            call czero(bmat,npt*nkept)
            chlrg=max(ch1,ch2)
            chsml=ch1+ch2-chlrg
            itri=chlrg*(chlrg-1)/2+chsml
            do 30 bfn=1,ngauss(ch2)
               bfnch2 = ngch(bfn,ch2)
               blocte=list(bfnch2)
               if (blocte.eq.0) then
                   call lnkerr('bound orbital out of range')
               endif
c----------------------------------------------------------------------c
               do 40 grpt=1,npt
                  bmat(grpt,blocte)=vpot(grpt,itri)*basis(grpt,blocte)
   40          continue
   30       continue
            call vcebc(scrc,cmat(1,1,1),bmat,scrr,nlm(ch1),npt,nkept)
            call cvadd(hpvb(1,1,ch1,ch2),maxlm,scrc,nlm(ch1),nlm(ch1),
     1                 nkept)
   20    continue
   10 continue
c----------------------------------------------------------------------c
c            compute matrix of (e - t)                                 c
c----------------------------------------------------------------------c
      do 70 ch1=1,nchan
         call czero(cmat,maxlm*npt*2)
         call filc12(cmat,ylm,grid,hp(1,1,ch1),hd(1,1,ch1),nlm(ch1),
     1               lch(1,ch1),mch(1,ch1),maxlm,npt,lmax,dimlm)
         call czero(bmat,npt*nkept)
         do 80 bfn=1,ngauss(ch1)
            bfnch1=ngch(bfn,ch1)
            blocte=list(bfnch1)
            if (blocte.eq.0) then
                call lnkerr('bound orbital out of range')
            endif
            do 90 grpt=1,npt
               bmat(grpt,blocte)=basis(grpt,blocte)
   90       continue
   80    continue
         call vcebc(scrc,cmat(1,1,1),bmat,scrr,nlm(ch1),npt,nkept)
         call cvadd(ovbf(1,1,ch1),maxlm,scrc,nlm(ch1),nlm(ch1),nkept)
         call vcebc(scrc,cmat(1,1,2),bmat,scrr,nlm(ch1),npt,nkept)
         call cvadd(hpvb(1,1,ch1,ch1),maxlm,scrc,nlm(ch1),nlm(ch1),
     1              nkept)
   70 continue
      return
      end
*deck tomobs
c***begin prologue     tomobs
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           transformation, kohn
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            bound-free ao to mo transformation
c***
c***references         none
c
c***routines called    cebc
c***end prologue       tomobs
 
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
*deck mkpow
c***begin prologue     mkpow
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn, integrals
c***author             rescigno, tom (llnl)
c***source             m6005
c***purpose            fill array ipow
c***references
c
c***routines called    none
c***end prologue       mkpow
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
*deck bfprnt
c***begin prologue     bfprnt
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print bound-free matrix elements
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
            title='hpvb matrix ch1-'//itoc(ch1)//' ch2-'//itoc(ch2)
            call cmprir(hpvb(1,1,ch1,ch2),rowv,colv,nlm(ch1),
     1                  ncon,maxlm,ncon,title,rowt,colt,iout)
   20    continue
   10 continue
      do 50 ch1=1,nchan
         title='bound free overlap matrix ch1-'//itoc(ch1)
         call cmprir(ovbf(1,1,ch1),rowv,colv,nlm(ch1),ncon,
     1               maxlm,ncon,title,rowt,colt,iout)
   50 continue
      return
      end
*deck ffprnt
c***begin prologue     ffprnt
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print free-free matrix elements
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
            title='hpvhp matrix ch1-'//itoc(ch1)//' ch2-'//itoc(ch2)
            rowt='lm index'
            colt=rowt
            rowv=-99.
            colv=-99
            call cmprir(hpvhp(1,1,ist),rowv,colv,nlm(ch1),nlm(ch2),
     1               maxlm,maxlm,title,rowt,colt,iout)
   20    continue
   10 continue
      do 70 ch1=1,nchan
         do 80 ch2=1,nchan
            title='hpvhm matrix ch1-'//itoc(ch1)//' ch2-'//itoc(ch2)
            rowt='lm index'
            colt=rowt
            rowv=-99.
            colv=-99
            call cmprir(hpvhm(1,1,ch1,ch2),rowv,colv,nlm(ch1),nlm(ch2),
     1               maxlm,maxlm,title,rowt,colt,iout)
   80    continue
   70 continue
      return
      end
*deck grdprn
c***begin prologue     grdprn
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, print grid
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print co-ordinate grid
c***references         none
c
c***routines called    none
c***end prologue       grdprn
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
   10 format (/,5x,'grid points for region',1x,i4)
   20 format(/,5x,'     x     ','     y     ','     z     ')
      end
*deck vprnt
c***begin prologue     vprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, print potential
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print static potentials on grid
c***references         none
c
c***routines called    none
c***end prologue       vprnt
      subroutine vprnt(v,npnts,ntri,reg)
      implicit integer (a-z)
      character *80 title
      character *3 itoc
      character *8 colt, rowt
      common /io/ inp, iout
      real *8 v, rowv
      dimension v(npnts,ntri)
      title='potential matrices for region-'//itoc(reg)
      rowv=-99.
      colv=-99
      colt='state'
      rowt='point'
      call mprir(v,rowv,colv,npnts,ntri,npnts,ntri,title,rowt,colt,
     1            iout)
      return
      end
*deck yprnt
c***begin prologue     yprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, spherical harmonic print
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print spherical harmonics
c***references         none
c
c***routines called    none
c***end prologue       yprnt
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
      rowt='point'
      colt='l'
      rowv=-99.
      write (iout,10) reg
      mcnt=0
      do 20 m=0,mumax
         iup=2
         if (m.eq.0) then
             iup=1
         endif
         do 40 mm=1,iup
            title='mu-'//itoc(m)//' comp-'//itoc(mm)
            call mprir(ylm(1,m,mcnt),rowv,lval(m),npnts,lmax+1-m,npnts,
     1                 lmax+1,title,rowt,colt,iout)
            mcnt=mcnt+1
   40    continue
   20 continue
   10 format(/,5x,'ylm for region',1x,i4)
      return
      end
*deck bndprn
c***begin prologue     bndprn
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, print
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print bound orbitals
c***references         none
c
c***routines called    none
c***end prologue       bndprn
      subroutine bndprn(orbs,npnts,ncon,reg)
      implicit integer (a-z)
      common /io/ inp, iout
      character *80 title
      character *8 rowt, colt
      character *3 itoc
      real *8 orbs, rowv
      dimension orbs(npnts,ncon)
      title='bound orbitals region-'//itoc(reg)
      rowv=-99.
      colv=-99
      colt='con ao'
      rowt='point'
      call mprir(orbs,rowv,colv,npnts,ncon,npnts,ncon,title,rowt,colt,
     1            iout)
      return
      end
*deck bsprnt
c***begin prologue     bsprnt
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, bessel print
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            print channel "bessel" functions
c***routines called    none
c***end prologue       bsprnt
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
      rowt='point'
      colt='l'
      write (iout,10) reg
      do 20 ch1=1,nchan
         title='hp functions channel-'//itoc(ch1)
         call cmprir(hp(1,1,ch1),rowv,lch(1,ch1),npnts,nlm(ch1),npnts,
     1               maxlm,title,rowt,colt,iout)
         title='hd functions channel-'//itoc(ch1)
         call cmprir(hd(1,1,ch1),rowv,lch(1,ch1),npnts,nlm(ch1),npnts,
     1               maxlm,title,rowt,colt,iout)
   20 continue
      return
   10 format(/,5x,'free functions for region',1x,i4)
      end
*deck scalev
c***begin prologue     scalev
c***date written       890512   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            multiply potential by integration weights
c
c***references
c
c***routines called   none
c***end prologue      scalev
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
*deck filcf1
c***begin prologue     filcf1
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       filcf1
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
*deck filcf2
c***begin prologue     filcf2
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       filcf2
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
*deck fildf1
c***begin prologue     fildf1
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       fildf1
      subroutine fildf1(d,a,b,nd,l,m,maxl,npt,lmax,dim,bcond)
      implicit integer (a-z)
      complex *16  b, d
      real *8 a
      character *(*) bcond
      dimension d(nd,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim)
c----------------------------------------------------------------------c
c                 test boundary conditions                             c
c----------------------------------------------------------------------c
      if (bcond.eq.'s-matrix') then
          do 10 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 20 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*conjg(b(grpt,nolm1))
   20        continue
   10     continue
      elseif (bcond.eq.'t-matrix') then
           do 30 nolm1=1,nd
              l1=l(nolm1)
              m1=m(nolm1)
              do 40 grpt=1,npt
                 d(nolm1,grpt)=a(grpt,l1,m1)*aimag(b(grpt,nolm1))
   40        continue
   30     continue
      else
          call lnkerr('boundary condition error in ffints')
      endif
      return
      end
*deck fildf2
c***begin prologue     fildf2
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       fildf2
      subroutine fildf2(d,a,b,v,nd,l,m,maxl,npt,lmax,dim,bcond)
      implicit integer (a-z)
      complex *16 b, d
      real *8 a, v
      character *(*) bcond
      dimension d(nd,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim), v(npt)
c----------------------------------------------------------------------c
c                 test boundary conditions                             c
c----------------------------------------------------------------------c
      if (bcond.eq.'s-matrix') then
          do 10 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 20 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*v(grpt)*conjg(b(grpt,nolm1))
   20        continue
   10     continue
      elseif (bcond.eq.'t-matrix') then
          do 30 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 40 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*v(grpt)*aimag(b(grpt,nolm1))
   40       continue
   30    continue
      else
         call lnkerr('boundary condition error in ffints')
      endif
      return
      end
*deck filc12
c***begin prologue     filc12
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       filc12
      subroutine filc12(c,a,b,h1,h2,nc,l,m,maxl,npt,lmax,dim)
      implicit integer (a-z)
      complex*16 c, h1, h2
      real *8 a, b
      dimension c(nc,npt,2), a(npt,0:lmax,0:2*lmax), b(4,npt)
      dimension h1(npt,maxl), h2(npt,maxl), l(dim), m(dim)
      do 10 nolm1=1,nc
         l1=l(nolm1)
         m1=m(nolm1)
         do 20 grpt=1,npt
            c(nolm1,grpt,1)=a(grpt,l1,m1)*b(4,grpt)*h1(grpt,nolm1)
            c(nolm1,grpt,2)=a(grpt,l1,m1)*b(4,grpt)*h2(grpt,nolm1)
   20    continue
   10 continue
      return
      end
