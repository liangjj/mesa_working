*deck m6004
c***begin prologue     m6004
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6004, link 6004, spline
c***author             rescigno, t. n.(llnl)
c***source             m6004
c***purpose            spline fit of bessel functions
c***description        bessel functions calculated and then converted
c***                   to appropriate cutoff functions which are then
c***                   spline fit, coefficients put on file. this is a
c***                   stand alone code.  
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6004
      program fitbes
      implicit integer (a-z)
      parameter (acc=30)
      character *8 cpass, bessfl, chrkey, filkne
      character *4096 ops
      character *800 card
      logical logkey, prnbes, prnspn, green
      real *8 z, fpkey, alpha, rmin, rmax, gamma, rdel, rd26
      real *8 rl
      dimension z(1)
      common a(1)
      common /io/ inp, iout
      common /memory/ ioff
      equivalence(a,z)
c----------------------------------------------------------------------c
c  lbig = biggest physical l value allowed                             c 
c  nblokk = # pts. in large argument region for splines - don't change c
c  nsblok = # pts. in small argument region for splines - don't change c
c  ltop = top l for backwards recursion (maximum value)                c
c----------------------------------------------------------------------c
      call drum
      call iosys ('read character "kohn data filename" from rwf',0,0,0,
     1             filkne)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('open kohndt as old',0,0,0,filkne)
      prnbes=logkey(ops,'print=m6004=bessel',.false.,' ')
      prnspn=logkey(ops,'print=m6004=spline',.false.,' ')
      green=logkey(ops,'m6004=no-greens-function',.false.,' ')
      call posinp('$bess',cpass)
      call cardin(card)
      bessfl=chrkey(card,'bessel-file-name','bessel',' ')
      call locase(bessfl,bessfl)
      call iosys ('open bessel as new on ssd',262144,0,0,bessfl)
      call iosys ('write character "bessel function filename" to rwf',
     1             0,0,0,bessfl)
c---------------------------------------------------------------------c
c   note that rmin is the first spline knot point. if the real        c
c    grid goes below rmin, the cubic spline is being extrapolated     c
c                         (should be o.k.)                            c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      generate spline data                           c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c             llmax = max l for splined bessel fcns                   c
c             alpha is for inner cut-off (y fcns)                     c
c             gamma and ncut are for outer cut-off                    c
c             rmax(rmin) = max(min) argument for                      c
c                          splined bessel fcns                        c
c             nper is the number of points to use per unit interval   c
c                         in computing the spline coefficients        c
c---------------------------------------------------------------------c  
      write (iout,100)
      lmax=intkey(card,'maximum-l-for-bessel',10,' ')
      alpha=fpkey(card,'inner-cutoff',1.0d0,' ')
      call iosys ('read real rmin from kohndt',1,rmin,0,' ')
      call iosys ('read real rmax from kohndt',1,rmax,0,' ')
      gamma=fpkey(card,'outer-exponential-cutoff',1.0d0,' ')
      ncut=intkey(card,'outer-n-cutoff',6,' ')
      nper=intkey(card,'points-per-interval-for-spline-coefficients',
     1            4,' ')      
      ltop=intkey(card,'top-recursion-l',100,' ')       
c---------------------------------------------------------------------c
c            estimate starting l                                      c
c---------------------------------------------------------------------c
      rl=lmax*acc
      strtl=lmax+sqrt(rl)
      strtl=max(strtl,lmax+1)
      ltop=min(ltop,strtl)
c-----------------------------------------------------------------------c
c compute total number of points,nr, and number of points in backward   c
c                       recursion call, ns.                             c
c-----------------------------------------------------------------------c
      nr=(rmax-rmin)*nper
      rdel=(rmax-rmin)/(nr-1)
      rd26=rdel*rdel/6.
      ns=lmax/rdel+1
      if(ns.gt.nr)ns=nr
      write(iout,10) rmin,rmax,rdel,lmax,alpha,gamma,ncut
      write (iout,20) nr, ns, ltop
      dimf1=nr*(ltop+1)
      dimf2=nr*(lmax+1)
      words=iadtwp(lmax+1)+6*nr+4*dimf1+8*dimf2
      write (iout,30) words
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      if (words.gt.canget) then
          call lnkerr ('cannot get required memory:will quit')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(words,z,ngot,'m6004',0)
      ipow=wpadti(ioff)
      x=iadtwp(ipow+lmax+1)
      xinv=x+nr
      j=xinv+nr
      jp=j+dimf1
      y=jp+dimf1
      yp=y+dimf1
      norm=yp+dimf1
      toim=norm
      toim2=toim+nr
      hs=toim2+nr
      hsder=hs+2*dimf2
      cj=hsder+2*dimf2
      cy=cj+2*dimf2
      d=cy+2*dimf2
      call mkpow(a(ipow),lmax)
c----------------------------------------------------------------------c
c                 make grid                                            c
c----------------------------------------------------------------------c
      call mkx(z(x),z(xinv),rmin,rdel,nr)
c----------------------------------------------------------------------c
c                write out necessary stuff to bessel                   c
c----------------------------------------------------------------------c      
      call iosys ('write integer maximum-l to bessel',1,lmax,0,' ')
      call iosys ('write integer "total pts" to bessel',1,nr,0,' ')
      call iosys ('write real points to bessel',nr,z(x),0,' ')
      call iosys ('write integer "no. small pts" to bessel',1,ns,
     1            0,' ')
      call iosys ('write real rmin to bessel',1,rmin,0,' ')
      call iosys ('write real rmax to bessel',1,rmax,0,' ')
      call iosys ('write real "r spacing" to bessel',1,rdel,0,' ')
      call iosys ('write real "inner cutoff" to bessel',1,alpha,0,' ')
      call iosys ('write real "outer exp cutoff" to bessel',1,
     1            gamma,0,' ')
      call iosys ('write integer "outer n cutoff" to bessel',1,
     1            ncut,0,' ')
c----------------------------------------------------------------------c
c             small r, backward recursion                              c
c----------------------------------------------------------------------c
      call rcbesb(z(x),z(xinv),z(j),z(jp),z(y),z(yp),z(norm),ns,nr,lmax,
     1            ltop,prnbes)
c----------------------------------------------------------------------c
c             large r, forward recursion                               c
c----------------------------------------------------------------------c
      if (nr.gt.ns) then
          nptlft=nr-ns
          call rcbesf(z(x+ns),z(xinv+ns),z(j+ns),z(jp+ns),z(y+ns),
     1                z(yp+ns),nptlft,nr,lmax,ltop,prnbes)
      endif
c----------------------------------------------------------------------c
c              make the spline coefficients                            c
c----------------------------------------------------------------------c
      call mkspln(z(x),z(xinv),z(j),z(jp),z(y),z(yp),z(hs),z(hsder),
     1            z(cj),z(cy),z(d),z(toim),z(toim2),rdel,alpha,gamma,
     2            ncut,a(ipow),nr,lmax,ltop)
      nowds=2*dimf2
      if (prnspn) then
          call prspln(z(x),z(hs),z(hsder),z(cj),z(cy),nr,lmax)
      endif
      call iosys('write real hs to bessel',nowds,z(hs),0,' ')            
      call iosys('write real hsder to bessel',nowds,z(hsder),0,' ')            
      call iosys('write real cj to bessel',nowds,z(cj),0,' ')            
      call iosys('write real cy to bessel',nowds,z(cy),0,' ')            
   10 format(/,5x,'rmin spline',1x,f10.5,1x,'rmax spline',1x,f10.5,1x,
     1            'step size',1x,f10.5,/,30x,'max l',1x,i3,/,5x,
     2            'cutoff parameters:     alpha',1x,f6.3,1x,'gamma',1x,
     3            f6.3,1x,'ncut',1x,i3)
   20 format(/,5x,'total no. points',1x,i5,1x,/,5x,
     1       'no. points small r region',1x,i4,/,5x,
     2       'top recursion l',1x,i4)
   30 format (/,5x,'get',1x,i8,1x,'words of memory')
  100 format(//,20x,'***** m6004:splined ricatti-bessel functions *****'
     1)
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)
      stop
      end
