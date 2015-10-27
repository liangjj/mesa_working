*deck @(#)pm505.f	3.1  11/20/92
      subroutine pm505(z,a)
c
c***begin prologue     m505
c***date written       840727   (yymmdd)
c***revision date      890728   (yymmdd)
c
c   1 august 1990      rlm at lanl
c       adding symmetry restricted scf capabilities.
c  28 july  1989       bhl at llnl
c       adding scf=punch option to write occupied orbitals
c       so they can me modified and then read by m401
c  05 oct 1988       bhl at llnl
c       adding general scf code
c  21 march 1987       pws at lanl
c       adding level shifting options: level-shift and level-shift1
c  21 february 1987    pws at lanl
c       reducing core requirements by not holding all the diis information
c       for long past iterations. to do this, i am introducing a new
c       option (maxdiis-->'mxds') which is the maximum number of
c       iterations to hold for diis purposes (default is 8). this also
c       adds a couple of helper arrays: 'dsptr' which is 'mxds' long
c       and holds the iteration number of the data in that diis storage
c       bin, and 'ptr', 'mxiter' long, which holds which diis bin a
c       particular iterations data is in.
c  17 february 1987    pws at lanl
c       fixing up core allocation for 32/64 bit machines, especially
c       computation of 'ntriang'.
c   3 december 1986    pws at lanl
c       changing 'namint' and the iosys open to character.
c
c***keywords           m505, link 505, hartree-fock, scf,
c                      print, vector, convergence, scfcyc, pseudocanonical,
c                      extrapolate, diagonals, startdiis
c***author             saxe, paul (lanl)
c***source             @(#)pm505.f	3.1   11/20/92
c***purpose            solves the hartree-fock equations.
c***description
c     m505 recognizes the options subtrings:
c     timing                 collect and print timing statistics
c     level-shift1=0.0        minimum value of (approximate) hessian elements
c     level-shift=0.0       amount to add to hessian elements.
c     print_energy           print the energy each iteration.
c     print_vector           print the scf vector
c     print_mo_lagrangian    print the lagrangian in the mo basis.
c     print_ao_lagrangian    print the lagrangian in the ao basis.
c     convergence=n          diis error for convergence, as 10**(-n)
c                             (default 10**-8)
c     scfcyc=n               maximum number of iterations
c     extrapolate=n          1 to extrapolate fock matrices
c                            0 not (default)
c     diagonals=n            1 to use f(i,i)
c                            0 to use i (default)
c     startdiis=n            convergence at which to initiate diis, as
c                            10**(-n)    (default 10**4)
c     maxdiis=n              number of previous iterations to consider
c                            in the diis procedure (default 8).
c
c
c***references
c
c***routines called
c***end prologue       m505
c
      implicit integer (a-z)
c
      parameter (maxnbf=2000,maxrep=14)
c
      character*8 prtflg
      character*16 bflabl(maxnbf)
      character*8 lirrep(maxrep)
      character*128 namint,namchk
      character*4096 ops
      character*8 calc
      character*3 answer
      real*8 z(*)
      real*8 fpkey,level,level1
      integer a(*)
      logical prnt,logkey,gscf
      logical usesym
      logical pulay
c
      common /io/     inp,iout
c
      data maxcor /1/
      data prnt/.true./
      data usesym/.false./
      save maxcor,prnt,usesym
c
    2 format(1x,'m505:')
    3 format(5x,'memory use',18x,i9)
    4 format(5x,'all integrals held in core.')
    5 format(5x,'# integral triangles in core',i4)
c
c     ----- open the check file ------
c
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as unknown',0,0,0,namchk)
c
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- has printing been turned off externally? -----
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
      if(prtflg.eq.'minimum') prnt=.false.
c
c     ----- get the dimensions, etc -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $           1,nbf,0,' ')
      call iosys('read integer "number of auxiliary basis functions"'
     $           //' from rwf',1,nbfx,0,' ')
      nnp=(nbf+1)*nbf/2
      nnpx=(nbfx+1)*nbfx/2
      call iosys('read integer "spin multiplicity" from rwf',
     $           1,multip,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $           1,natoms,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $           1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $           1,nbe,0,' ')
c
      if (nbf.gt.maxnbf) then
         call lnkerr('character core not long enough')
      end if
c
c     ----- set up some parameters depending on multip -----
c
      gscf=logkey(ops,'scf=gscf',.false.,' ')
      pulay=logkey(ops,'scf=pulay',.false.,' ')
c
      if (multip.eq.1 .and. (.not.gscf)) then
         if (logkey(ops,'hf=tcscf',.false.,' ')) then
c
c           ----- two-configuration scf (gvb) -----
c
            calc='gvb'
            nfock=3
            ncoul=3
            nexch=3
            nshell=4
            ndmat=3
            if (pulay) then
               nfock=4
            end if
         else
c
c           ----- closed-shell rhf -----
c
            calc='closed'
            nfock=1
            ncoul=1
            nexch=1
            nshell=2
            ndmat=1
         end if
      else if ((pulay).and.(nbe.ne.0)) then
         calc='open'
         nfock=2
         ncoul=2
         nexch=2
         nshell=3
         ndmat=2
      else if(gscf.or.nbe.eq.0) then
c
c        ----- general scf -----
c
         calc='general'
         if(nbe.ne.0) then
            nshell=2
         else
            nshell=1
         end if
         nshell=intkey(ops,'scf=nshell',nshell,' ')
         nfock=intkey(ops,'scf=nfock',nshell,' ')
         ncoul=intkey(ops,'scf=ncoul',nshell,' ')
         nexch=intkey(ops,'scf=nexch',nshell,' ')
         ndmat=intkey(ops,'scf=ndmat',nshell,' ')
c
c        nshell is the total number of shells including the virtuals
c        input value is the number of occupied shells
         nshell=nshell+1
      else
         calc='open'
         nfock=2
         ncoul=2
         nexch=2
         nshell=3
         ndmat=2
      end if
c
      call iosys('write integer "number of shells" to rwf',1,nshell,
     $            0,' ')
      call iosys('write integer "number of hf density matrices" '//
     $            'to rwf',1,ndmat,0,' ')
c
c     ----- check on options and set defaults -----
c
      mxiter=intkey(ops,'scf=cycles',30,' ')
      extrap=intkey(ops,'scf=extrapolate',0,' ')
      diagnl=intkey(ops,'scf=diagonals',0,' ')
      stdiis=intkey(ops,'scf=startdiis',-4,' ')
      finish=intkey(ops,'scf=convergence',8,' ')
      level=fpkey(ops,'scf=level-shift',0.0d+00,' ')
      level1=fpkey(ops,'scf=level-shift1',0.0d+00,' ')
      if (level.ne.0.0d+00) level1=0.0d+00
      mxds=intkey(ops,'scf=maxdiis',5,' ')
      usesym=logkey(ops,'scf=symmetry',.false.,' ')
      if(usesym) then
         call iosys('does "number of irreducible representations"'
     $              //' exist on rwf',0,0,0,answer)
         if(answer.eq.'yes') then
            call iosys('read integer "number of irreducible'
     $              //'representations" from rwf',1,nirrep,0,' ')
         else
            write(iout,440)
  440       format(5x,'symmetry has been turned off in m202')
            nirrep=1
            usesym=.false.
         endif
      endif
c
c     ----- allocate core -----
c
      ptr=1
      dsptr=ptr+mxiter
      shlnbf=dsptr+mxds
      shlmin=shlnbf+nshell
      shlmax=shlmin+nshell
      fcoef=iadtwp(shlmax+nshell)
      alpha=fcoef+nshell
      beta=alpha+nshell**2
      s=beta+nshell**2
      t=s+nnp
      v=t+nnp
      h=v+nnp
      smhalf=h+nnp
      f=smhalf+nnp
      u=f+nnp
      d=u+nbf**2
      cdiis=d+nnp*ndmat
      c=cdiis+nbf**2
      eigval=c+nbf**2
      flast=eigval+nbf
      dlast=flast+nnp
      salc=dlast+nnp
      numso=wpadti(salc+nbf*nbf)
      lambda=numso+nirrep
      orbsym=lambda+nirrep
      occsym=orbsym+nbf
      kept=occsym+nirrep*nshell
      temp=kept+nshell
      t1=iadtwp(temp+2*nbf)
      t2=t1+max(nbf**2,nnp*ncoul)
      t3=t2+max(nbf**2,nnp*nexch)
      t4=t3+nbf**2
      t5=t4+nbf**2*nexch
      triang=max(t5+nbf**2,t1+mxds**2+mxds)
      jmat=triang+nnp
      kmat=jmat+nnp*ncoul
      error=kmat+nnp*nexch
      focksv=error+nbf**2*mxds
      energs=focksv+nnp*mxds
      diissv=energs+mxiter
      values=diissv+mxiter
c
c     ----- find out how much core is available -----
c
      call getscm(0,z,canget,'m505:',0)
c
      wish=wpadti(values+nnp**2)
      if (wish.le.canget) then
         ntriang=nnp
      else
         left=min(30000,iadtwp(canget)-values-nnp)
         ntriang=left/nnp+1
      end if
      need=wpadti(values+nnp*ntriang)
 
c     ----- get the core needed -----
c
      call getscm(need,z,maxcor,'m505',1)
c
      lenbuf=ntriang*nnp
c
      if(prnt) then
         write(iout,2)
         write(iout,3) maxcor
         if(ntriang.eq.nnp) then
            write(iout,4)
         else
            write(iout,5) ntriang
         endif
c
         if (level.ne.0.0d+00) write (iout,45) level
   45    format (t5,'level-shift used:',f5.2)
         if (level1.ne.0.0d+00) write (iout,46) level1
   46    format (t5,'minimum level-shift used:',f5.2)
      endif
      if (ntriang.lt.1) call lnkerr('not enough core in m505')
c
c     ----- open the integral file, copying to the ssd if needed -----
c
      call iosys('read character "integral filename" from rwf',
     $     0,0,0,namint)
c
      if (logkey(ops,'scf=ssd',.false.,' ')) then
         if (prnt) then
            write (iout,430)
 430        format (t5,'copying integrals to the ssd')
         end if
         call iosys('open oldints as old',0,0,0,namint)
         call iosys('open ints as scratch on ssd',nnp**2/50,0,0,' ')
         call iosys('copy "sorted ao integrals" from oldints '//
     $        'to ints',need,a,0,' ')
         call iosys('close oldints',0,0,0,' ')
      else
         call iosys('open ints as old',0,0,0,namint)
      end if
c
c     ----- read in s, t and v one-electron integrals -----
c
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $           nnp,z(t),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $           nnp,z(v),0,' ')
c
c     ----- form s**(1/2) and s**(-1/2) -----
c
      call sinv(z(s),z(smhalf),z(u),z(eigval),z(t1),z(t2),
     $          nbf,nnp,z(triang),iprint)
c
c     ----- initialize the arrays -----
c
      call setup(a(shlnbf),a(shlmin),a(shlmax),nshell,nbf,z(alpha),
     $     z(beta),z(fcoef),calc,ops,z(t1),z(t2),nocc,a(occsym),
     $     a(numso),a(lambda),lirrep,nirrep,usesym,a(orbsym),z(salc))
c
c     ----- do the scf iterations -----
c
      call scf(z(s),z(t),z(v),z(h),z(f),z(u),z(eigval),z(t1),z(t2),
     $         z(t3),z(triang),z(values),z(c),z(smhalf),
     $         nbf,nnp,lenbuf,itap44,z(d),z(t4),z(t5),z(error),
     $         z(focksv),mxiter,finish,stdiis,nocc,iprint,z(dlast),
     $         z(flast),ntriang,z(t1),nshell,a(shlmin),
     $         a(shlmax),z(alpha),z(beta),pulay,calc,ncoul,nexch,ndmat,
     $         nfock,extrap,diagnl,z(energs),z(diissv),prnt,ops,
     $         mxds,a(ptr),a(dsptr),level,level1,z(cdiis),z(fcoef),
     $         z(jmat),z(kmat),bflabl,usesym,z(salc),a(numso),a(lambda),
     $         lirrep,nirrep,a(orbsym),a(occsym),a(temp),a(kept),
     $         nbfx,nnpx)
c
c     ----- form the lagrangian -----
c           set up some parameters depending on multip -----
c
      if (calc.eq.'closed') then
         nfock=1
         ncoul=1
         nexch=1
         nshell=2
         ndmat=1
      else if (calc.eq.'open') then
         nfock=2
         ncoul=2
         nexch=2
         nshell=3
         ndmat=2
      else if (calc.eq.'gvb') then
         nfock=3
         ncoul=3
         nexch=3
         nshell=4
         ndmat=3
      end if
c
c     ----- allocate core -----
c
      h=beta+nshell**2
      c=h+nnp
      d=c+nbf**2
      j=d+nnp*ndmat
      k=j+nnp*ncoul
      t1=k+nnp*nexch
      t2=t1+nbf**2
      lag=t2+nbf**2*ndmat
      values=lag+nbf**2
c
c     ----- find out how much core is available -----
c
      call getscm(0,z,canget,'m505: how much core',0)
c
      wish=wpadti(values+nnp**2)
      if (wish.le.canget) then
         ntriang=nnp
      else
         left=min(30000,iadtwp(canget)-values-nnp)
         ntriang=left/nnp+1
      end if
      need=wpadti(values+nnp*ntriang)
 
c     ----- get the core needed -----
c
      call getscm(need,z,maxcor,'m505',1)
c
      lenbuf=ntriang*nnp
c
      if (ntriang.lt.1) call lnkerr('not enough core in l505')
c
c     ----- read in t and v one-electron integrals -----
c
      call iosys('read real "kinetic integrals" from rwf',
     $           nnp,z(values),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $           nnp,z(h),0,' ')
      call vadd(z(h),z(h),z(values),nnp)
c
c     ----- form the lagrangian -----
c
      call lagrng(a(shlnbf),a(shlmin),a(shlmax),z(fcoef),
     $            z(alpha),z(beta),
     $            z(h),z(c),z(d),z(j),z(k),z(values),nbf,nnp,nshell,
     $            ndmat,ncoul,nexch,ntriang,z(t1),z(t2),
     $            z(lag),ops,nbfx,nnpx)
c
c     ----- delete the temporary integral file -----
c
      if (logkey(ops,'scf=ssd',.false.,' ')) then
         call iosys('destroy ints',0,0,0,' ')
      end if
c
c     ----- and exit gracefully -----
c
      call chainx(0)
c
c
      stop
      end
