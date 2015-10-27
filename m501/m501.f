*deck  @(#)pm501.f	5.1 11/6/94
      program m501
c
c***begin prologue     m501
c***date written       840727   (yymmdd)
c***revision date      890728   (yymmdd)
c   1 august 1990      rlm at lanl
c       adding symmetry restricted scf capabilities.
c
c  28 july  1989       bhl at llnl
c       adding scf=punch option to write occupied orbitals
c       so they can me modified and then read by m401
c
c  05 oct 1988       bhl at llnl
c       adding general scf code
c
c  21 march 1987       pws at lanl
c       adding level shifting options: level-shift and level-shift1
c
c  21 february 1987    pws at lanl
c       reducing core requirements by not holding all the diis information
c       for long past iterations. to do this, i am introducing a new
c       option (maxdiis-->'mxds') which is the maximum number of
c       iterations to hold for diis purposes (default is 8). this also
c       adds a couple of helper arrays: 'dsptr' which is 'mxds' long
c       and holds the iteration number of the data in that diis storage
c       bin, and 'ptr', 'mxiter' long, which holds which diis bin a
c       particular iterations data is in.
c
c  17 february 1987    pws at lanl
c       fixing up core allocation for 32/64 bit machines, especially
c       computation of 'ntriang'.
c
c   3 december 1986    pws at lanl
c       changing 'namint' and the iosys open to character.
c
c***keywords           m501, link 501, hartree-fock, scf,
c                      print, vector, convergence, scfcyc, pseudocanonical,
c                      extrapolate, diagonals, startdiis
c***author             saxe, paul (lanl)
c***source              @(#)pm501.f	5.1 11/6/94
c***purpose            solves the hartree-fock equations.
c***description
c     m501 recognizes the options subtrings:
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
c     pseudocanonical=n      1 to use pseudocanonical orbitals (default)
c                            0 not
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
c***end prologue       m501
c
      implicit integer (a-z)
c
      parameter (maxnbf=2000)
c
      character*8 prtflg
      character*16 bflabl(maxnbf)
      character*128 namint,namchk
      character*4096 ops
      real*8 z1, z2
      real*8 fpkey,level,level1
      integer a1, a2, ngot(2)
      logical prnt,logkey,gscf
c
      common /io/     inp,iout
c
c
      data maxcor /1/
      data prnt/.true./
      data memtri/1000000/
      save maxcor,prnt
      pointer (p1,z1(1)), (p1,a1(1))
      pointer (p2,z2(1)), (p2,a2(1))
c
    2 format(1x,'m501:')
    3 format(5x,'memory use',18x,i9)
    4 format(5x,'all integrals held in core.')
    5 format(5x,'# integral triangles in core',i4)
c
      call drum
c     ----- recover the options string -----
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- open the check file ------
c
      if(logkey(ops,'scf=symmetry',.false.,' ')) then
         call iosys('read character "checkpoint filename" from rwf',
     $               0,0,0,namchk)
         call iosys('open chk as unknown',0,0,0,namchk)
      end if
c
c     has printing been turned off externally?
c
      call iosys('read character "print flag" from rwf',-1,0,0,prtflg)
c
      if(prtflg.eq.'minimum') prnt=.false.
c
c     ----- get the dimensions, etc -----
c
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbf,0,' ')
      nnp=(nbf+1)*nbf/2
      call iosys('read integer "spin multiplicity" from rwf',
     $     1,multip,0,' ')
      call iosys('read integer "number of atoms" from rwf',
     $     1,natoms,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $     1,nae,0,' ')
      call iosys('read integer "number of beta electrons" from rwf',
     $     1,nbe,0,' ')
c
      if (nbf.gt.maxnbf) then
         call lnkerr('character core not long enough')
      end if
c
c     ----- set up some parameters depending on multip -----
c
c
      gscf=logkey(ops,'scf=gscf',.false.,' ')
c
      if (multip.eq.1 .and. (.not.gscf)) then
         if (logkey(ops,'hf=tcscf',.false.,' ')) then
c
c           ----- two-configuration scf (gvb) -----
c
            scfnum=2
            nfock=3
            ncoul=3
            nexch=3
            nshell=4
            ndmat=3
            if (logkey(ops,'scf=pulay',.false.,' ')) then
               scfnum=3
               nfock=4
            end if
         else
c
c           ----- closed-shell rhf -----
c
            scfnum=0
            nfock=1
            ncoul=1
            nexch=1
            nshell=2
            ndmat=1
            if (logkey(ops,'scf=pulay',.false.,' ')) scfnum=5
         end if
      else if (logkey(ops,'scf=pulay',.false.,' ') .and.
     #          nbe.ne.0) then
         scfnum=4
         nfock=2
         ncoul=2
         nexch=2
         nshell=3
         ndmat=2
      else if(gscf.or.nbe.eq.0) then
c
c        general scf
c
c
         scfnum=6
c
         if(nbe.ne.0) then
            nshell=2
         else
            nshell=1
         end if
c
         nshell=intkey(ops,'scf=nshell',nshell,' ')
         nfock=intkey(ops,'scf=nfock',nshell,' ')
         ncoul=intkey(ops,'scf=ncoul',nshell,' ')
         nexch=intkey(ops,'scf=nexch',nshell,' ')
         ndmat=intkey(ops,'scf=ndmat',nshell,' ')
c
c        nshell is the total number of shells including the virtuals
c        input value is the number of occupied shells
c
         nshell=nshell+1
c
      else
         scfnum=1
         nfock=2
         ncoul=2
         nexch=2
         nshell=3
         ndmat=2
      end if
c
      call iosys('write integer "number of shells" to rwf',1,nshell,
     $     0,' ')
      call iosys('write integer "number of hf density matrices" '//
     $     'to rwf',1,ndmat,0,' ')
c
c     ----- check on options and set defaults -----
c
      mxiter=intkey(ops,'scf=cycles',30,' ')
      pseudo=intkey(ops,'scf=pseudocanonical',1,' ')
      extrap=intkey(ops,'scf=extrapolate',0,' ')
      diagnl=intkey(ops,'scf=diagonals',0,' ')
      stdiis=intkey(ops,'scf=startdiis',-4,' ')
      finish=intkey(ops,'scf=convergence',8,' ')
      level=fpkey(ops,'scf=level-shift',0.0d+00,' ')
      level1=fpkey(ops,'scf=level-shift1',0.0d+00,' ')
      if (level.ne.0.0d+00) level1=0.0d+00
      mxds=intkey(ops,'scf=maxdiis',5,' ')
c
c     ----- allocate core -----
c
c
      ptr=1
      dsptr=ptr+mxiter
      shlnbf=dsptr+mxds
      shlmin=shlnbf+nshell
      shlmax=shlmin+nshell
      fcoef=iadtwp(shlmax+nshell)
      alpha=fcoef+nshell
      beta=alpha+nshell**2
      need1=wpadti(beta+nshell*nshell)
      call getmem(need1,p1,ngot(1),'m501',0)
      s=1
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
c..bhl
      csym=dlast+nnp
      t1=csym+nbf**2
c..bhl
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
c      call getscm(0,z,canget,'m501: how much core',0)
c
      wish=wpadti(values+nnp**2)
      call getmem(wish,p2,ngot(2),'m501',1)
      if (wish.le.ngot(2)) then
         ntriang=nnp
      else
         left=min(memtri,iadtwp(ngot(2))-values-nnp)
         ntriang=left/nnp+1
      end if
      need2=wpadti(values+nnp*ntriang)
 
c     ----- get the core needed -----
c
c      call getscm(need,z,maxcor,'m501',1)
c
      lenbuf=ntriang*nnp
c
      if(prnt) then
         write(iout,2)
         write(iout,3) ngot(2)
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
      if (ntriang.lt.1) call lnkerr('not enough core in m501')
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
     $        'to ints',need2,a2,0,' ')
         call iosys('close oldints',0,0,0,' ')
      else
         call iosys('open ints as old',0,0,0,namint)
      end if
c
c     ----- read in s, t and v one-electron integrals -----
c
      call iosys('read real "overlap integrals" from rwf',
     $           nnp,z2(s),0,' ')
      call iosys('read real "kinetic integrals" from rwf',
     $           nnp,z2(t),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $           nnp,z2(v),0,' ')
c
c     ----- form s**(1/2) and s**(-1/2) -----
c
      call sinv(z2(s),z2(smhalf),z2(u),z2(eigval),z2(t1),z2(t2),
     $          nbf,nnp,z2(triang),iprint)
c
c     ----- initialize the arrays -----
c
      call setup(a1(shlnbf),a1(shlmin),a1(shlmax),nshell,nbf,z1(alpha),
     $     z1(beta),z1(fcoef),scfnum,ops,z2(t1),z2(t2),nocc)
c
c     ----- do the scf iterations -----
c
      call scf(z2(s),z2(t),z2(v),z2(h),z2(f),z2(u),z2(eigval),
     $         z2(t1),z2(t2),z2(t3),z2(triang),z2(values),
     $         z2(c),z2(smhalf),nbf,nnp,lenbuf,itap44,z2(d),z2(t4),
     $         z2(t5),z2(error),z2(focksv),mxiter,finish,stdiis,
     $         nocc,iprint,z2(dlast),z2(flast),ntriang,z2(t1),
     $         nshell,a1(shlmin),a1(shlmax),z1(alpha),z1(beta),
     $         scfnum,ncoul,nexch,ndmat,nfock,pseudo,extrap,diagnl,
     $         z2(energs),z2(diissv),prnt,ops,mxds,a1(ptr),a1(dsptr),
     $         level,level1,z2(cdiis),z1(fcoef),z2(jmat),z2(kmat),
     $         bflabl,z2(csym))
      call getmem(-ngot(2),p2,idum,'m501',idum)
c
c     ----- form the lagrangian -----
c
c
c     ----- set up some parameters depending on multip -----
c
      if (scfnum.eq.0.or.scfnum.eq.5) then
         nfock=1
         ncoul=1
         nexch=1
         nshell=2
         ndmat=1
      else if (scfnum.eq.1.or.scfnum.eq.4) then
         nfock=2
         ncoul=2
         nexch=2
         nshell=3
         ndmat=2
      else if (scfnum.eq.2.or.scfnum.eq.3) then
         nfock=3
         ncoul=3
         nexch=3
         nshell=4
         ndmat=3
      end if
c
c     ----- allocate core -----
c
      h=1
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
c      call getscm(0,z,canget,'m501: how much core',0)
c
      wish=wpadti(values+nnp**2)
      call getmem(wish,p2,ngot(2),'m501',0)
      if (wish.le.ngot(2)) then
         ntriang=nnp
      else
         left=min(memtri,iadtwp(ngot(2))-values-nnp)
         ntriang=left/nnp+1
      end if
c      need=wpadti(values+nnp*ntriang)
 
c     ----- get the core needed -----
c
c      call getscm(need,z,maxcor,'m501',1)
c
      lenbuf=ntriang*nnp
c
      if (ntriang.lt.1) call lnkerr('not enough core in m501')
c
c     ----- read in t and v one-electron integrals -----
c
      call iosys('read real "kinetic integrals" from rwf',
     $           nnp,z2(values),0,' ')
      call iosys('read real "potential integrals" from rwf',
     $            nnp,z2(h),0,' ')
      call vadd(z2(h),z2(h),z2(values),nnp)
c
c     ----- form the lagrangian -----
c
      call lagrng(a1(shlnbf),a1(shlmin),a1(shlmax),z1(fcoef),
     $            z1(alpha),z1(beta),
     $            z2(h),z2(c),z2(d),z2(j),z2(k),z2(values),nbf,
     $            nnp,nshell,ndmat,ncoul,nexch,scfnum,ntriang,
     $            z2(t1),z2(t2),z2(lag),ops)
      call getmem(-ngot(2),p2,idum,'m501',idum)
      call getmem(-ngot(1),p1,idum,'m501',idum)
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
