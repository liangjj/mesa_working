*deck trap3d.f 
c***begin prologue     trap3d
c***date written       970810   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bec, trap
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             trap3d
c***purpose            three dimensional bec trap
c
c***description        the one, two or three dimensional eigenvalue
c***                   problem defined by the gross-pitaevski equation
c***                   is solved using a discrete variable representation. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       trap3d
      program trap3d
c
      implicit integer (a-z)
      parameter ( nrmax=20 )
      character*4096 ops
      character*8 prtflg
      character*1 itoc
      character*80 cpass, title, chrkey, task
      character*1600 card
      character*24 search
      character*128 filbec, filham, becin, becout
      character*2 atom
      character*16 coord, type
      character*8 qtyp, qdtyp, ptype
      logical dollar, logkey, prncof, prnply, prnwpt, prnh 
      logical prnht, prall, chko, check, fix, nfix, parity, prnh0
      logical toau, useau, hamd, incore, cgrid
      logical itdiag, dodiis, immed, prbufh
      logical prdvd, dvdall, prdiis, prntds, restrt, trunc, plot
      logical usetf, dstroy
      real*8 z, fpkey, omegat, amass, natmin, natmax, left, right, scatl 
      real*8 pi, u0, der, todiis
      real*8 thresh, cnverg, hbar, diiscn, mindiis
      real*8 massau, lenau, timau, dum, scale, sfac, norm0
      real*8 nrestr, lftept, rtept
      dimension der(2,3), nmax(3), npt(3)
      dimension lftbc(3), rtbc(3), dim(3), neven(3), nodd(3)
      dimension prncof(3), prnply(3), prnwpt(3), prnh(3), check(3)
      dimension q(3), wt(3), p(3), dp(3), ddp(3)
      dimension pn(3), dpn(3), ddpn(3), eigc(3), wtc(3), work(4)
      dimension eig(3), tr(3), omegat(3), left(3), right(3)
      dimension ham(3), v(3), u(3), tpr(3)
      dimension prdvd(11), qtyp(3), qdtyp(3), ncon(3)
      dimension sfac(3)
      dimension nreg(3), nptsr(nrmax,3), lftept(nrmax,3), rtept(nrmax,3)
      dimension fix(3), nfix(2,nrmax,3)
      pointer(pnt,z(1)),(pnt,ia(1))
      common/io/inp, iout      
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
c          mass in Kilograms length in meters
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.418884d-17 /
      call drum
      call manmem(0,idum,idum,'trap3d',idum)
      write(iout,*)
      write(iout,*) '    three dimensional mean-field code          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c                       set print options
c
      prall=logkey(ops,'print=m7000=all',.false.,' ')
c----------------------------------------------------------------------c
c   1. if task = diagonalize-h0, the code will diagonalize a hamiltonian
c      consisting of the kinetic energy plus one and two body local 
c      interactions.  the interaction terms are not assumed to depend on
c      the wavefunction.  parity ( see below ) will work if the potential
c      is an even function of the co-ordinates.
c
c   2. if task = diagonalize-gross-pitaevski, the code will diagonalize
c      a hamiltonian as given above plus a local, non-linear term depending
c      on the square of the wavefunction.  in order to converge this eigenvalue
c      problem it is necessary to use an accelerator; straight iteration will
c      not converge.  we have implimented the diis procedure of peter pulay
c      to deal with this problem.  it works very well in the cases tried.
c      there is an option to use parity if the interaction potential does
c      not change sign when the coordinate is changed to minus its value.
c      this is useful in reducing the basis set size for the non-linear
c      problem by up to a factor of eight in the three-dimensional case.
c      there is no need to transform the potential interaction because the
c      potential must be invariant to replacing x by -x if this whole parity
c      business is to work.
c
c   3. clearly the code could be modified easily to handle more general
c      cases.
c
c   4. the user has the option to diagonalize the hamiltonian directly or
c      iteratively using the davidson technique. a non-diagonal 
c      preconditioning of the hamiltonian is needed to make this efficient.
c      this preconditioning is based on the ability to efficiently solve
c      a separable hamiltonian and to transform vectors efficiently to this
c      new representation.  for large eigenvalue problems, the iterative 
c      approach is the only feasible technique.
c
c----------------------------------------------------------------------c     
      atom=chrkey(ops,'atom','Cs',' ')
      task=chrkey(ops,'calculation','diagonalize-h0',' ')
      coord=chrkey(ops,'coordinate-system','cartesian',' ')
      parity=logkey(ops,'parity=on',.false.,' ')
      if(parity) then
         ptype=chrkey(ops,'parity-type','even',' ')
      endif
      plot=logkey(ops,'plot=on',.false.,' ')
      if(plot) then
         dstroy=logkey(ops,'destroy-old-plot-file',.false.,' ')
         becin=chrkey(ops,'input-plot-filename','"in plot file"',' ')
         becout=chrkey(ops,'output-plot-filename','"out plot file"',' ')
      endif
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      if(dstroy) then
         call iosys('destroy',0,0,0,filbec)
         call iosys ('open bec as new',0,0,0,filbec)
      else
         call iosys ('open bec as unknown',0,0,0,filbec)      
      endif   
      write(iout,1) coord
      if(task.eq.'diagonalize-h0') then 
         write(iout,2)      
      elseif(task.eq.'diagonalize-gross-pitaevski') then
         write(iout,3) atom      
      endif
      numdim=intkey(ops,'number-of-dimensions',1,' ')
      call iosys('write integer "number of dimensions" to bec',
     1            1,numdim,0,' ')
      restrt=logkey(ops,'restart',.false.,' ')
      if(prall) then
         prncof(1)=.true.
         prncof(2)=.true.
         prncof(3)=.true.
         prnply(1)=.true.
         prnply(2)=.true.
         prnply(3)=.true.
         prnwpt(1)=.true. 
         prnwpt(2)=.true.
         prnwpt(3)=.true.
         prnh0=.true. 
         prnht=.true. 
         prnh(1)=.true.
         prnh(2)=.true.
         prnh(3)=.true.
      else         
         prncof(1)=logkey(ops,'print=m7000=q1-polynomial-coefficients',
     1                    .false.,' ')
         prncof(2)=logkey(ops,'print=m7000=q2-polynomial-coefficients',
     1                    .false.,' ')
         prncof(3)=logkey(ops,'print=m7000=q3-polynomial-coefficients',
     1                    .false.,' ')
         prnply(1)=logkey(ops,'print=m7000=q1-polynomials',.false.,' ')
         prnply(2)=logkey(ops,'print=m7000=q2-polynomials',.false.,' ')
         prnply(3)=logkey(ops,'print=m7000=q3-polynomials',.false.,' ')
         prnwpt(1)=logkey(ops,'print=m7000=q1-points/weights',
     1                    .false.,' ')
         prnwpt(2)=logkey(ops,'print=m7000=q2-points/weights',
     1                    .false.,' ')
         prnwpt(3)=logkey(ops,'print=m7000=q3-points/weights',
     1                    .false.,' ')
         prnh0=logkey(ops,'print=m7000=h0',.false.,' ')
         prnht=logkey(ops,'print=m7000=h',.false.,' ')
         prnh(1)=logkey(ops,'print=m7000=q1-hamiltonian',.false.,' ')
         prnh(2)=logkey(ops,'print=m7000=q2-hamiltonian',.false.,' ')
         prnh(3)=logkey(ops,'print=m7000=q3-hamiltonian',.false.,' ')
      endif
      chko=logkey(ops,'m7000=check-orthogonality',.false.,' ')
      if(chko) then
         check(1)=.true.
         check(2)=.true.
         check(3)=.true.
      else                  
         check(1)=logkey(ops,'check-q1-orthogonality',.false.,' ')
         check(2)=logkey(ops,'check-q2-orthogonality',.false.,' ')
         check(3)=logkey(ops,'check-q3-orthogonality',.false.,' ')
      endif   
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
c     incosp=logkey(ops,'in-coordinate-space',.false.,' ')
      if(numdim.eq.1) then
         write(iout,4) prncof(1),
     1                 prnply(1),
     2                 prnwpt(1), prnh0,
     3                 prnh(1), 
     4                 check(1)
      endif
      if(numdim.eq.2) then
         write(iout,5) prncof(1), prncof(2),
     1                 prnply(1), prnply(2),
     2                 prnwpt(1), prnwpt(2), prnh0,
     3                 prnh(1), prnh(2), 
     4                 check(1), check(2)
      endif     
      if(numdim.eq.3) then
         write(iout,6) prncof(1), prncof(2), prncof(3),
     1                 prnply(1), prnply(2), prnply(3),
     2                 prnwpt(1), prnwpt(2), prnwpt(3), prnh0,
     3                 prnh(1), prnh(2), prnh(3), 
     4                 check(1), check(2), check(3)
      endif                     
c
c                set diis options
c
      if(task.eq.'diagonalize-gross-pitaevski') then     
         dodiis=.false.
         if( dollar('$diis',card,cpass,inp) ) then
             dodiis=.true.
             ntrials=1
             diiscy=intkey(card,'number-of-diis-cycles',1,' ')
             diisit=intkey(card,'number-of-diis-iterations-per-cycle',
     1                           10,' ')
             immed=logkey(card,'diis-switch=off',.false.,' ')
             trunc=logkey(card,'truncate',.false.,' ')
             diisze=intkey(card,'truncation-size-of-diis-space',
     1                           diisit,' ')
             todiis=1.d0
             if(immed) then
                todiis=1.d+10
             endif
             diiscn=fpkey(card,'diis-convergence-tight',1.d-08,' ')
             mindiis=fpkey(card,'diis-convergence-loose',1.d-08,' ')
             prdiis=logkey(card,'print=diis',.false.,' ')
             prntds=logkey(card,'print=iteration-information',
     1                     .false.,' ')
             write(iout,7) diiscy, diisit, todiis, diiscn, trunc, diisze
          endif
          if( dollar('$iteration',card,cpass,inp) ) then
             ntrials=1
             diiscy=intkey(card,'number-of-iteration-cycles',1,' ')
             diisit=intkey(card,'number-of-iterations-per-cycle',
     1                           10,' ')
             immed=.false.
             todiis=1.d0
             diiscn=fpkey(card,'iteration-convergence',1.d-08,' ')
             write(iout,17) diisit, diiscn
          endif
      endif    
c
c              set various program options
c      
c
c         set trap configuration, numbers of atoms, various constants
c         and other physical parameters
c       
      if ( dollar('$trap',card,cpass,inp) ) then
           natmin=fpkey(card,'minimum-number-of-atoms',0.d+00,' ')
           natmax=fpkey(card,'maximum-number-of-atoms',1.d+05,' ')
           nstep=intkey(card,'number-of-steps',10,' ')
           call iosys('write real "minimum number of atoms" to bec',
     1                 1,natmin,0,' ')
           call iosys('write real "maximum number of atoms" to bec',
     1                 1,natmax,0,' ')
           call iosys('write integer "number of steps" to bec',
     1                 1,nstep,0,' ')
           call atmdat(atom,amass,scatl,u0,omegat,natmax,
     1                 sfac,type)
           usetf=logkey(card,'use-thomas-fermi',.false.,' ')
      endif
c
c               spatial basis set information
c       
      if ( dollar('$q1-polynomials',card,cpass,inp) ) then
           qtyp(1)=chrkey(card,'coordinate-type','x',' ')
           qdtyp(1)=chrkey(card,'quadrature-type','legendre',' ')
           qdtyp(1)='legendre'
           nreg(1)=intkey(card,'number-of-regions',1,' ')
           call intarr(card,'number-of-points-per-region',nptsr(1,1),
     1                 nreg(1),' ')
           call totpt(npt(1),nptsr(1,1),nreg(1),1,.false.) 
           call fparr(card,'left-end-point-for-each-region',
     1                lftept(1,1),nreg(1),' ')
           call fparr(card,'right-end-point-for-each-region',
     1                rtept(1,1),nreg(1),' ')
           nmax(1)=npt(1)
           fix(1)=logkey(card,'fix-end-points',.false.,' ')
           nr=1
           left(1)=lftept(nr,1)
           nr=nreg(1)
           np=nptsr(nr,1)
           right(1)=rtept(nr,1)
           der(1,1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2,1)=fpkey(card,'right-derivative',0.d0,' ')
           call truth(nfix(1,1,1),.false.,nr)
           if (fix(1)) then
               call fxpt(nfix(1,1,1),nr,card)               
           endif
           lftbc(1)=intkey(card,'left-boundary-condition',0,' ')
           rtbc(1)=intkey(card,'right-boundary-condition',0,' ')
           if(lftbc(1).eq.0) then
              nmax(1)=nmax(1)-1
           endif
           if(rtbc(1).eq.0) then
              nmax(1)=nmax(1)-1
           endif
           ncon(1)=nmax(1)
      endif
      if(numdim.gt.1) then
         if ( dollar('$q2-polynomials',card,cpass,inp) ) then
              qtyp(2)=chrkey(card,'coordinate-type','y',' ')
              qdtyp(2)=chrkey(card,'quadrature-type','legendre',' ')
              qdtyp(2)='legendre'
              nreg(2)=intkey(card,'number-of-regions',1,' ')
              call intarr(card,'number-of-points-per-region',nptsr(1,2),
     1                    nreg(2),' ')
              call totpt(npt(2),nptsr(1,2),nreg(2),2,.false.) 
              call fparr(card,'left-end-point-for-each-region',
     1                   lftept(1,2),nreg(2),' ')
              call fparr(card,'right-end-point-for-each-region',
     1                   rtept(1,2),nreg(2),' ')
              nmax(2)=npt(2)
              fix(2)=logkey(card,'fix-end-points',.false.,' ')
              nr=1
              left(2)=lftept(nr,2)
              nr=nreg(2)
              np=nptsr(nr,2)
              right(2)=rtept(nr,2)
              der(1,2)=fpkey(card,'left-derivative',0.d0,' ')
              der(2,2)=fpkey(card,'right-derivative',0.d0,' ')
              call truth(nfix(1,1,2),.false.,nr)
              if (fix(2)) then
                  call fxpt(nfix(1,1,2),nr,card)               
              endif
              lftbc(2)=intkey(card,'left-boundary-condition',0,' ')
              rtbc(2)=intkey(card,'right-boundary-condition',0,' ')
              if(lftbc(2).eq.0) then
                 nmax(2)=nmax(2)-1
              endif
              if(rtbc(2).eq.0) then
                 nmax(2)=nmax(2)-1
              endif
              ncon(2)=nmax(2)
         endif
      endif
      if(numdim.gt.2) then
         if ( dollar('$q3-polynomials',card,cpass,inp) ) then
              qtyp(3)=chrkey(card,'coordinate-type','z',' ')
              qdtyp(3)=chrkey(card,'quadrature-type','legendre',' ')
              qdtyp(3)='legendre'
              nreg(3)=intkey(card,'number-of-regions',1,' ')
              call intarr(card,'number-of-points-per-region',nptsr(1,3),
     1                    nreg(3),' ')
              call totpt(npt(3),nptsr(1,3),nreg(3),3,.false.) 
              call fparr(card,'left-end-point-for-each-region',
     1                   lftept(1,3),nreg(3),' ')
              call fparr(card,'right-end-point-for-each-region',
     1                   rtept(1,3),nreg(3),' ')
              nmax(3)=npt(3)
              fix(3)=logkey(card,'fix-end-points',.false.,' ')
              nr=1
              left(3)=lftept(nr,3)
              nr=nreg(3)
              np=nptsr(nr,3)
              right(3)=rtept(nr,3)
              der(1,3)=fpkey(card,'left-derivative',0.d0,' ')
              der(2,3)=fpkey(card,'right-derivative',0.d0,' ')
              call truth(nfix(1,1,3),.false.,nr)
              if (fix(3)) then
                  call fxpt(nfix(1,1,3),nr,card)               
              endif
              lftbc(3)=intkey(card,'left-boundary-condition',0,' ')
              rtbc(3)=intkey(card,'right-boundary-condition',0,' ')
              if(rtbc(3).eq.0) then
                 nmax(3)=nmax(3)-1
              endif
              ncon(3)=nmax(3)
         endif
      endif
      if(toau) then
c
c        converts everything to atomic units
c
         write(iout,*) 'converting to atomic units'
         hbar=1.d0
         scatl=scatl/lenau
         do 10 i=1,numdim
            left(i)=left(i)/lenau
            right(i)=right(i)/lenau
            omegat(i)=omegat(i)*timau
 10      continue            
         amass=amass/massau        
      endif
      if(useau) then
c
c        if this option is used then all is entered in atomic units
c
         write(iout,*) 'assuming atomic units'
         hbar=1.d0
         amass=1.d0
         do 20 i=1,numdim
            omegat(i)=1.d0
 20      continue            
      endif
      if(numdim.ge.1) then
         title='polynomial information for first dimension'
         write(iout,12) title, nmax(1), npt(1),
     1                  left(1), right(1), 
     2                  der(1,1), der(2,1) , lftbc(1), 
     3                  rtbc(1)
         title='fixed point information for first dimension'
         write(iout,13) title
         nr=nreg(1)
         do 30 i=1,nr
            if(nfix(1,i,1)) then
               write(iout,14) i, lftept(i,1)
            endif
            if(nfix(2,i,1)) then
               write(iout,14) i, rtept(i,1)
            endif
 30      continue   
      endif
      if(numdim.ge.2) then
         title='polynomial information for second dimension'
         write(iout,12) title, nmax(2), npt(2),
     1                  left(2), right(2), 
     2                  der(1,2), der(2,2) , lftbc(2), 
     3                  rtbc(2)
         title='fixed point information for second dimension'
         write(iout,13) title
         nr=nreg(2)
         do 40 i=1,nr
            if(nfix(1,i,2)) then
               write(iout,14) i, lftept(i,2)
            endif
            if(nfix(2,i,2)) then
               write(iout,14) i, rtept(i,2)
            endif
 40      continue   
      endif
      if(numdim.eq.3) then
         title='polynomial information for third dimension'
         write(iout,12) title, nmax(3), npt(3),
     1                  left(3), right(3), 
     2                  der(1,3), der(2,3) , lftbc(3), 
     3                  rtbc(3)
         title='fixed point information for third dimension'
         write(iout,13) title
         nr=nreg(3)
         do 50 i=1,nr
            if(nfix(1,i,3)) then
               write(iout,14) i, lftept(i,3)
            endif
            if(nfix(2,i,3)) then
               write(iout,14) i, rtept(i,3)
            endif
 50      continue   
      endif  
      if(numdim.eq.1) then
         write(iout,8) amass, omegat(1), natmin, natmax, 
     1                 scatl
      elseif(numdim.eq.2) then
         write(iout,9) amass, (omegat(i),i=1,2), natmin, natmax, 
     1                 scatl
      elseif(numdim.eq.3) then
         write(iout,11) amass, (omegat(i),i=1,3), natmin, natmax, 
     1                 scatl
      endif
      if(parity) then
         ncon(1)=nmax(1)/2 + 1
         ncon(2)=nmax(2)/2 + 1
         ncon(3)=nmax(3)/2 + 1
      else
         ncon(1)=nmax(1)
         ncon(2)=nmax(2)
         ncon(3)=nmax(3)   
      endif     
      dim(1)=max(nmax(1),npt(1))
      dim(2)=max(nmax(2),npt(2))
      dim(3)=max(nmax(3),npt(3))
      maxd=max(dim(1),dim(2),dim(3))
      call iosys('write integer maxd to bec',1,maxd,0,' ')
      n3d=1
      do 60 i=1,numdim
         n3d=n3d*ncon(i)
 60   continue 
      nroots=n3d         
c
c             set diagonalization procedures
c      
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      if(itdiag) then
         if(dollar('$davidson',card,cpass,inp)) then
            nroots=intkey(card,'number-of-roots',1,' ')
            nroots=min(nroots,n3d)
            ntrials=intkey(card,'number-of-trial-vectors',nroots,' ')
            cnverg=fpkey(card,'convergence',1.d-08,' ')
            thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
            niter=intkey(card,'maximum-number-of-iterations',
     1                   n3d,' ')
            lenbuf=intkey(card,'hamiltonian-buffer',
     1                    min(1000000,n3d*n3d),' ')
            cgrid=logkey(card,'use-coarse-grid-vectors',.false.,' ')
            prbufh=logkey(card,'print=buffered-hamiltonian',.false.,' ') 
            hamd=logkey(card,'to-h0',.false.,' ')
            n0=intkey(card,'to-h0=nkeep',n3d,' ')
            prdvd(1)=logkey(card,'print=davidson=trials',.false.,' ')
            prdvd(2)=logkey(card,'print=davidson=hamiltonian',
     1                      .false.,' ')
            prdvd(3)=logkey(card,'print=davidson=overlaps',
     1                      .false.,' ')
            prdvd(4)=logkey(card,'print=davidson=vectors',
     1                      .false.,' ')
            prdvd(5)=logkey(card,'print=davidson=h-on-vectors',
     1                      .false.,' ')
            prdvd(6)=logkey(card,'print=davidson=small-matrix',
     1                      .false.,' ')
            prdvd(7)=logkey(card,'print=davidson=transformed-vectors',
     1                     .false.,' ')
            prdvd(8)=logkey(card,'print=davidson=transformed-h-'//
     1                      'on-vectors',
     1                      .false.,' ')
            prdvd(9)=logkey(card,'print=davidson=residuals',
     1                      .false.,' ')
            prdvd(10)=logkey(card,'print=davidson=new-trial-vectors',
     1                      .false.,' ')
            prdvd(11)=logkey(card,'print=davidson=iteration-'//
     1                            'information',.false.,' ')
            dvdall=logkey(card,'print=davidson=all',.false.,' ')
            if(dvdall) then
               do 70 i=1,10
                  prdvd(i)=.true.
 70            continue
            endif
         endif             
         call iosys ('read character "hamiltonian filename" from rwf',
     1                -1,0,0,filham)
         call iosys('open ham as new',0,0,0,filham)
         call iosys('write integer "buffer size" to ham',1,
     1               lenbuf,0,' ')           
         call iosys('write integer "davidson print options" to ham',
     1               11,prdvd,0,' ')
         write(iout,15) nroots, thresh, cnverg, niter
      else
         nroots=intkey(ops,'number-of-roots',n3d,' ')
         nroots=min(nroots,n3d)
         write(iout,16) n3d, nroots
      endif
      call iosys('write integer "size of basis" to bec',1,
     1            n3d,0,' ')
      ioff=1
c
c     note the alignment of variables.  it is important in passing these 
c     arrays later.
c 
      if(numdim.eq.1) then
         q(1)=ioff
         wt(1)=q(1)+dim(1)
         p(1)=wt(1)+dim(1)
         dp(1)=p(1)+dim(1)*dim(1)
         ddp(1)=dp(1)+dim(1)*dim(1)
         pn(1)=ddp(1)+dim(1)*dim(1)
         dpn(1)=pn(1)+dim(1)*dim(1)
         ddpn(1)=dpn(1)+dim(1)*dim(1)
         eigc(1)=ddpn(1)+dim(1)*dim(1)
         wtc(1)=eigc(1)+dim(1)
         tr(1)=wtc(1)+dim(1)
         ham(1)=tr(1)+dim(1)*dim(1)
         eig(1)=ham(1)+dim(1)*dim(1)
         v(1)=eig(1)+dim(1)
         u(1)=v(1)+dim(1)*dim(1)
         words=u(1)
         if(parity) then
            tpr(1)=words+dim(1)*dim(1)
            words=tpr(1)
         endif           
         ind=wpadti(words+dim(1)*dim(1))
         idum=ind+2*n3d
         words=iadtwp(idum+2*n3d)            
      elseif(numdim.eq.2) then 
         q(1)=ioff
         q(2)=q(1)+dim(1)
         wt(1)=q(2)+dim(2)
         wt(2)=wt(1)+dim(1)          
         p(1)=wt(2)+npt(2)
         p(2)=p(1)+dim(1)*dim(1)            
         dp(1)=p(2)+dim(2)*dim(2)
         dp(2)=dp(1)+dim(1)*dim(1)
         ddp(1)=dp(2)+dim(2)*dim(2)
         ddp(2)=ddp(1)+dim(1)*dim(1)
         pn(1)=ddp(2)+dim(2)*dim(2)
         pn(2)=pn(1)+dim(1)*dim(1)
         dpn(1)=pn(2)+dim(2)*dim(2)
         dpn(2)=dpn(1)+dim(1)*dim(1)
         ddpn(1)=dpn(2)+dim(2)*dim(2)
         ddpn(2)=ddpn(1)+dim(1)*dim(1)
         eigc(1)=ddpn(2)+dim(2)*dim(2)
         eigc(2)=eigc(1)+dim(1)
         wtc(1)=eigc(2)+dim(2)
         wtc(2)=wtc(1)+dim(1)
         tr(1)=wtc(2)+dim(2)
         tr(2)=tr(1)+dim(1)*dim(1)
         ham(1)=tr(2)+dim(2)*dim(2)
         ham(2)=ham(1)+dim(1)*dim(1)
         eig(1)=ham(2)+dim(2)*dim(2) 
         eig(2)=eig(1)+dim(1)
         v(1)=eig(2)+dim(2) 
         v(2)=v(1)+dim(1)*dim(1) 
         u(1)=v(2)+dim(2)*dim(2)
         u(2)=u(1)+dim(1)*dim(1)
         words=u(2)
         if(parity) then
            tpr(1)=words+dim(2)*dim(2)
            tpr(2)=tpr(1)+dim(1)*dim(1)
            words=tpr(2)
         endif    
         ind=wpadti(words+dim(2)*dim(2))
         idum=ind+3*n3d
         words=iadtwp(idum+3*n3d)
      elseif(numdim.eq.3) then
         q(1)=ioff
         q(2)=q(1)+dim(1)
         q(3)=q(2)+dim(2)
         wt(1)=q(3)+dim(3)
         wt(2)=wt(1)+dim(1)          
         wt(3)=wt(2)+dim(2)          
         p(1)=wt(3)+dim(3)+1
         p(2)=p(1)+dim(1)*dim(1)            
         p(3)=p(2)+dim(2)*dim(2)            
         dp(1)=p(3)+dim(3)*dim(3)
         dp(2)=dp(1)+dim(1)*dim(1)
         dp(3)=dp(2)+dim(2)*dim(2)
         ddp(1)=dp(3)+dim(3)*dim(3)
         ddp(2)=ddp(1)+dim(1)*dim(1)
         ddp(3)=ddp(2)+dim(2)*dim(2)
         pn(1)=ddp(3)+dim(3)*dim(3)
         pn(2)=pn(1)+dim(1)*dim(1)
         pn(3)=pn(2)+dim(2)*dim(2)
         dpn(1)=pn(3)+dim(3)*dim(3)
         dpn(2)=dpn(1)+dim(1)*dim(1)
         dpn(3)=dpn(2)+dim(2)*dim(2)
         ddpn(1)=dpn(3)+dim(3)*dim(3)
         ddpn(2)=ddpn(1)+dim(1)*dim(1)
         ddpn(3)=ddpn(2)+dim(2)*dim(2)
         eigc(1)=ddpn(3)+dim(3)*dim(3)
         eigc(2)=eigc(1)+dim(1)
         eigc(3)=eigc(2)+dim(2)
         wtc(1)=eigc(3)+dim(3)
         wtc(2)=wtc(1)+dim(1)
         wtc(3)=wtc(2)+dim(2)
         tr(1)=wtc(3)+dim(3)
         tr(2)=tr(1)+dim(1)*dim(1)
         tr(3)=tr(2)+dim(2)*dim(2)
         ham(1)=tr(3)+dim(3)*dim(3)
         ham(2)=ham(1)+dim(1)*dim(1)
         ham(3)=ham(2)+dim(2)*dim(2)
         eig(1)=ham(3)+dim(3)*dim(3) 
         eig(2)=eig(1)+dim(1) 
         eig(3)=eig(2)+dim(2) 
         v(1)=eig(3)+dim(3)
         v(2)=v(1)+dim(1)*dim(1)
         v(3)=v(2)+dim(2)*dim(2)
         u(1)=v(3)+dim(3)*dim(3)
         u(2)=u(1)+dim(1)*dim(1)
         u(3)=u(2)+dim(2)*dim(2)
         words=u(3)
         if(parity) then
            tpr(1)=words+dim(3)*dim(3)
            tpr(2)=tpr(1)+dim(1)*dim(1)
            tpr(3)=tpr(2)+dim(2)*dim(2)
            words=tpr(3)
         endif                  
         ind=wpadti(words+dim(3)*dim(3))
         idum=ind+4*n3d
         words=iadtwp(idum+4*n3d)
      endif
      space=words    
      work(1)=words
      work(2)=work(1)+max(2*maxd*maxd,3*n3d,3*maxd)
      words=work(2)+2*maxd*maxd
      if(task.eq.'diagonalize-h0') then
         if(.not.itdiag) then
             hamtot=words
             eigtot=hamtot+n3d*n3d
             vtot=eigtot+n3d
             words=vtot+n3d
         else
             hbuf=words
             ibuf=wpadti(hbuf+lenbuf)
             diag=iadtwp(ibuf+2*lenbuf)
             vtot=diag+n3d
             etrial=vtot+n3d
             trials=etrial+n3d
             pvec=trials+n3d*ntrials
             hpvec=pvec+n3d*niter
             vec=hpvec+niter*n3d
             bmat=vec+n3d*niter
             bmatm=bmat+niter*niter
             eigtot=bmatm+niter*niter
             work(3)=eigtot+n3d
             work(4)=work(3)+n3d
             words=work(4)+n3d
         endif
      elseif(task.eq.'diagonalize-gross-pitaevski') then
         psi=words
         vnl=psi+n3d*(diisit+1)
         vnlcur=vnl+n3d*(diisit+1)
         fpsi=vnlcur+n3d*(diisit+1)
         bb=fpsi+n3d*(diisit+1)
         bbtmp=bb+(diisit+1)*(diisit+1)
         sol=bbtmp+(diisit+1)*(diisit+1)
         hh=sol+max(diisit+1,n3d)
         pp=hh+diisit*diisit
         ph=pp+diisit*diisit
         hp=ph+diisit*diisit
         etrial=hp+diisit*diisit
         trials=etrial+n3d
         eigtot=trials+n3d*ntrials
         vtot=eigtot+n3d
         hbuf=vtot+n3d
         fmat=hbuf       
         ipvt=idum
         ibuf=idum
         diag=hbuf
         pvec=hbuf
         hpvec=hbuf
         vec=hbuf
         bmat=hbuf
         bmatm=hbuf
         if(.not.itdiag) then
            fmat=hbuf+n3d*n3d
            ipvt=wpadti(fmat+n3d*n3d)
            words=iadtwp(ipvt+n3d) 
         else
            ibuf=wpadti(hbuf+lenbuf)
            diag=iadtwp(ibuf+2*lenbuf)
            pvec=diag+n3d
            hpvec=pvec+n3d*niter
            vec=hpvec+niter*n3d
            bmat=vec+n3d*niter
            bmatm=bmat+niter*niter
            work(3)=bmat+niter*niter
            work(4)=work(3)+n3d
            words=work(4)+n3d
         endif
      endif
      words=wpadti(words)
      call manmem(words,pnt,ngot,'poly',0)
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      do 90 i=1,numdim
c      
         call getqpt(z(q(i)),z(wt(i)),lftept(1,i),rtept(1,i),
     1               qdtyp(i),'before',z(work(1)),
     2               nfix(1,1,i),nptsr(1,i),npt(i),nreg(i),.false.)
         if(qdtyp(i).eq.'hermite') then
            call smul(z(q(i)),z(q(i)),sfac(i),npt(i))
            call smul(z(wt(i)),z(wt(i)),sfac(i),npt(i))
         endif   
         call lgrply(z(p(i)),z(dp(i)),z(ddp(i)),z(q(i)),
     1               z(q(i)),npt(i),npt(i),.false.)
         call prepfn(z(p(i)),z(dp(i)),z(ddp(i)),z(q(i)),
     1               z(wt(i)),npt(i),sfac(i),qdtyp(i))
         call copy(z(p(i)),z(pn(i)),npt(i)*npt(i))
         call copy(z(dp(i)),z(dpn(i)),npt(i)*npt(i))
         call copy(z(ddp(i)),z(ddpn(i)),npt(i)*npt(i))
         if(prnply(i)) then                  
            title='dvr polynomials '//itoc(i)//' dimension'
            call prntfm(title,z(pn(i)),npt(i),npt(i),npt(i),
     1                  npt(i),iout)
            title='first derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(dpn(i)),npt(i),npt(i),npt(i),
     1                  npt(i),iout)
            title='second derivative of dvr polynomials '
     1             //itoc(i)//' dimension'
            call prntfm(title,z(ddpn(i)),npt(i),npt(i),npt(i),
     1                  npt(i),iout)
         endif
         if(check(i)) then
            call chk(z(pn(i)),z(wt(i)),z(work(1)),z(work(2)),
     1               npt(i),npt(i))
            title='overlap matrix dvr polynomials '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(work(2)),npt(i),npt(i),npt(i),
     1                                npt(i),iout)
         endif
         if(lftbc(i).eq.0.or.rtbc(i).eq.0) then
            call dropfn(z(p(i)),z(dp(i)),z(ddp(i)),
     1                  z(pn(i)),z(dpn(i)),z(ddpn(i)),z(q(i)),
     2                  z(wt(i)),z(eigc(i)),z(wtc(i)),
     3                  lftbc(i),rtbc(i),npt(i),nmax(i))
         endif
         search='look for v0'
         call vmat(z(eigc(i)),dum,dum,z(v(i)),z(work(1)),hbar,amass,
     1             scale,n3d,nmax(1),nmax(2),nmax(3),1,prnh0,
     2             search,ops,parity)
c
c        the hamiltonian matrix contains the kinetic energy and any
c        one body potential.
c
         call ham0(z(pn(i)),z(dpn(i)),z(ddpn(i)),z(eigc(i)),z(wtc(i)),
     1             z(ham(i)),z(v(i)),z(u(i)),z(eig(i)),z(work(1)),
     2             z(tpr(i)),hbar,amass,nmax(i),lftbc(i),rtbc(i),
     3             coord,qtyp(i),prnh0,parity,ptype,neven(i),nodd(i))
         if(parity) then
            if(ptype.eq.'even') then
               ncon(i)=neven(i)
            elseif(ptype.eq.'odd') then
               ncon(i)=nodd(i)
            endif
         else
            ncon(i)=nmax(i)
         endif
         call vscale(z(eig(i)),z(eig(i)),scale,ncon(i)) 
         if(prnh0) then
            title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(eig(i)),ncon(i),1,ncon(i),1,iout)
         endif
         call vscale(z(eig(i)),z(eig(i)),1.d0/scale,ncon(1))
         if (prnh0) then
             title='eigenvectors of unperturbed hamiltonian '//itoc(i)
     1          //' dimension'
             call prntfm(title,z(u(i)),ncon(i),ncon(i),ncon(i),
     1                   ncon(i),iout)
         endif
         call iosys('write integer "number of points for '//
     1              'dimension = '//itoc(i)//'" to bec',1,
     2               nmax(i),0,' ') 
         call iosys('write real "points for dimension = '//itoc(i)
     1              //'" to bec',nmax(i),z(eigc(i)),0,' ')
         call iosys('write real "weights for dimension = '//itoc(i)
     1              //'" to bec',nmax(i),z(wtc(i)),0,' ')
         if (prnply(i)) then
             title='dvr polynomials '//itoc(i)//' dimension'         
             call prntfm(title,z(pn(i)),nmax(i),ncon(i),nmax(i),
     1                   ncon(i),iout)
             title='first derivative of dvr polynomials '//
     1              itoc(i)//' dimension'
             call prntfm(title,z(dpn(i)),nmax(i),ncon(i),nmax(i),
     1                   ncon(i),iout)
             title='second derivative of dvr polynomials '//
     1              itoc(i)//' dimension'
             call prntfm(title,z(ddpn(i)),nmax(i),ncon(i),nmax(i),
     1                   ncon(i),iout)
         endif
 90   continue
      n3d=1
      do 100 i=1,numdim
         n3d=n3d*ncon(i)
 100  continue
      if(.not.itdiag) then 
         nroots=n3d
      endif          
      call setind(ia(ind),ncon,n3d,numdim)
      search='look for v1'
      if(numdim.eq.1) then
         call vmat(z(eigc(1)),dum,dum,z(vtot),z(work(1)),
     1             hbar,amass,scale,n3d,ncon(1),ncon(2),ncon(3),
     2             numdim,prnh0,search,ops,parity)   
      elseif(numdim.eq.2) then
         call vmat(z(eigc(1)),z(eigc(2)),dum,z(vtot),z(work(1)),
     1             hbar,amass,scale,n3d,ncon(1),ncon(2),ncon(3),
     2             numdim,prnh0,search,ops,parity)   
      elseif(numdim.eq.3) then
         call vmat(z(eigc(1)),z(eigc(2)),z(eigc(3)),z(vtot),
     1             z(work(1)),hbar,amass,scale,n3d,
     2             ncon(1),ncon(2),ncon(3),numdim,prnh0,search,
     3             ops,parity)  
      endif
c
c     subtract the one body interations included in the zeroth
c     order hamiltonian.
c
      call vpert(z(vtot),z(v(1)),z(v(2)),z(v(3)),n3d,
     1           ncon(1),ncon(2),ncon(3),numdim)
      if(task.eq.'diagonalize-gross-pitaevski') then
         if(numdim.eq.1) then
            call vtrial(z(ham(1)),dum,dum,
     1                  z(u(1)),dum,dum,
     2                  z(eig(1)),dum,dum,
     3                  z(trials),z(etrial),ia(idum),
     4                  ncon(1),ncon(2),ncon(3),numdim,n3d,
     5                  ntrials,.false.)
            if(itdiag) then
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),dum,dum,z(vtot),z(hbuf),ia(ibuf),
     1                    ia(ind),z(diag),lenbuf,
     2                    ncon(1),ncon(2),ncon(3),numdim,n3d,
     3                    prbufh,nel,incore,title,.true.)
            else
               call lschr(z(ham(1)),dum,dum,z(vtot),z(hbuf),
     1                    dum,ndum,n3d,ncon(1),ncon(2),ncon(3),
     2                    numdim,.true.,.false.,prnht)    
            endif
            norm0=1.d0
            if(coord.eq.'spherical') then
                     norm0=sqrt(1.d0/(4.d0*pi))
            endif     
            call nlschr(z(eigtot),z(psi),z(fpsi),z(hbuf),
     1                  ia(ibuf),z(diag),z(fmat),z(trials),z(etrial),
     2                  z(pvec),z(hpvec),z(vec),
     3                  z(vnl),z(vnlcur),
     4                  z(pn(1)),dum,dum,
     5                  z(eigc(1)),dum,dum,
     6                  z(eig(1)),dum,dum,
     7                  z(u(1)),dum,dum,
     8                  z(bb),z(bbtmp),z(hh),z(pp),z(ph),z(hp),
     9                  z(sol),z(bmat),z(bmatm),
     x                  z(work(3)),z(work(4)),
     x                  ia(ipvt),natmin,natmax,
     x                  norm0,u0,scale,todiis,diiscn,mindiis,cnverg,
     x                  thresh,n3d,diiscy,diisit,niter,nroots,ntrials,
     x                  nstep,numdim,nmax(1),nmax(2),nmax(3),
     x                  ncon(1),ncon(2),ncon(3),idum,incore,prdiis,
     x                  prntds,dodiis,hamd,itdiag,trunc,diisze,
     x                  restrt,nrestr,atom,usetf,plot,becin,becout)
         elseif(numdim.eq.2) then
            norm0=1.d0
            call vtrial(z(ham(1)),z(ham(2)),dum,
     1                  z(u(1)),z(u(2)),dum,
     2                  z(eig(1)),z(eig(2)),dum,
     3                  z(trials),z(etrial),ia(idum),
     4                  ncon(1),ncon(2),ncon(3),numdim,n3d,
     5                  ntrials,.false.)
            if(itdiag) then
               title='buffered hamiltonian for '//itoc(numdim)
     1               //'d dimensional hamiltonian'
               call hamil(z(ham(1)),z(ham(2)),dum,z(vtot),z(hbuf),
     1                    ia(ibuf),ia(ind),z(diag),lenbuf,
     2                    ncon(1),ncon(2),ncon(3),
     3                    numdim,n3d,prbufh,nel,incore,title,
     4                   .true.)
            else                     
               call lschr(z(ham(1)),z(ham(2)),dum,z(vtot),z(hbuf),
     1                    dum,ia(ind),n3d,ncon(1),ncon(2),ncon(3),
     2                    numdim,.true.,.false.,prnht)
            endif
            call nlschr(z(eigtot),z(psi),z(fpsi),z(hbuf),
     1                  ia(ibuf),z(diag),z(fmat),z(trials),z(etrial),
     2                  z(pvec),z(hpvec),z(vec),
     3                  z(vnl),z(vnlcur),
     4                  z(pn(1)),z(pn(2)),dum,
     5                  z(eigc(1)),z(eigc(2)),dum,
     6                  z(eig(1)),z(eig(1)),dum,
     7                  z(u(1)),z(u(2)),dum,
     8                  z(bb),z(bbtmp),z(hh),z(pp),z(ph),z(hp),
     9                  z(sol),z(bmat),z(bmatm),
     x                  z(work(3)),z(work(4)),
     x                  ia(ipvt),natmin,natmax,
     x                  norm0,u0,scale,todiis,diiscn,mindiis,cnverg,
     x                  thresh,n3d,diiscy,diisit,niter,nroots,ntrials,
     x                  nstep,numdim,nmax(1),nmax(2),nmax(3),
     x                  ncon(1),ncon(2),ncon(3),idum,incore,prdiis,
     x                  prntds,dodiis,hamd,itdiag,trunc,diisze,
     x                  restrt,nrestr,atom,usetf,plot,becin,becout)
         elseif(numdim.eq.3) then
            norm0=1.d0
            call vtrial(z(ham(1)),z(ham(2)),z(ham(3)),
     1                  z(u(1)),z(u(2)),z(u(3)),
     2                  z(eig(1)),z(eig(2)),z(eig(3)),
     3                  z(trials),z(etrial),ia(idum),
     4                  ncon(1),ncon(2),ncon(3),numdim,n3d,
     5                  ntrials,.false.)
            if(itdiag) then
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),
     1                    z(hbuf),ia(ibuf),ia(ind),z(diag),
     2                    lenbuf,ncon(1),ncon(2),ncon(3),numdim,
     3                    n3d,prbufh,nel,incore,title,.true.)
            else
               call lschr(z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),
     1                    z(hbuf),dum,ia(ind),n3d,
     2                    ncon(1),ncon(2),ncon(3),
     2                    numdim,.true.,.false.,prnht)
            endif                                             
            call nlschr(z(eigtot),z(psi),z(fpsi),z(hbuf),
     1                  ia(ibuf),z(diag),z(fmat),z(trials),z(etrial),
     2                  z(pvec),z(hpvec),z(vec),
     3                  z(vnl),z(vnlcur),
     4                  z(pn(1)),z(pn(2)),z(pn(3)),
     5                  z(eigc(1)),z(eigc(2)),z(eigc(3)),
     6                  z(eig(1)),z(eig(1)),z(eig(3)),
     7                  z(u(1)),z(u(2)),z(u(3)),
     8                  z(bb),z(bbtmp),z(hh),z(pp),z(ph),z(hp),
     9                  z(sol),z(bmat),z(bmatm),
     x                  z(work(3)),z(work(4)),
     x                  ia(ipvt),natmin,natmax,
     x                  norm0,u0,scale,todiis,diiscn,mindiis,cnverg,
     x                  thresh,n3d,diiscy,diisit,niter,nroots,ntrials,
     x                  nstep,numdim,nmax(1),nmax(2),nmax(3),
     x                  ncon(1),ncon(2),ncon(3),idum,incore,prdiis,
     x                  prntds,dodiis,hamd,itdiag,trunc,diisze,
     x                  restrt,nrestr,atom,usetf,plot,becin,becout)
         endif
      elseif(task.eq.'diagonalize-h0') then
         if(itdiag) then
            call iosys('write integer "length of davidson vector" '//
     1                 'to ham',1,n3d,0,' ')
            if(numdim.eq.1) then
               call vtrial(z(ham(1)),dum,dum,
     1                     z(u(1)),dum,dum,
     2                     z(eig(1)),dum,dum,
     3                     z(trials),z(etrial),ia(idum),
     4                     ncon(1),ncon(2),ncon(3),numdim,n3d,
     5                     ntrials,.false.)
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),dum,dum,z(vtot),z(hbuf),ia(ibuf),
     1                    ia(ind),z(diag),lenbuf,
     2                    ncon(1),ncon(2),ncon(3),numdim,n3d,
     3                    prbufh,nel,incore,title,.true.)
               call iosys('write integer "size of primitive '//
     1                    'basis for dimension 1" to ham',1,
     2                     nmax(1),0,' ')
               call iosys('write integer "size of contracted '//
     1                    'basis for dimension 1" to ham',1,
     2                     ncon(1),0,' ')
               call iosys('write real "lagrange basis functions '//
     1                    'for dimension 1" to ham',ncon(1)*nmax(1),
     2                     z(pn(1)),0,' ')
               call iosys('write real "first derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 1" to ham',ncon(1)*nmax(1),
     3                     z(dpn(1)),0,' ')
               call iosys('write real "second derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 1" to ham',ncon(1)*nmax(1),
     3                     z(ddpn(1)),0,' ')
               call dvd(z(hbuf),ia(ibuf),z(diag),z(trials),z(etrial),
     1                  z(eig(1)),dum,dum,
     2                  z(u(1)),dum,dum,
     3                  z(pvec),z(hpvec),z(vec),
     4                  z(bmat),z(bmatm),z(vec),z(eigtot),
     5                  z(work(1)),z(work(2)),
     6                  scale,cnverg,thresh,ncon(1),ncon(2),ncon(3),
     7                  numdim,n3d,nroots,ntrials,niter,lenbuf,
     8                  hamd,incore,nel,prdvd)
            elseif(numdim.eq.2) then
               call vtrial(z(ham(1)),z(ham(2)),dum,
     1                     z(u(1)),z(u(2)),dum,
     2                     z(eig(1)),z(eig(2)),dum,
     3                     z(trials),z(etrial),ia(idum),
     4                     ncon(1),ncon(2),ncon(3),numdim,n3d,
     5                     ntrials,.false.)            
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),z(ham(2)),dum,z(vtot),z(hbuf),
     1                    ia(ibuf),ia(ind),z(diag),lenbuf,
     2                    ncon(1),ncon(2),ncon(3),numdim,n3d,prbufh,
     3                    nel,incore,title,.true.)
               call iosys('write integer "size of primitive '//
     1                    'basis for dimension 1" to ham',1,
     2                     nmax(1),0,' ')
               call iosys('write integer "size of contracted '//
     1                    'basis for dimension 1" to ham',1,
     2                     ncon(1),0,' ')
               call iosys('write integer "size of primitive '//
     1                    'basis for dimension 2" to ham',1,
     2                     nmax(2),0,' ')
               call iosys('write integer "size of contracted '//
     1                    'basis for dimension 2" to ham',1,
     2                     ncon(2),0,' ')
               call iosys('write real "lagrange basis functions '//
     1                    'for dimension 1" to ham',ncon(1)*nmax(1),
     2                     z(pn(1)),0,' ')
               call iosys('write real "first derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 1" to ham',ncon(1)*nmax(1),
     3                     z(dpn(1)),0,' ')
               call iosys('write real "second derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 1" to ham',ncon(1)*nmax(1),
     3                     z(ddpn(1)),0,' ')
               call iosys('write real "lagrange basis functions '//
     1                    'for dimension 2" to ham',ncon(2)*nmax(2),
     2                     z(pn(2)),0,' ')
               call iosys('write real "first derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 2" to ham',ncon(2)*nmax(2),
     3                     z(dpn(2)),0,' ')
               call iosys('write real "second derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 2" to ham',ncon(2)*nmax(2),
     3                     z(ddpn(2)),0,' ')
               call dvd(z(hbuf),ia(ibuf),z(diag),z(trials),z(etrial),
     1                  z(eig(1)),z(eig(2)),dum,
     2                  z(u(1)),z(u(2)),dum,
     3                  z(pvec),z(hpvec),z(vec),
     4                  z(bmat),z(bmatm),z(vec),z(eigtot),
     5                  z(work(1)),z(work(2)),scale,
     6                  cnverg,thresh,ncon(1),ncon(2),ncon(3),
     7                  numdim,n3d,nroots,ntrials,niter,lenbuf,
     8                  hamd,incore,nel,prdvd)     
            elseif(numdim.eq.3) then
               call vtrial(z(ham(1)),z(ham(2)),z(ham(3)),
     1                     z(u(1)),z(u(2)),z(u(3)),
     2                     z(eig(1)),z(eig(2)),z(eig(3)),
     3                     z(trials),z(etrial),ia(idum),
     4                     ncon(1),ncon(2),ncon(3),numdim,n3d,
     5                     ntrials,.false.)
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),
     1                    z(hbuf),ia(ibuf),ia(ind),z(diag),
     2                    lenbuf,ncon(1),ncon(2),ncon(3),numdim,
     3                    n3d,prbufh,nel,incore,title,.true.)
               call iosys('write integer "size of primitive '//
     1                    'basis for dimension 1" to ham',1,
     2                     nmax(1),0,' ')
               call iosys('write integer "size of contracted '//
     1                    'basis for dimension 1" to ham',1,
     2                     ncon(1),0,' ')
               call iosys('write integer "size of primitive '//
     1                    'basis for dimension 2" to ham',1,
     2                     nmax(2),0,' ')
               call iosys('write integer "size of contracted '//
     1                    'basis for dimension 2" to ham',1,
     2                     ncon(2),0,' ')
               call iosys('write integer "size of primitive '//
     1                    'basis for dimension 3" to ham',1,
     2                     nmax(3),0,' ')
               call iosys('write integer "size of contracted '//
     1                    'basis for dimension 3" to ham',1,
     2                     ncon(3),0,' ')
               call iosys('write real "lagrange basis functions '//
     1                    'for dimension 1" to ham',ncon(1)*nmax(1),
     2                     z(pn(1)),0,' ')
               call iosys('write real "first derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 1" to ham',ncon(1)*nmax(1),
     3                     z(dpn(1)),0,' ')
               call iosys('write real "second derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 1" to ham',ncon(1)*nmax(1),
     3                     z(ddpn(1)),0,' ')
               call iosys('write real "lagrange basis functions '//
     1                    'for dimension 2" to ham',ncon(2)*nmax(2),
     2                     z(pn(2)),0,' ')
               call iosys('write real "first derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 2" to ham',ncon(2)*nmax(2),
     3                     z(dpn(2)),0,' ')
               call iosys('write real "second derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 2" to ham',ncon(2)*nmax(2),
     3                     z(ddpn(2)),0,' ')
               call iosys('write real "lagrange basis functions '//
     1                    'for dimension 3" to ham',ncon(3)*nmax(3),
     2                     z(pn(3)),0,' ')
               call iosys('write real "first derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 3" to ham',ncon(3)*nmax(3),
     3                     z(dpn(3)),0,' ')
               call iosys('write real "second derivative of '//
     1                    'lagrange basis functions for '//
     2                    'dimension 3" to ham',ncon(3)*nmax(3),
     3                     z(ddpn(3)),0,' ')
               call dvd(z(hbuf),ia(ibuf),z(diag),z(trials),z(etrial),
     1                  z(eig(1)),z(eig(2)),z(eig(2)),
     2                  z(u(1)),z(u(2)),z(u(3)),
     3                  z(pvec),z(hpvec),z(vec),
     4                  z(bmat),z(bmatm),z(vec),z(eigtot),
     5                  z(work(1)),z(work(2)),
     6                  scale,cnverg,thresh,ncon(1),ncon(2),ncon(3),
     7                  numdim,n3d,nroots,ntrials,niter,lenbuf,
     8                  hamd,incore,nel,prdvd)     
            endif
         else    
            if(numdim.eq.1) then
               call lschr(z(ham(1)),dum,dum,z(vtot),z(hamtot),z(eigtot),
     1                    ndum,n3d,ncon(1),ncon(2),ncon(3),
     2                    numdim,.true.,.true.,prnht)
            elseif(numdim.eq.2) then
               call lschr(z(ham(1)),z(ham(2)),dum,z(vtot),z(hamtot),
     1                    z(eigtot),ia(ind),n3d,ncon(1),ncon(2),ncon(3),
     2                    numdim,.true.,.true.,prnht)
            elseif(numdim.eq.3) then
               call lschr(z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),
     1                    z(hamtot),z(eigtot),ia(ind),n3d,
     2                    ncon(1),ncon(2),ncon(3),
     2                    numdim,.true.,.true.,prnht)
            endif
            call vscale(z(eigtot),z(eigtot),scale,nroots)
            title='eigenvalues of hamiltonian'
            call prntfm(title,z(eigtot),nroots,1,nroots,1,iout)
         endif
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call manmem(-ngot,pnt,idum,'trap3d',idum)
      call chainx(0)               
      stop
 1    format(/,15x,'coordinate system = ',a16)
 2    format(/,15x,'program options',/,/,5x,
     1             'diagonalize zeroth-order hamiltonian')  
 3    format(/,15x,'program options',/,/,5x,
     1             'diagonalize gross-pitaevski hamiltonian',/,/,5x,
     2             'input data is for the atom = ',a4)
 4    format(/,15x,'code options',/,
     1       /,5x,'print q1 coefficients                  = ',l1,
     2       /,5x,'print q1 polynomials                   = ',l1,
     3       /,5x,'print q1 points/weights                = ',l1,
     4       /,5x,'print spatial hamiltonian information  = ',l1, 
     5       /,5x,'print q1 eigenvectors                  = ',l1, 
     6       /,5x,'check q1 orthonormality                = ',l1)
 5    format(/,15x,'code options',/,
     1       /,5x,'print q1 coefficients                  = ',l1,
     2       /,5x,'print q2 coefficients                  = ',l1,
     3       /,5x,'print q1 polynomials                   = ',l1,
     4       /,5x,'print q2 polynomials                   = ',l1,
     5       /,5x,'print q1 points/weights                = ',l1,
     6       /,5x,'print q2 points/weights                = ',l1,
     7       /,5x,'print spatial hamiltonian information  = ',l1, 
     8       /,5x,'print q1 eigenvectors                  = ',l1, 
     9       /,5x,'print q2 eigenvectors                  = ',l1, 
     x       /,5x,'check q1 orthonormality                = ',l1, 
     x       /,5x,'check q2 orthonormality                = ',l1) 
 6    format(/,15x,'code options',/,
     1       /,5x,'print q1 coefficients                  = ',l1,
     2       /,5x,'print q2 coefficients                  = ',l1,
     3       /,5x,'print q3 coefficients                  = ',l1,
     4       /,5x,'print q1 polynomials                   = ',l1,
     5       /,5x,'print q2 polynomials                   = ',l1,
     6       /,5x,'print q3 polynomials                   = ',l1,
     7       /,5x,'print q1 points/weights                = ',l1,
     8       /,5x,'print q2 points/weights                = ',l1,
     9       /,5x,'print q3 points/weights                = ',l1,
     x       /,5x,'print spatial hamiltonian information  = ',l1, 
     x       /,5x,'print q1 eigenvectors                  = ',l1, 
     x       /,5x,'print q2 eigenvectors                  = ',l1, 
     x       /,5x,'print q3 eigenvectors                  = ',l1, 
     x       /,5x,'check q1 orthonormality                = ',l1, 
     x       /,5x,'check q2 orthonormality                = ',l1, 
     x       /,5x,'check q3 orthonormality                = ',l1) 
 7    format(/,15x,'diis options',/,/,5x,
     1             'number of diis cycles         = ',i4,/,5x,
     2             'maximum number of iterations  = ',i4,/,5x,
     3             'diis begins when error        = ',e15.8,/,5x,
     4             'diis convergence criterion    = ',e15.8,/,5x,
     5             'truncate diis space           = ',l1,/,5x,
     6             'size of diis truncation space = ',i3)
 8    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap q1 frequency                = ',e15.8,/,5x,
     3             'minimum number of atoms in trap  = ',e15.8,/,5x,
     4             'maximum number of atoms in trap  = ',e15.8,/,5x,
     5             'scattering length                = ',e15.8)
 9    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8,/,5x,
     4             'minimum number of atoms in trap   = ',e15.8,/,5x,
     5             'maximum number of atoms in trap   = ',e15.8,/,5x,
     6             'scattering length                 = ',e15.8)
 11   format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                       = ',e15.8,/,5x,
     2             'trap q1 frequency                 = ',e15.8,/,5x,
     3             'trap q2 frequency                 = ',e15.8,/,5x,
     4             'trap q3 frequency                 = ',e15.8,/,5x,
     5             'minimum number of atoms in trap   = ',e15.8,/,5x,
     6             'maximum number of atoms in trap   = ',e15.8,/,5x,
     7             'scattering length                 = ',e15.8)     
 12   format(/,15x,a60,/,/,5x,
     1             'order of polynomials           = ',i3,/,5x,
     2             'number of integration points   = ',i3,/,5x,
     3             'left boundary                  = ',e15.8,/,5x,
     4             'right boundary                 = ',e15.8,/,5x,
     5             'derivative at left boundary    = ',e15.8,/,5x,
     6             'derivative at right boundary   = ',e15.8,/,5x,
     7             'left-boundary-condition        = ',i2,/,5x,
     8             'right-boundary-condition       = ',i2)
 13   format(/,5x,a60)
 14   format(/,5x,
     1             'region        = ',i3,/,5x,
     2             'fixed point   = ',e15.8)
 15   format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'overlap tolerance                  = ',e15.8,/,5x,
     3             'convergence criterion              = ',e15.8,/,5x,
     4             'maximum number of iterations       = ',i6)
 16   format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix  = ',i3,/,5x,
     2             'number of roots = ',i3)
 17   format(/,15x,'iteration options',/,/,5x,
     1             'maximum number of iterations = ',i4,/,5x,
     2             'convergence criterion   = ',e15.8)
      end






