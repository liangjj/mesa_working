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
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       trap3d
      program trap3d
c
      implicit integer (a-z)
      parameter ( number=500 )
      character*4096 ops
      character*8 prtflg
      character*1 itoc
      character*80 cpass, title, chrkey
      character*1600 card
      character*24 search
      character*128 fillam, filham
      character*2 timtyp
      logical posinp, logkey, prncof, prnply, prnwpt, prnh, prnh0 
      logical prnht, prall, chko, check, fix, calhmo, caltpr, calnls
      logical prgues, prvec, toau, useau, hamd, incore, icdiis
      logical prres, itsol, itdiag, dodiis, immed, prbufh, prtrn
      logical precon, prdvd, dvdall, prdiis
c     logical incosp
      real*8 z, fpkey, omegat, amass, natmin, natmax, left, right, scatl 
      real*8 pi, u0, tsize, endpts, tmax, der, todiis
      real*8 norm0, temp, thresh, cnverg, hbar, omega, cnvrg, table
      real*8 massau, lenau, timau, dum, scale
      common z(1)
      dimension ia(1), endpts(2,3), der(2,3), nmax(3), npt(3), fix(3)
      dimension nfix(3), temp(3), pleft(3), pright(3), dim(3), norm0(3)
      dimension prncof(3), prnply(3), prnwpt(3), prnh(3), check(3)
      dimension q(3), wt(3), a(3), b(3), p(3), dp(3), ddp(3)
      dimension pn(3), dpn(3), ddpn(3), eigc(3), work(4), eig(3)
      dimension nply(3), tr(3), omegat(3), left(3), right(3)
      dimension table(number,3), ham(3), v(3), u(3), ncon(3)
      dimension prdvd(10)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c
c          mass in Kilograms   length in Kilometers
      data massau, lenau, timau / 9.109558d-31, 5.291771d-11, 
     1                            2.418884d-17 /
      call drum
      write(iout,*)
      write(iout,*) '    three dimensional mean-field code          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c                       set print options
c
      prall=logkey(ops,'print=m7000=all',.false.,' ')
      calhmo=logkey(ops,'diagonalize-h0',.false.,' ')
      calnls=logkey(ops,'diagonalize-gross-pitaevski',.false.,' ')
      if(calhmo) then 
         write(iout,1) calhmo      
      else
         write(iout,2) calnls      
      endif
      numdim=intkey(ops,'number-of-dimensions',1,' ')
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
         prncof(1)=logkey(ops,'print=m7000=x-polynomial-coefficients',
     1                    .false.,' ')
         prncof(2)=logkey(ops,'print=m7000=y-polynomial-coefficients',
     1                    .false.,' ')
         prncof(3)=logkey(ops,'print=m7000=z-polynomial-coefficients',
     1                    .false.,' ')
         prnply(1)=logkey(ops,'print=m7000=x-polynomials',.false.,' ')
         prnply(2)=logkey(ops,'print=m7000y-polynomials',.false.,' ')
         prnply(3)=logkey(ops,'print=m7000=z-polynomials',.false.,' ')
         prnwpt(1)=logkey(ops,'print=m7000=x-points/weights',
     1                    .false.,' ')
         prnwpt(2)=logkey(ops,'print=m7000=y-points/weights',
     1                    .false.,' ')
         prnwpt(3)=logkey(ops,'print=m7000=z-points/weights',
     1                    .false.,' ')
         prnh0=logkey(ops,'print=m7000=h0',.false.,' ')
         prnht=logkey(ops,'print=m7000=h',.false.,' ')
         prnh(1)=logkey(ops,'print=m7000=x-hamiltonian',.false.,' ')
         prnh(2)=logkey(ops,'print=m7000=y-hamiltonian',.false.,' ')
         prnh(3)=logkey(ops,'print=m7000=z-hamiltonian',.false.,' ')
      endif
      chko=logkey(ops,'m7000=check-orthogonality',.false.,' ')
      if(chko) then
         check(1)=.true.
         check(2)=.true.
         check(3)=.true.
      else                  
         check(1)=logkey(ops,'check-x-orthogonality',.false.,' ')
         check(2)=logkey(ops,'check-y-orthogonality',.false.,' ')
         check(3)=logkey(ops,'check-z-orthogonality',.false.,' ')
      endif   
      toau=logkey(ops,'to-atomic-units',.false.,' ')
      useau=logkey(ops,'use-atomic-units',.false.,' ')
c     incosp=logkey(ops,'in-coordinate-space',.false.,' ')
      if(numdim.eq.1) then
         write(iout,3) prncof(1),
     1                 prnply(1),
     2                 prnwpt(1), prnh0,
     3                 prnh(1), 
     4                 check(1)
      endif
      if(numdim.eq.2) then
         write(iout,4) prncof(1), prncof(2),
     1                 prnply(1), prnply(2),
     2                 prnwpt(1), prnwpt(2), prnh0,
     3                 prnh(1), prnh(2), 
     4                 check(1), check(2)
      endif     
      if(numdim.eq.3) then
         write(iout,5) prncof(1), prncof(2), prncof(3),
     1                 prnply(1), prnply(2), prnply(3),
     2                 prnwpt(1), prnwpt(2), prnwpt(3), prnh0,
     3                 prnh(1), prnh(2), prnh(3), 
     4                 check(1), check(2), check(3)
      endif                     
c
c                set diis options
c      
      dodiis=.false.
      maxit=intkey(ops,'maximum-number-of-iterations',10,' ')
      cnvrg=fpkey(ops,'convergence',1.d-08,' ')
      if( posinp('$diis',cpass) ) then
         dodiis=.true.
         call cardin(card)
         maxit=intkey(card,'maximum-number-of-diis-iterations',10,' ')
         trunc=intkey(card,'diis-restart-value',5,' ') 
         immed=logkey(card,'diis-switch=off',.false.,' ')
         todiis=1.d0
         if(immed) then
            todiis=1.d+10
         endif
         cnvrg=fpkey(card,'diis-convergence',1.d-08,' ')
         prdiis=logkey(card,'print=diis',.false.,' ')
         onoff=logkey(card,'new-diis',.false.,' ')
         write(iout,6) maxit, trunc, todiis, cnvrg
      endif
c
c              set various program options
c      

      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
c
c         set trap configuration, numbers of atoms, various constants
c         and other physical parameters
c       
      if ( posinp('$trap',cpass) ) then
           call cardin(card)
c       frequency of trap in Hz           
           omegat(1)=fpkey(card,'x-trap-frequency',10.d0,' ')
           omegat(1)=omegat(1)*2.d0*pi
           omegat(2)=fpkey(card,'y-trap-frequency',10.d0,' ')
           omegat(2)=omegat(2)*2.d0*pi
           omegat(3)=fpkey(card,'z-trap-frequency',10.d0,' ')
           omegat(3)=omegat(3)*2.d0*pi
c       mass in kilograms of atom.  default is Cs.           
           amass=fpkey(card,'atomic-mass',2.2d-25,' ')
           natmin=fpkey(card,'minimum-number-of-atoms',0.d+00,' ')
           natmax=fpkey(card,'maximum-number-of-atoms',1.d+05,' ')
           nstep=intkey(card,'number-of-steps',10,' ')
c       scattering length in meters                              
           u0=fpkey(card,'potential-strength',2.0d-51,' ')
           scatl=amass*u0/(4.d0*hbar*hbar*pi)
      endif
c
c               spatial basis set information
c       
      if ( posinp('$x-polynomials',cpass) ) then
           call cardin(card)
           nmax(1)=intkey(card,'order-of-polynomials',10,' ')
           ncon(1)=nmax(1)
           npt(1)=intkey(card,'number-of-points',nmax(1),' ')
           fix(1)=logkey(card,'fix-end-points',.false.,' ')
           left(1)=fpkey(card,'left',-2.5d0,' ')
           right(1)=fpkey(card,'right',2.5d0,' ')
           der(1,1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2,1)=fpkey(card,'right-derivative',0.d0,' ')
           if (fix(1)) then
               nfix(1)=intkey(card,'number-of-fixed-points',
     1                      2,' ')
               endpts(1,1)=left(1)
               endpts(2,1)=right(1)
               call fparr(card,'end-points',endpts(1,1),
     1                    nfix(1),' ')
           endif
           pleft(1)=intkey(card,'order-of-leading-left-'//
     1                          'polynomials',1,' ')
           pright(1)=intkey(card,'order-of-leading-right-'//
     1                           'polynomials',0,' ')
      endif
      if(numdim.gt.1) then
         if ( posinp('$y-polynomials',cpass) ) then
              call cardin(card)
              nmax(2)=intkey(card,'order-of-polynomials',10,' ')
              ncon(2)=nmax(2)
              npt(2)=intkey(card,'number-of-points',nmax(2),' ')
              fix(2)=logkey(card,'fix-end-points',.false.,' ')
              left(2)=fpkey(card,'left',-2.5d0,' ')
              right(2)=fpkey(card,'right',2.5d0,' ')
              der(1,2)=fpkey(card,'left-derivative',0.d0,' ')
              der(2,2)=fpkey(card,'right-derivative',0.d0,' ')
              if (fix(2)) then
                  nfix(2)=intkey(card,'number-of-fixed-points',
     1                         2,' ')
                  endpts(1,2)=left(2)
                  endpts(2,2)=right(2)
                  call fparr(card,'end-points',endpts(1,2),
     1                       nfix(2),' ')
              endif
              pleft(2)=intkey(card,'order-of-leading-left-'//
     1                             'polynomials',1,' ')
              pright(2)=intkey(card,'order-of-leading-right-'//
     1                              'polynomials',0,' ')
         endif
      endif
      if(numdim.gt.2) then
         if ( posinp('$z-polynomials',cpass) ) then
              call cardin(card)
              nmax(3)=intkey(card,'order-of-polynomials',10,' ')
              ncon(3)=nmax(3)
              npt(3)=intkey(card,'number-of-points',nmax(3),' ')
              fix(3)=logkey(card,'fix-end-points',.false.,' ')
              left(3)=fpkey(card,'left',-2.5d0,' ')
              right(3)=fpkey(card,'right',2.5d0,' ')
              der(1,3)=fpkey(card,'left-derivative',0.d0,' ')
              der(2,3)=fpkey(card,'right-derivative',0.d0,' ')
              if (fix(3)) then
                  nfix(3)=intkey(card,'number-of-fixed-points',
     1                         2,' ')
                  endpts(1,3)=left(3)
                  endpts(2,3)=right(3)
                  call fparr(card,'end-points',endpts(1,3),
     1                       nfix(3),' ')
              endif
              pleft(3)=intkey(card,'order-of-leading-left-'//
     1                             'polynomials',1,' ')
              pright(3)=intkey(card,'order-of-leading-right-'//
     1                              'polynomials',0,' ')
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
            endpts(1,i)=left(i)
            endpts(2,i)=right(i)            
            omegat(i)=omegat(i)*timau
 10      continue            
         amass=amass/massau        
      endif
      if(useau) then
c
c        if this option is used then tsize is entered in atomic units
c
         write(iout,*) 'assuming atomic units'
         hbar=1.d0
         amass=1.d0
         do 20 i=1,numdim
            omegat(i)=1.d0
            endpts(1,i)=left(i)
            endpts(2,i)=right(i)            
 20      continue            
      endif
      if(numdim.ge.1) then
         write(iout,7) amass, omegat(1), natmin, natmax, 
     1                 scatl
         title='polynomial information for first dimension'
         if (fix(1)) then
             if(nfix(1).eq.1) then
                write(iout,11) title, nmax(1), npt(1),
     1                         left(1), right(1), 
     2                         der(1,1), der(2,1) , pleft(1), 
     3                         pright(1), nfix(1), endpts(1,1)
             endif
             if(nfix(1).eq.2) then
               write(iout,12) title, nmax(1), npt(1), 
     1                        left(1), right(1), 
     2                        der(1,1), der(2,1) , pleft(1), 
     3                        pright(1), nfix(1), endpts(1,1),
     4                        endpts(2,1)
             endif 
         endif
      endif
      if(numdim.ge.2) then
         write(iout,8) amass, (omegat(i),i=1,2), natmin, natmax, 
     1                 scatl
         title='polynomial information for second dimension'
         if (fix(2)) then
             if(nfix(2).eq.1) then
                write(iout,11) title, nmax(2), npt(2),
     1                         left(2), right(2), 
     2                         der(1,2), der(2,2) , pleft(2), 
     3                        pright(2), nfix(2), endpts(1,2)
             endif
             if(nfix(2).eq.2) then
                write(iout,12) title, nmax(2), npt(2),
     1                         left(2), right(2), 
     2                         der(1,2), der(2,2) , pleft(2), 
     3                         pright(2), nfix(2), endpts(1,2),
     4                         endpts(2,2)
             endif 
         endif
      endif
      if(numdim.eq.3) then
         write(iout,9) amass, (omegat(i),i=1,3), natmin, natmax, 
     1                 scatl
         title='polynomial information for third dimension'
         if (fix(3)) then
             if(nfix(3).eq.1) then
                write(iout,11) title, nmax(3), npt(3),
     1                         left(3), right(3), 
     2                         der(1,3), der(2,3) , pleft(3), 
     3                        pright(3), nfix(3), endpts(1,3)
             endif
             if(nfix(3).eq.2) then
                write(iout,12) title, nmax(3), npt(3),
     1                         left(3), right(3), 
     2                         der(1,3), der(2,3) , pleft(3), 
     3                         pright(3), nfix(3), endpts(1,3),
     4                         endpts(2,3)
             endif 
         endif
      endif  
      if( posinp('$3dim',cpass) ) then
          call cardin(card)
      endif
      ncon(1)=intkey(card,'number-of-x-contracted-functions',
     1               nmax(1),' ')
      ncon(2)=intkey(card,'number-of-y-contracted-functions',
     1               nmax(2),' ')
      ncon(3)=intkey(card,'number-of-z-contracted-functions',
     1               nmax(3),' ')
      dim(1)=max(nmax(1),npt(1))
      dim(2)=max(nmax(2),npt(2))
      dim(3)=max(nmax(3),npt(3))
      maxd=max(dim(1),dim(2),dim(3))
      nply(1)=nmax(1)
      nply(2)=nmax(2)
      nply(3)=nmax(3)
      n3d=1
      do 30 i=1,numdim
         n3d=n3d*ncon(i)
 30   continue 
      nroots=n3d         
c
c             set diagonalization procedures
c      
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      if(itdiag) then
         if(posinp('$davidson',cpass)) then
            call cardin(card)
         endif 
         nroots=intkey(card,'number-of-roots',1,' ')
         ntrials=intkey(card,'number-of-trial-vectors',nroots,' ')
         cnverg=fpkey(card,'convergence',1.d-08,' ')
         thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
         niter=intkey(card,'maximum-number-of-iterations',
     1                n3d,' ')
         lenbuf=intkey(card,'hamiltonian-buffer',
     1                 min(1000000,n3d*n3d),' ')
         prbufh=logkey(card,'print=buffered-hamiltonian',.false.,' ') 
         prtrl=logkey(card,'print=trials',.false.,' ') 
         hamd=logkey(card,'to-h0',.false.,' ')
         precon=logkey(card,'precondition',.false.,' ')
         n0=intkey(card,'to-h0=nkeep',n3d,' ')
         prdvd(1)=logkey(card,'print=davidson=trials',.false.,' ')
         prdvd(2)=logkey(card,'print=davidson=hamiltonian',.false.,' ')
         prdvd(3)=logkey(card,'print=davidson=overlaps',.false.,' ')
         prdvd(4)=logkey(card,'print=davidson=vectors',.false.,' ')
         prdvd(5)=logkey(card,'print=davidson=h-on-vectors',.false.,' ')
         prdvd(6)=logkey(card,'print=davidson=small-matrix',.false.,' ')
         prdvd(7)=logkey(card,'print=davidson=transformed-vectors',
     1                  .false.,' ')
         prdvd(8)=logkey(card,'print=davidson=transformed-h-on-vectors',
     1                   .false.,' ')
         prdvd(9)=logkey(card,'print=davidson=residuals',.false.,' ')
         prdvd(10)=logkey(card,'print=davidson=new-trial-vectors',
     1                   .false.,' ')
         dvdall=logkey(card,'print=davidson=all',.false.,' ')
         if(dvdall) then
            do 40 i=1,10
               prdvd(i)=.true.
 40         continue
         endif   
         call iosys ('read character "hamiltonian filename" from rwf',
     1                -1,0,0,filham)
         call iosys('open ham as new',0,0,0,filham)
         call iosys('write integer "buffer size" to ham',1,
     1               lenbuf,0,' ')           
         write(iout,13) nroots, thresh, cnverg, niter
      else
         nroots=intkey(ops,'number-of-roots',n3d,' ')
         write(iout,14) n3d, nroots
      endif
      ioff=1
c
c     note the alignment of variables.  it is important in passing these 
c     arrays later.
c 
      do 50 i=1,2
         if(numdim.eq.1) then
            q(1)=ioff
            wt(1)=q(1)+dim(1)
            a(1)=wt(1)+dim(1)
            b(1)=a(1)+nmax(1)+1
            p(1)=b(1)+nmax(1)+1
            dp(1)=p(1)+nmax(1)*dim(1)
            ddp(1)=dp(1)+nmax(1)*dim(1)
            pn(1)=ddp(1)+nmax(1)*dim(1)
            dpn(1)=pn(1)+nmax(1)*dim(1)
            ddpn(1)=dpn(1)+nmax(1)*dim(1)
            eigc(1)=ddpn(1)+nmax(1)*dim(1)
            tr(1)=eigc(1)+nmax(1)
            ham(1)=tr(1)+nmax(1)*nmax(1)
            eig(1)=ham(1)+nmax(1)*nmax(1) 
            v(1)=eig(1)+nmax(1)
            u(1)=v(1)+nmax(1)*nmax(1)
            ind=wpadti(u(1)+nmax(1)*nmax(1))
            idum=ind+2*n3d
            words=iadtwp(idum+2*n3d)            
         elseif(numdim.eq.2) then 
            q(1)=ioff
            q(2)=q(1)+dim(1)
            wt(1)=q(2)+dim(2)
            wt(2)=wt(1)+dim(1)          
            a(1)=wt(2)+dim(2)
            a(2)=a(1)+nmax(1)+1
            b(1)=a(2)+nmax(2)+1
            b(2)=b(1)+nmax(1)+1
            p(1)=b(2)+nmax(2)+1
            p(2)=p(1)+nmax(1)*dim(1)            
            dp(1)=p(2)+nmax(2)*dim(2)
            dp(2)=dp(1)+nmax(1)*dim(1)
            ddp(1)=dp(2)+nmax(2)*dim(2)
            ddp(2)=ddp(1)+nmax(1)*dim(1)
            pn(1)=ddp(2)+nmax(2)*dim(2)
            pn(2)=pn(1)+nmax(1)*dim(1)
            dpn(1)=pn(2)+nmax(2)*dim(2)
            dpn(2)=dpn(1)+nmax(1)*dim(1)
            ddpn(1)=dpn(2)+nmax(2)*dim(2)
            ddpn(2)=ddpn(1)+nmax(1)*dim(1)
            eigc(1)=ddpn(2)+nmax(2)*dim(2)
            eigc(2)=eigc(1)+nmax(1)
            tr(1)=eigc(2)+nmax(2)
            tr(2)=tr(1)+nmax(1)*nmax(1)
            ham(1)=tr(2)+nmax(2)*nmax(2)
            ham(2)=ham(1)+nmax(1)*nmax(1)
            eig(1)=ham(2)+nmax(2)*nmax(2) 
            eig(2)=eig(1)+nmax(1)
            v(1)=eig(2)+nmax(2) 
            v(2)=v(1)+nmax(1)*nmax(1) 
            u(1)=v(2)+nmax(2)*nmax(2)
            u(2)=u(1)+nmax(1)*nmax(1)
            ind=wpadti(u(2)+nmax(2)*nmax(2))
            idum=ind+3*n3d
            words=iadtwp(idum+3*n3d)
         elseif(numdim.eq.3) then
            q(1)=ioff
            q(2)=q(1)+dim(1)
            q(3)=q(2)+dim(2)
            wt(1)=q(3)+dim(3)
            wt(2)=wt(1)+dim(1)          
            wt(3)=wt(2)+dim(2)          
            a(1)=wt(3)+dim(3)
            a(2)=a(1)+nmax(1)+1
            a(3)=a(2)+nmax(2)+1
            b(1)=a(3)+nmax(3)+1
            b(2)=b(1)+nmax(1)+1
            b(3)=b(2)+nmax(2)+1
            p(1)=b(3)+nmax(3)+1
            p(2)=p(1)+nmax(1)*dim(1)            
            p(3)=p(2)+nmax(2)*dim(2)            
            dp(1)=p(3)+nmax(3)*dim(3)
            dp(2)=dp(1)+nmax(1)*dim(1)
            dp(3)=dp(2)+nmax(2)*dim(2)
            ddp(1)=dp(3)+nmax(3)*dim(3)
            ddp(2)=ddp(1)+nmax(1)*dim(1)
            ddp(3)=ddp(2)+nmax(2)*dim(2)
            pn(1)=ddp(3)+nmax(3)*dim(3)
            pn(2)=pn(1)+nmax(1)*dim(1)
            pn(3)=pn(2)+nmax(2)*dim(2)
            dpn(1)=pn(3)+nmax(3)*dim(3)
            dpn(2)=dpn(1)+nmax(1)*dim(1)
            dpn(3)=dpn(2)+nmax(2)*dim(2)
            ddpn(1)=dpn(3)+nmax(3)*dim(3)
            ddpn(2)=ddpn(1)+nmax(1)*dim(1)
            ddpn(3)=ddpn(2)+nmax(2)*dim(2)
            eigc(1)=ddpn(3)+nmax(3)*dim(3)
            eigc(2)=eigc(1)+nmax(1)
            eigc(3)=eigc(2)+nmax(2)
            tr(1)=eigc(3)+nmax(3)
            tr(2)=tr(1)+nmax(1)*nmax(1)
            tr(3)=tr(2)+nmax(2)*nmax(2)
            ham(1)=tr(3)+nmax(3)*nmax(3)
            ham(2)=ham(1)+nmax(1)*nmax(1)
            ham(3)=ham(2)+nmax(2)*nmax(2)
            eig(1)=ham(3)+nmax(3)*nmax(3) 
            eig(2)=eig(1)+nmax(1) 
            eig(3)=eig(2)+nmax(2) 
            v(1)=eig(3)+nmax(3)
            v(2)=v(1)+nmax(1)*nmax(1)
            v(3)=v(2)+nmax(2)*nmax(2)
            u(1)=v(3)+nmax(3)*nmax(3)
            u(2)=u(1)+nmax(1)*nmax(1)
            u(3)=u(2)+nmax(2)*nmax(2)
            ind=wpadti(u(3)+nmax(3)*nmax(3))
            idum=ind+4*n3d
            words=iadtwp(idum+4*n3d)
         endif    
         work(1)=words
         work(2)=work(1)+2*maxd*maxd
         work(3)=work(2)+2*maxd*maxd
         words=work(3)+max(3*n3d,2*maxd*maxd)
         if(.not.itdiag) then
            hamtot=words
            eigtot=hamtot+n3d*n3d
            vtot=eigtot+n3d
            words=vtot+n3d
            if(calnls) then
               rho=words
               psi=rho+n3d*n3d
               fmat=psi+n3d*(maxit+1)
               err=fmat+n3d*n3d*(maxit+1)
               bb=err+n3d*n3d*(maxit+1)
               bbtmp=bb+(maxit+1)*(maxit+1)
               sol=bbtmp+(maxit+1)*(maxit+1)
               vnl=sol+maxit+1
               ipvt=wpadti(vnl+n3d*(maxit+1))
               words=iadtwp(ipvt+n3d)
            endif
         else
            diag=words
            vtot=diag+n3d
            etrial=vtot+n3d
            trials=etrial+n3d
            pvec=trials+n3d*ntrials
            hpvec=pvec+n3d*niter
            vec=hpvec+niter*n3d
            bmat=vec+n3d*niter
            bmatm=bmat+niter*niter
            resid=bmatm+niter*niter
            eigtot=resid+n3d*niter
            hbuf=eigtot+n3d
            ibuf=wpadti(hbuf+lenbuf)
            words=iadtwp(ibuf+2*lenbuf)
            if(calnls) then
               vnl=words
               psi=vnl+n3d*niter
               words=psi+n3d
            endif
            if(precon) then
               bigham=words
               words=wpadti(bigham+n3d*n3d)
            endif               
         endif
         words=wpadti(words)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'poly',0)
         endif
   50 continue
      icdiis=.true.
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      do 60 i=1,numdim
         temp(1)=-1.d0
         temp(2)=1.d0
         if(nfix(i).eq.1) then
            temp(1)=1.d0
         endif
         call gaussq('legendre',npt(i),0.d0,0.d0,nfix(i),temp,
     1                z(work(1)),z(q(i)),z(wt(i)))
c     
c           convert to points and weights on (left,right) if needed
c
         call convt(z(q(i)),z(wt(i)),endpts(1,i),endpts(2,i),
     1              norm0(i),npt(i),prnwpt(i))
c
c           calculate the recursion coefficients
c
         call cpoly(z(p(i)),z(q(i)),z(wt(i)),z(a(i)),z(b(i)),
     1              endpts(1,i),endpts(2,i),z(work(1)),nmax(i),npt(i),
     2              pleft(i),pright(i))
         if (prncof(i)) then
             title='a coefficients for '//itoc(i)//' dimension'
             call prntfm(title,z(a(i)+1),nmax(i),1,nmax(i),1,iout)
             title='b coefficients for '//itoc(i)//' dimension'
             call prntfm(title,z(b(i)+1),nmax(i)-1,1,nmax(i)-1,1,iout)
         endif
c
c           calculate the polynomials and their first and 
c                        second derivatives.
c
         call gpoly(z(p(i)),z(dp(i)),z(ddp(i)),z(q(i)),z(a(i)),
     1              z(b(i)),endpts(1,i),endpts(2,i),
     2              pleft(i),pright(i),nmax(i),npt(i),.false.)
         call diagx(z(p(i)),z(dp(i)),z(ddp(i)),z(a(i)+1),
     1              z(b(i)+1),z(pn(i)),z(dpn(i)),z(ddpn(i)),
     2              z(eigc(i)),z(work(1)),z(tr(i)),nmax(i),
     3              npt(i),prnh(i))
         if (prnply(i)) then
             title='polynomials '//itoc(i)//' dimension'
             call prntfm(title,z(p(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='first derivative of polynomials '
     1              //itoc(i)//' dimension'
             call prntfm(title,z(dp(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='second derivative of polynomials '
     1              //itoc(i)//' dimension'
             call prntfm(title,z(ddp(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='dvr polynomials '//itoc(i)//' dimension'
             call prntfm(title,z(pn(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='first derivative of dvr polynomials '
     1              //itoc(i)//' dimension'
             call prntfm(title,z(dpn(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='second derivative of dvr polynomials '
     1              //itoc(i)//' dimension'
             call prntfm(title,z(ddpn(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
         endif
         if(check(i)) then
            call chk(z(p(i)),z(wt(i)),z(work(1)),z(work(2)),
     1               nmax(i),npt(i))
            title='overlap matrix polynomials '//itoc(i)//' dimension'
            call prntfm(title,z(work(2)),nmax(i),nmax(i),nmax(i),
     1                                nmax(i),iout)
            call chk(z(pn(i)),z(wt(i)),z(work(1)),z(work(2)),
     1               nmax(i),npt(i))
            title='overlap matrix dvr polynomials '//itoc(i)
     1             //' dimension'
            call prntfm(title,z(work(2)),nmax(i),nmax(i),nmax(i),
     1                                nmax(i),iout)
         endif
         search='look for v0'
c        if(incosp) then
c           call vmcosp(z(pn(i)),z(q(i)),z(wt(i)),z(v(i)),z(work(1)),
c    1                  z(work(2)),hbar,amass,scale,nmax(i),npt(i),
c    2                  prnh0,ops)
c           title='potential matrix in FBR representation'
c           call prntrm(title,z(v(i)),nmax(i),nmax(i),nmax(i),
c    1                                nmax(i),iout)     
c           call ham0(z(pn(i)),z(dpn(i)),z(ddpn(i)),z(wt(i)),z(ham(i)),
c    1                z(v(i)),z(u(i)),z(eig(i)),z(work(1)),hbar,amass,
c    2                nmax(i),npt(i),prnh0)
c                                         ,incosp)
c        else
         call vmat(z(eigc(i)),dum,dum,z(v(i)),z(work(1)),hbar,amass,
     1             scale,nmax(1),nmax,numdim,prnh0,search,ops)
c        title='potential matrix in DVR representation'
c           call prntrm(title,z(v(i)),nmax(i),1,nmax(i),1,iout)     
         call ham0(z(pn(i)),z(dpn(i)),z(ddpn(i)),z(wt(i)),z(ham(i)),
     1             z(v(i)),z(u(i)),z(eig(i)),z(work(1)),hbar,amass,
     2             nmax(i),npt(i),prnh0)
c                                      ,incosp)
c        endif
         call vscale(z(eig(i)),z(eig(i)),scale,nmax(1)) 
         title='eigenvalues of unperturbed hamiltonian '//itoc(i)
     1          //' dimension'
         call prntfm(title,z(eig(i)),nmax(i),1,nmax(i),1,iout)
         call vscale(z(eig(i)),z(eig(i)),1.d0/scale,nmax(1)) 
         if (prnh0) then
             title='eigenvectors of unperturbed hamiltonian '//itoc(i)
     1          //' dimension'
             call prntfm(title,z(u(i)),nmax(i),nmax(i),nmax(i),
     1                   nmax(i),iout)
         endif         
c        evaluate dvr polynomials at dvr eigenvalues     
         call newply(z(p(i)),z(dp(i)),z(ddp(i)),z(a(i)),z(b(i)),
     1               z(pn(i)),z(dpn(i)),z(ddpn(i)),z(eigc(i)),
     2               z(tr(i)),endpts(1,i),endpts(2,i),
     3               pleft(i),pright(i),nmax(i))
         if (prnply(i)) then
             title='dvr polynomials '//itoc(i)//' dimension'         
             call prntfm(title,z(pn(i)),nmax(i),nmax(i),nmax(i),
     1                   nmax(i),iout)
             title='first derivative of dvr polynomials '//
     1              itoc(i)//' dimension'
             call prntfm(title,z(dpn(i)),nmax(i),nmax(i),nmax(i),
     1                   nmax(i),iout)
             title='second derivative of dvr polynomials '//
     1              itoc(i)//' dimension'
             call prntfm(title,z(ddpn(i)),nmax(i),nmax(i),nmax(i),
     1                   nmax(i),iout)
         endif            
 60   continue
      call setind(ia(ind),nmax,n3d,numdim)
      search='look for v1'
      if(numdim.eq.1) then
         call vmat(z(eigc(1)),dum,dum,z(vtot),z(work(3)),
     1             hbar,amass,scale,n3d,nmax,numdim,prnh0,
     2             search,ops)   
      elseif(numdim.eq.2) then
         call vmat(z(eigc(1)),z(eigc(2)),dum,z(vtot),z(work(3)),
     1             hbar,amass,scale,n3d,nmax,numdim,prnh0,
     2             search,ops)   
      elseif(numdim.eq.3) then
         call vmat(z(eigc(1)),z(eigc(2)),z(eigc(3)),z(vtot),
     1             z(work(3)),hbar,amass,scale,n3d,nmax,
     2             numdim,prnh0,search,ops)  
      endif
      if(calnls) then
         if(.not.itdiag) then
            if(numdim.eq.1) then
               call nlschr(z(fmat),z(eigtot),z(rho),z(psi),z(err),
     1                     z(vnl),z(pn(1)),dum,dum,
     2                     z(eigc(1)),dum,dum,
     3                     z(ham(1)),dum,dum,z(vtot),
     4                     z(bb),z(bbtmp),z(sol),z(work(1)),z(work(2)),
     5                     z(work(3)),table,ia(ipvt),natmin,natmax,
     6                     hbar,amass,omegat,scatl,todiis,cnvrg,n3d,
     7                     maxit,nstep,number,numdim,nmax,trunc,
     8                     prdiis,dodiis,icdiis,onoff)
            elseif(numdim.eq.2) then
               call nlschr(z(fmat),z(eigtot),z(rho),z(psi),z(err),
     1                     z(vnl),z(pn(1)),z(pn(2)),dum,
     2                     z(eigc(1)),z(eigc(2)),dum,
     3                     z(ham(1)),z(ham(2)),dum,z(vtot),
     4                     z(bb),z(bbtmp),z(sol),z(work(1)),z(work(2)),
     5                     z(work(3)),table,ia(ipvt),natmin,natmax,
     6                     hbar,amass,omegat,scatl,todiis,cnvrg,n3d,
     7                     maxit,nstep,number,numdim,nmax,trunc,
     8                     prdiis,dodiis,icdiis,onoff)
            elseif(numdim.eq.3) then
               call nlschr(z(fmat),z(eigtot),z(rho),z(psi),z(err),
     1                     z(vnl),z(pn(1)),z(pn(2)),z(pn(3)),
     2                     z(eigc(1)),z(eigc(2)),z(eigc(3)),
     3                     z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),
     4                     z(bb),z(bbtmp),z(sol),z(work(1)),z(work(2)),
     5                     z(work(3)),table,ia(ipvt),natmin,natmax,
     6                     hbar,amass,omegat,scatl,todiis,cnvrg,n3d,
     7                     maxit,nstep,number,numdim,nmax,trunc,
     8                     prdiis,dodiis,icdiis,onoff)
            endif
         else
            if(numdim.eq.1) then
               call vtrial(z(u(1)),dum,dum,z(eig(1)),dum,dum,z(trials),
     1                     z(etrial),ia(idum),nmax,numdim,n3d,
     2                     ntrials,prtrl)
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),dum,dum,z(vtot),z(hbuf),ia(ibuf),
     1                    ia(ind),z(diag),lenbuf,nmax,numdim,nmax(1),
     2                    prbufh,nel,incore,title)
               if(precon) then
                  call trmat(z(hbuf),ia(ibuf),z(diag),z(vtot),
     1                       z(bigham),z(eigtot),z(bigham),z(resid),
     2                       z(trials),z(etrial),scale,numdim,n3d,
     3                       ntrials,lenbuf,incore)
               endif
               call nldvd(z(hbuf),ia(ibuf),z(diag),z(vnl),z(trials),
     1                    z(etrial),z(pvec),z(hpvec),z(vec),z(resid),
     2                    z(psi),z(pn(1)),dum,dum,
     3                    z(eigc(1)),dum,dum,
     4                    z(ham(1)),dum,dum,
     5                    z(eig(1)),dum,dum,
     6                    z(u(1)),dum,dum,
     7                    z(bmat),z(bmatm),z(eigtot),
     8                    z(work(1)),z(work(2)),z(work(3)),table,natmin,
     9                    natmax,hbar,amass,omegat,scatl,cnverg,
     x                    thresh,nstep,number,nmax,numdim,n3d,nroots,
     x                    ntrials,niter,lenbuf,hamd,incore,prdvd) 
            else
               call lnkerr('iterative NLSE not yet implimented')
            endif
         endif
      else
         if(.not.itdiag) then    
            if(numdim.eq.1) then
               call lschr1(z(ham(1)),z(eigtot),z(vtot),nmax(1))
            elseif(numdim.eq.2) then
               call lschr2(z(ham(1)),z(ham(2)),z(vtot),z(hamtot),
     1                     z(eigtot),ia(ind),n3d,nmax,numdim,
     2                     prnht,ops)
            elseif(numdim.eq.3) then
               call lschr3(z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),
     1                     z(hamtot),z(eigtot),ia(ind),n3d,nmax,
     2                     numdim,prnht,ops)
            endif
         else
            if(numdim.eq.1) then
               call vtrial(z(u(1)),dum,dum,z(eig(1)),dum,dum,z(trials),
     1                     z(etrial),ia(idum),nmax,numdim,n3d,
     2                     ntrials,prtrl)
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),dum,dum,z(vtot),z(hbuf),ia(ibuf),
     1                    ia(ind),z(diag),lenbuf,nmax,numdim,nmax(1),
     2                    prbufh,nel,incore,title)
               if(precon) then
                  call trmat(z(hbuf),ia(ibuf),z(diag),z(vtot),
     1                       z(bigham),z(eigtot),z(bigham),z(resid),
     2                       z(trials),z(etrial),scale,numdim,n3d,
     3                       ntrials,lenbuf,incore)
               endif
               call dvd(z(hbuf),ia(ibuf),z(diag),z(trials),z(eig(1)),
     1                 z(pvec),z(hpvec),z(vec),z(bmat),z(bmatm),
     2                 z(eigtot),z(eig(1)),dum,dum,z(u(1)),dum,dum,
     3                 z(resid),z(work(1)),z(work(2)),z(work(3)),
     4                 scale,cnverg,thresh,nmax,numdim,nmax(1),nroots,
     5                 ntrials,niter,lenbuf,hamd,incore,prdvd) 
            elseif(numdim.eq.2) then
               call vtrial(z(u(1)),z(u(2)),dum,z(eig(1)),z(eig(2)),
     1                     dum,z(trials),z(etrial),ia(idum),nmax,
     2                     numdim,n3d,ntrials,prtrl)
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),z(ham(2)),dum,z(vtot),z(hbuf),
     1                    ia(ibuf),ia(ind),z(diag),lenbuf,nmax,
     2                    numdim,n3d,prbufh,nel,incore,title)
               if(precon) then
                  call trmat(z(hbuf),ia(ibuf),z(diag),z(vtot),
     1                       z(bigham),z(eigtot),z(bigham),z(resid),
     2                       z(trials),z(etrial),scale,numdim,n3d,
     3                       ntrials,lenbuf,incore)
               endif
               call dvd(z(hbuf),ia(ibuf),z(diag),z(trials),z(etrial),
     1                 z(pvec),z(hpvec),z(vec),z(bmat),z(bmatm),
     2                 z(eigtot),z(eig(1)),z(eig(2)),dum,z(u(1)),
     3                 z(u(2)),dum,z(resid),z(work(1)),z(work(2)),
     4                 z(work(3)),scale,cnverg,thresh,nmax,numdim,
     5                 n3d,nroots,ntrials,niter,lenbuf,hamd,incore,
     6                 prdvd)
            elseif(numdim.eq.3) then
               call vtrial(z(u(1)),z(u(2)),z(u(3)),z(eig(1)),z(eig(2)),
     1                     z(eig(3)),z(trials),z(etrial),ia(idum),nmax,
     2                     numdim,n3d,ntrials,prtrl)
               title='buffered hamiltonian for '//itoc(numdim)
     1                //'d dimensional hamiltonian'
               call hamil(z(ham(1)),z(ham(2)),z(ham(3)),z(vtot),z(hbuf),
     1                    ia(ibuf),ia(ind),z(diag),lenbuf,nmax,
     2                    numdim,n3d,prbufh,nel,incore,title)
               if(precon) then
                  call trmat(z(hbuf),ia(ibuf),z(diag),z(vtot),
     1                       z(bigham),z(eigtot),z(bigham),z(resid),
     2                       z(trials),z(etrial),scale,numdim,n3d,
     3                       ntrials,lenbuf,incore)
                endif
               call dvd(z(hbuf),ia(ibuf),z(diag),z(trials),z(etrial),
     1                 z(pvec),z(hpvec),z(vec),z(bmat),z(bmatm),
     2                 z(eigtot),z(eig(1)),z(eig(2)),z(eig(3)),z(u(1)),
     3                 z(u(2)),z(u(3)),z(resid),z(work(1)),z(work(2)),
     4                 z(work(3)),scale,cnverg,thresh,nmax,numdim,
     5                 n3d,nroots,ntrials,niter,lenbuf,hamd,incore,
     6                 prdvd)
            endif
         endif
      endif
      call vscale(z(eigtot),z(eigtot),scale,nroots)
      title='eigenvalues of hamiltonian'
      call prntfm(title,z(eigtot),nroots,1,nroots,1,iout)
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,15x,'program options',/,/,5x,
     1             'diagonalize zeroth-order hamiltonian = ',l1)  
 2    format(/,15x,'program options',/,/,5x,
     1             'diagonalize gross-pitaevski hamiltonian = ',l1)
 3    format(/,15x,'code options',/,
     1       /,5x,'print x coefficients                  = ',l1,
     2       /,5x,'print x polynomials                   = ',l1,
     3       /,5x,'print x points/weights                = ',l1,
     4       /,5x,'print spatial hamiltonian information = ',l1, 
     5       /,5x,'print x eigenvectors                  = ',l1, 
     6       /,5x,'check x orthonormality                = ',l1)
 4    format(/,15x,'code options',/,
     1       /,5x,'print x coefficients                  = ',l1,
     2       /,5x,'print y coefficients                  = ',l1,
     3       /,5x,'print x polynomials                   = ',l1,
     4       /,5x,'print y polynomials                   = ',l1,
     5       /,5x,'print x points/weights                = ',l1,
     6       /,5x,'print y points/weights                = ',l1,
     7       /,5x,'print spatial hamiltonian information = ',l1, 
     8       /,5x,'print x eigenvectors                  = ',l1, 
     9       /,5x,'print y eigenvectors                  = ',l1, 
     x       /,5x,'check x orthonormality                = ',l1, 
     x       /,5x,'check y orthonormality                = ',l1) 
 5    format(/,15x,'code options',/,
     1       /,5x,'print x coefficients                  = ',l1,
     2       /,5x,'print y coefficients                  = ',l1,
     3       /,5x,'print z coefficients                  = ',l1,
     4       /,5x,'print x polynomials                   = ',l1,
     5       /,5x,'print y polynomials                   = ',l1,
     6       /,5x,'print z polynomials                   = ',l1,
     7       /,5x,'print x points/weights                = ',l1,
     8       /,5x,'print y points/weights                = ',l1,
     9       /,5x,'print z points/weights                = ',l1,
     x       /,5x,'print spatial hamiltonian information = ',l1, 
     x       /,5x,'print x eigenvectors                  = ',l1, 
     x       /,5x,'print y eigenvectors                  = ',l1, 
     x       /,5x,'print z eigenvectors                  = ',l1, 
     x       /,5x,'check x orthonormality                = ',l1, 
     x       /,5x,'check y orthonormality                = ',l1, 
     x       /,5x,'check z orthonormality                = ',l1) 
 6    format(/,15x,'diis options',/,/,5x,
     1             'maximum number of iterations = ',i4,/,5x,
     2             'restart diis at iteration    = ',i2,/,5x,
     3             'diis begins when error       = ',e15.8,/,5x,
     4             'diis convergence criterion   = ',e15.8)
 7    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap x frequency                 = ',e15.8,/,5x,
     3             'minimum number of atoms in trap  = ',e15.8,/,5x,
     4             'maximum number of atoms in trap  = ',e15.8,/,5x,
     5             'scattering length                = ',e15.8)
 8    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap x frequency                 = ',e15.8,/,5x,
     3             'trap y frequency                 = ',e15.8,/,5x,
     4             'minimum number of atoms in trap  = ',e15.8,/,5x,
     5             'maximum number of atoms in trap  = ',e15.8,/,5x,
     6             'scattering length                = ',e15.8)
 9    format(/,15x,'trap configuration and atomic parameters',/,/,5x,
     1             'atomic mass                      = ',e15.8,/,5x,
     2             'trap x frequency                 = ',e15.8,/,5x,
     3             'trap y frequency                 = ',e15.8,/,5x,
     4             'trap z frequency                 = ',e15.8,/,5x,
     5             'minimum number of atoms in trap  = ',e15.8,/,5x,
     6             'maximum number of atoms in trap  = ',e15.8,/,5x,
     7             'scattering length                = ',e15.8)     
 11   format(/,15x,a60,/,/,5x,
     1             'order of polynomials           = ',i3,/,5x,
     2             'number of integration points   = ',i3,/,5x,
     3             'left boundary                  = ',e15.8,/,5x,
     4             'right boundary                 = ',e15.8,/,5x,
     5             'derivative at left boundary    = ',e15.8,/,5x,
     6             'derivative at right boundary   = ',e15.8,/,5x,
     7             'order-of-left-polynomials      = ',i2,/,5x,
     8             'order-of-right-polynomials     = ',i2,/,5x,
     9             'number fixed end points        = ',i2,/,5x,     
     x             'end point                      = ',e15.8)
 12   format(/,15x,a60,/,/,5x,
     1             'order of polynomials           = ',i3,/,5x,
     2             'number of integration points   = ',i3,/,5x,
     3             'left boundary                  = ',e15.8,/,5x,
     4             'right boundary                 = ',e15.8,/,5x,
     5             'derivative at left boundary    = ',e15.8,/,5x,
     6             'derivative at right boundary   = ',e15.8,/,5x,
     7             'order-of-left-polynomials      = ',i2,/,5x,
     8             'order-of-right-polynomials     = ',i2,/,5x,
     9             'number fixed end points        = ',i2,/,5x,     
     x             'left end point                 = ',e15.8,/,5x,
     x             'right end point                = ',e15.8)
 13   format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'overlap tolerance                  = ',e15.8,/,5x,
     3             'convergence criterion              = ',e15.8,/,5x,
     4             'maximum number of iterations       = ',i6)
 14   format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix  = ',i3,/,5x,
     2             'number of roots = ',i3)

      end
c               call nldvd(z(hbuf),ia(ibuf),z(diag),z(vnl),z(trials),
c     1                    z(etrial),z(pvec),z(hpvec),z(vec),z(resid),
c     2                    z(psi),z(pn(1)),z(pn(2)),z(pn(3)),
c     3                    z(eigc(1)),z(eigc(2)),z(eigc(3)),
c     4                    z(ham(1)),z(ham(2)),z(ham(3)),
c     5                    z(eig(1)),z(eig(2)),z(eig(3)),
c     6                    z(u(1)),z(u(2)),z(u(3)),z(v),
c     7                    z(bmat),z(bmatm),z(eigtot),
c     8                    z(work(1)),z(work(2)),z(work(3)),table,natmin,
c     9                    natmax,hbar,amass,omegat,scatl,scale,cnverg,
c     x                    thresh,nstep,number,nmax,numdim,n3d,nroots,
c     x                    ntrials,niter,lenbuf,hamd,incore,ops) 



