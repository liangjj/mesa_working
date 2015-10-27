*deck tprop.f 
c***begin prologue     tprop
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, propoagation, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            driver for time propoagation
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       tprop
      program tprop
c
      implicit integer (a-z)
      character*4096 ops
      character*8 prtflg
      character*80 cpass, title, chrkey
      character*1600 card
      character*16 type, trials
      character*24 pottyp
      character*128 fillam, filham
      character*2 itoc, timtyp
      logical posinp, logkey, prncof, prnply, prnwpt, prnh, prnhmo 
      logical check, fix, decomp, calhmo, caltpr
      logical prham, prgues, prvec
      logical prres, prnch, itsolv
      common z(1)
      dimension ia(1), endpts(2,2), der(2,2), nmax(2), npt(2), fix(2)
      dimension nfix(2), temp(2), pleft(2), pright(2), dim(2), norm0(2)
      dimension prncof(2), prnply(2), prnwpt(2), prnh(2), check(2)
      dimension decomp(2), q(2), wt(2), a(2), b(2), p(2), dp(2), ddp(2)
      dimension pn(2), dpn(2), ddpn(2), eigc(2), work(4), ham(2), eig(2)
      dimension nply(2)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      data pi/3.14159265358979323844d0/
      data hbar/1.05450d-34/
      real*8 z, fpkey, omegat, amass, atomx, left, right, scatl 
      real*8 pi, u0, tsize, endpts, tmax, der
      real*8 norm0, temp, thresh, cnverg, hbar, omega
c
      call drum
      write(iout,*)
      write(iout,*) '    time dependent mean-field code          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      prncof(1)=logkey(ops,'print=m6285=radial-polynomial-coefficients',
     1                 .false.,' ')
      prncof(2)=logkey(ops,'print=m6285=time-polynomial-coefficients',
     1                 .false.,' ')
      prnply(1)=logkey(ops,'print=m6285=radial-polynomials',.false.,' ')
      prnply(2)=logkey(ops,'print=m6285=time-polynomials',.false.,' ')
      prnwpt(1)=logkey(ops,'print=m6285=radial-points/weights',
     1                 .false.,' ')
      prnwpt(2)=logkey(ops,'print=m6285=time-points/weights',
     1                 .false.,' ')
      prnhmo=logkey(ops,'print=m6285=harmonic-oscillator',.false.,' ')
      prntim=logkey(ops,'print=m6285=time-dependent-results',
     1              .false.,' ')
      prham=logkey(ops,'print=m6285=hamiltonian',.false.,' ')
      prnh(1)=logkey(ops,'print=m6285=r-hamiltonian',.false.,' ')
      prnh(1)=logkey(ops,'print=m6285=t-hamiltonian',.false.,' ')
      prnch=logkey(ops,'print=m6285=complex-hamiltonian',.false.,' ')
      check(1)=logkey(ops,'check-radial-orthogonality',.false.,' ')
      check(2)=logkey(ops,'check-time-orthogonality',.false.,' ')
      itsolv=logkey(ops,'iterative-solution',.false.,' ')
      calhmo=logkey(ops,'diagonalize-harmonic-oscillator',.false.,' ')
      caltpr=logkey(ops,'solve-time-dependent-equation',.false.,' ')
      if(caltpr) then
         timtyp=chrkey(ops,'time-potential','0',' ')
      endif
      pottyp=chrkey(ops,'potential-type','none',' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$trap',cpass) ) then
           call cardin(card)
           omegat=fpkey(card,'trap-frequency',10.d0,' ')
           amass=fpkey(card,'atomic-mass',23.d0,' ')
           atomx=intkey(card,'log-base-10-of-number-of-atoms',5,' ')
           natoms=10.d0**atomx            
           scatl=fpkey(card,'scattering-length',-20.d0,' ')
           tsize=fpkey(card,'trap-size',10.d0,' ')
           u0=4.d0*pi*scatl/amass
           state=intkey(card,'initial-wavepacket',0,' ')
           write(iout,1) amass, omegat, natoms, scatl, tsize 
      endif 
      if ( posinp('$radial-polynomials',cpass) ) then
           call cardin(card)
           nmax(1)=intkey(card,'order-of-polynomials',10,' ')
           npt(1)=intkey(card,'number-of-points',nmax(1),' ')
           fix(1)=logkey(card,'fix-end-points',.false.,' ')
           decomp(1)=logkey(card,'expand',.false.,' ')
           left=0.d0
           right=tsize
           der(1,1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2,1)=fpkey(card,'right-derivative',0.d0,' ')
           if (fix(1)) then
               nfix(1)=intkey(card,'number-of-fixed-points',
     1                      2,' ')
               endpts(1,1)=left
               endpts(2,1)=right
               call fparr(card,'end-points',endpts(1,1),
     1                    nfix(1),' ')
           endif
           pleft(1)=intkey(card,'order-of-leading-left-'//
     1                          'polynomials',1,' ')
           pright(1)=intkey(card,'order-of-leading-right-'//
     1                           'polynomials',0,' ')
           write(iout,2)  nmax(1), npt(1), left, right, 
     1                    der(1,1), der(2,1) , pleft(1), pright(1)
           if (fix(1)) then
               if(nfix(1).eq.1) then
                  write(iout,3) nfix(1), endpts(1,1)
               endif
               if(nfix(1).eq.2) then
                  write(iout,4) nfix(1), endpts(1,1), endpts(2,1)
               endif 
           endif
      endif
      if ( posinp('$time-polynomials',cpass) ) then
           call cardin(card)
           tmax=fpkey(card,'maximum-propagation-time',5.d0,' ')
           omega=fpkey(card,'electric-field-frequency',1.d0,' ')
           left=0.d0
           right=tmax 
           nmax(2)=intkey(card,'order-of-polynomials',10,' ')
           npt(2)=intkey(card,'number-of-points',nmax(2),' ')
           fix(2)=logkey(card,'fix-end-points',.false.,' ')
           decomp(2)=logkey(card,'expand',.false.,' ')
           der(1,2)=fpkey(card,'left-derivative',0.d0,' ')
           der(2,2)=fpkey(card,'right-derivative',0.d0,' ')
           if (fix(2)) then
               nfix(2)=intkey(card,'number-of-fixed-points',
     1                        2,' ')
               endpts(1,2)=0.d0
               endpts(2,2)=tmax
               call fparr(card,'end-points',endpts(1,2),
     1                    nfix(2),' ')
           endif
           pleft(2)=intkey(card,'order-of-leading-left-'//
     1                          'polynomials',1,' ')
           pright(2)=intkey(card,'order-of-leading-right-'//
     1                           'polynomials',0,' ')
           write(iout,5)  nmax(2), npt(2), left, right, 
     1                    der(1,2), der(2,2) ,pleft(2), pright(2)
           if (fix(2)) then
               if(nfix(2).eq.1) then
                  write(iout,6) nfix(2), endpts(1,2)
               endif
               if(nfix(2).eq.2) then
                  write(iout,7) nfix(2), endpts(1,2), endpts(2,2)
               endif 
           endif
      endif
      dim(1)=max(nmax(1),npt(1))
      dim(2)=max(nmax(2),npt(2))
      maxd=max(dim(1),dim(2))
      nply(1)=nmax(2)
      nply(2)=nmax(1)
      n2d=nply(1)*nply(2) 
      if(itsolv) then
         lenbuf=intkey(ops,'linslv=buffer',min(1000000,4*n2d*n2d),' ')
         nrhs=intkey(ops,'linslv=number-of-right-hand-sides',1,' ')
         thresh=fpkey(ops,'linslv=tolerance',1.0d-06,' ')
         cnverg=fpkey(ops,'linslv=convergence',1.d-08,' ')
         nattim=intkey(ops,'linslv=number-of-solutions-at-a-time',
     1                 nrhs,' ')
         maxvec=intkey(ops,'linslv=maximum-number-of-vectors',
     1                 n2d,' ')
         niter=intkey(ops,'linslv=maximum-number-of-iterations',
     1                n2d,' ')
         prgues=logkey(ops,'print=linslv=guess',.false.,' ')
         prvec=logkey(ops,'print=linslv=final-solutions',.false.,' ')
         prres=logkey(ops,'print=linslv=residuals',.false.,' ')
         maxvec=min(maxvec,n2d)
         niter=min(niter,n2d)
         write(iout,8) nrhs, nattim, lenbuf, thresh, cnverg, maxvec,
     1                 niter
         call iosys ('read character "hamiltonian filename" from rwf',
     1                -1,0,0,filham)
         call iosys('open ham as new',0,0,0,filham)
      endif
      ioff=1
      do 10 i=1,2
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
         q(2)=eigc(1)+nmax(1)
         wt(2)=q(2)+dim(2)
         a(2)=wt(2)+dim(2)
         b(2)=a(2)+nmax(2)+1
         p(2)=b(2)+nmax(2)+1
         dp(2)=p(2)+nmax(2)*dim(2)
         ddp(2)=dp(2)+nmax(2)*dim(2)
         pn(2)=ddp(2)+nmax(2)*dim(2)
         dpn(2)=pn(2)+nmax(2)*dim(2)
         ddpn(2)=dpn(2)+nmax(2)*dim(2)
         eigc(2)=ddpn(2)+nmax(2)*dim(2)
         f=eigc(2)+nmax(2)
         df=f+2*maxd
         ddf=df+2*maxd
         driver=ddf+2*maxd
         ind=wpadti(driver+2*maxd)
         ipvt=ind+4*nmax(1)*nmax(2)
         ham(1)=iadtwp(ipvt+maxd) 
         eig(1)=ham(1)+maxd*maxd
         v=eig(1)+maxd
         ham(2)=v+maxd*maxd
         eig(2)=ham(2)+2*maxd*maxd
         work(1)=eig(2)+2*maxd
         work(2)=work(1)+2*maxd*maxd
         work(3)=work(2)+2*maxd*maxd
         work(4)=work(3)+2*maxd*maxd
         wds0=work(4)+2*maxd*maxd
         if(.not.itsolv) then
             hamc=wds0
             rhsc=hamc+2*nmax(1)*nmax(1)*nmax(2)*nmax(2)
             wds0=rhsc+2*nmax(1)*nmax(2)
         else
             diag=wds0
             hbuf=diag+2*n2d
             ibuf=wpadti(hbuf+2*lenbuf)
             rhsc=iadtwp(ibuf+2*lenbuf)
             wds0=rhsc+2*n2d
         endif
         words=wds0
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
   10 continue
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      do 20 i=1,2
         temp(1)=-1.d0
         temp(2)=1.d0
         if(nfix(1).eq.1) then
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
             title='a coefficients'
             call prntrm(title,z(a(i)+1),nmax(i),1,nmax(i),1,iout)
             title='b coefficients'
             call prntrm(title,z(b(i)+1),nmax(i)-1,1,nmax(i)-1,1,iout)
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
     2              z(eigc(i)),z(work(1)),z(work(2)),nmax(i),
     3              npt(i),prnh(i))
         if (prnply(i)) then
             title='polynomials'
             call prntrm(title,z(p(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='first derivative of polynomials'
             call prntrm(title,z(dp(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='second derivative of polynomials'
             call prntrm(title,z(ddp(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='dvr polynomials'
             call prntrm(title,z(pn(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='first derivative of dvr polynomials'
             call prntrm(title,z(dpn(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
             title='second derivative of dvr polynomials'
             call prntrm(title,z(ddpn(i)),npt(i),nmax(i),npt(i),
     1                   nmax(i),iout)
         endif
         if(check(i)) then
            call chk(z(p(i)),z(wt(i)),z(work(1)),z(work(2)),
     1               nmax(i),npt(i))
            call chk(z(pn(i)),z(wt(i)),z(work(1)),z(work(2)),
     1               nmax(i),npt(i))
         endif   
         if(decomp(i)) then
            call fcoef(z(f),z(df),z(ddf),z(pn(i)),z(dpn(i)),z(ddpn(i)),
     1                 z(q(i)),z(wt(i)),z(work(1)),z(work(2)),nmax(i),
     2                 npt(i))
         endif 
 20   continue
      if(calhmo) then
         call ham0(z(pn(1)),z(dpn(1)),z(ddpn(1)),z(wt(1)),z(eigc(1)),
     1             z(ham(1)),z(v),nmax(1),npt(1),pottyp,prnhmo)
         call rdiag(z(v),z(eig(1)),z(work(1)),nmax(1))
         call iosys('write real "H0 eigenvalues" to lamdat',nmax(1),
     1               z(eig(1)),0,' ')
         call iosys('write real "H0 eigenfunctions" to lamdat',
     1               nmax(1)*nmax(1),z(v),0,' ')
         title='eigenvalues of unperturbed hamiltonian'
         call prntrm(title,z(eig(1)),nmax(1),1,nmax(1),1,iout)
         if (prnhmo) then
             title='eigenvectors of unperturbed hamiltonian'
             call prntrm(title,z(v),nmax(1),nmax(1),nmax(1),
     1                   nmax(1),iout)
         endif
         call mkpsi0(z(v),z(driver),nmax(1),state)
      endif 
      if(caltpr) then
         call cfcoef(z(f),z(df),z(ddf),z(pn(2)),z(dpn(2)),z(ddpn(2)),
     1               z(q(2)),z(wt(2)),z(work(1)),z(work(2)),nmax(2),
     2               npt(2))
         call ham1(z(pn(2)),z(dpn(2)),z(ddpn(2)),z(wt(2)),
     1             z(eigc(2)),z(ham(2)),nmax(2),npt(2),timtyp,prntim)
         call rhs1(z(pn(2)),z(q(2)),z(wt(2)),z(eig(2)),nmax(2),npt(2))
         call cmpslv(z(ham(2)),z(eig(2)),ia(ipvt),nmax(2),1)
         call tstsol(z(pn(2)),z(q(2)),z(eig(2)),nmax(2),npt(2),timtyp)
      endif
c     calculate the space/time matrix elements of the harmonic oscillator 
c     plus external field.
      call setind(ia(ind),nply,n2d,2)
      if(.not.itsolv) then
          call ehamxt(z(ham(1)),z(ham(2)),z(eigc(1)),z(eigc(2)),
     1                z(driver),z(hamc),z(rhsc),hbar,omegat,omega,
     2                ia(ind),nply,n2d,prnch)
      else
          call ihamxt(z(ham(1)),z(ham(2)),z(eigc(1)),z(eigc(2)),
     1                z(driver),z(hbuf),ia(ibuf),ia(ind),z(diag),
     2                z(rhsc),hbar,omegat,omega,lenbuf,nply,2,
     3                n2d,pottyp,prnch,ntot,incore)
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'atomic mass              = ',e15.8,/,5x,
     1            'trap frequency           = ',e15.8,/,5x,
     2            'number of atoms in trap  = ',e15.8,/,5x,
     3            'scattering length        = ',e15.8,/,5x,
     4            'trap size                = ',e15.8)
 2    format(/,5x,'order of radial polynomials  = ',i3,/,5x,
     1            'number of integration points = ',i3,/,5x,
     2            'left boundary                = ',e15.8,/,5x,
     3            'right boundary               = ',e15.8,/,5x,
     4            'derivative at left boundary  = ',e15.8,/,5x,
     5            'derivative at right boundary = ',e15.8,/,5x,
     6            'order-of-left-polynomials    = ',i2,/,5x,
     7            'order-of-right-polynomials   = ',i2)
 3    format(/,5x,'number fixed radial end points = ',i1,/,5x,     
     1            'end point                      = ',e15.8)
 4    format(/,5x,'number fixed radial end points = ',i1,/,5x,     
     1            'left end point                 = ',e15.8,/,5x,
     2            'right end point                = ',e15.8)
 5    format(/,5x,'order of time polynomials    = ',i3,/,5x,
     1            'number of integration points = ',i3,/,5x,
     2            'left boundary                = ',e15.8,/,5x,
     3            'right boundary               = ',e15.8,/,5x,
     4            'derivative at left boundary  = ',e15.8,/,5x,
     5            'derivative at right boundary = ',e15.8,/,5x,
     6            'order-of-left-polynomials    = ',i2,/,5x,
     7            'order-of-right-polynomials   = ',i2)
 6    format(/,5x,'number fixed time end points = ',i1,/,5x,     
     1            'end point                      = ',e15.8)
 7    format(/,5x,'number fixed time end points = ',i1,/,5x,     
     1            'left end point               = ',e15.8,/,5x,
     2            'right end point              = ',e15.8)
 8    format(/,5x,'number of right hand sides         = ',i3,/,5x,
     1            'number solved at a time            = ',i3,/,5x,
     2            'buffer length                      = ',i8,/,5x,
     3            'overlap tolerance                  = ',e15.8,/,5x,
     4            'convergence criterion              = ',e15.8,/,5x,
     5            'maximum size of projected subspace = ',i6,/,5x,
     6            'maximum number of iterations       = ',i6)
      end




