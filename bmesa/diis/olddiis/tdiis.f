*deck tdiis.f 
c***begin prologue     tdiis
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            test diis code
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       tprop
      program tdiis
c
      implicit integer (a-z)
      parameter ( number=500 )
      character*4096 ops
      character*8 prtflg
      character*80 cpass, title, chrkey
      character*1600 card
      character*24 pottyp
      character*16 objtyp
      character*128 fillam
      real*8 z, fpkey, left, right, endpts, der, norm0, temp, cnverg
      real*8 pi, hbar, omega, amass, natmin, natmax, u0, tsize
      real*8 scatl, muapp, todiis, suminfo, energy
      logical posinp, logkey, fix, doall, drctv, tfermi, prnt, tcalc
      logical switch, type
      common z(1)
      dimension ia(1), endpts(2), der(2), temp(2), work(2)
      dimension suminfo(number,3), energy(number)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
      call drum
      write(iout,*)
      write(iout,*) '    test iteration/diis          '
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      doall=logkey(ops,'all-tests-on',.false.,' ')
      drctv=logkey(ops,'diis=on',.false.,' ')
      type=logkey(ops,'diis=on=jacobi',.false.,' ')
      prnt=logkey(ops,'print=diis=on',.false.,' ')
      tcalc=logkey(ops,'type-calculation=linear-system',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1             -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$diis',cpass) ) then
           call cardin(card)
           pottyp=chrkey(card,'potential-type','none',' ')
           nmax=intkey(card,'order-of-polynomials',30,' ')
           npt=intkey(card,'number-of-points',nmax+2,' ')
           cnverg=fpkey(card,'convergence',1.d-08,' ')
           maxit=intkey(card,'maximum-number-of-iterations',10,' ')
           trunc=intkey(card,'restart-value',5,' ') 
           keep=intkey(card,'iterations-retained',4,' ')
           todiis=fpkey(card,'switch-to-diis',1.d0,' ')
           switch=logkey(card,'no-iteration',.false.,' ')
           objtyp=chrkey(card,'object-type','fock matrix',' ')
           tfermi=logkey(card,'thomas-fermi',.false.,' ')
           tsize=fpkey(card,'trap-size',5.0d0,' ')
           omega=fpkey(card,'trap-frequency',10.d0,' ')
           omega=omega*2.d0*pi
           if(.not.tcalc) then
              write(iout,1)
c          mass in kilograms of atom.  default is Cs.           
              amass=fpkey(card,'atomic-mass',2.2d-25,' ')
              natmin=fpkey(card,'minimum-number-of-atoms',0.d+00,' ')
              natmax=fpkey(card,'maximum-number-of-atoms',1.d+05,' ')
              nstep=intkey(card,'number-of-steps',10,' ')
c          scattering length in meters                              
              u0=fpkey(card,'potential-strength',2.0d-51,' ')
              scatl=amass*u0/(4.d0*hbar*hbar*pi)
           else
              nrhs=intkey(card,'number-of-right-hand-sides',1,' ')
              nen=intkey(card,'number-of-energies',1,' ')
              call fparr(card,'energies',energy,nen,' ')
              write(iout,2) nrhs, nen
           endif
           fix=logkey(card,'fix-end-points',.false.,' ')
           left=0.d0
           right=tsize
           der(1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2)=fpkey(card,'right-derivative',0.d0,' ')
           if (fix) then
               nfix=intkey(card,'number-of-fixed-points',2,' ')
               endpts(1)=left
               endpts(2)=right
               call fparr(card,'end-points',endpts,nfix,' ')
           endif
           pleft=intkey(card,'order-of-leading-left-polynomials',1,' ')
           pright=intkey(card,'order-of-leading-right-polynomials',
     1                   0,' ')
           write(iout,3)  nmax, npt, left, right, der(1), der(2), 
     1                    pleft, pright, nrhs, maxit, nen
           if (fix) then
               if(nfix.eq.1) then
                  write(iout,4) nfix, endpts(1)
               endif
               if(nfix.eq.2) then
                  write(iout,5) nfix, endpts(1), endpts(2)
               endif 
           endif
      endif
      if(drctv) then
         write(iout,6)
      else
         write(iout,7)
      endif
      maxd=max(nmax,npt)
      ioff=1
      do 10 i=1,2
         q=ioff
         wt=q+maxd
         a=wt+maxd
         b=a+nmax+1
         px=b+nmax+1
         dpx=px+nmax*maxd
         ddpx=dpx+nmax*maxd
         pnx=ddpx+nmax*maxd
         dpnx=pnx+nmax*maxd
         ddpnx=dpnx+nmax*maxd
         pdvr=ddpnx+nmax*maxd
         dpdvr=pdvr+nmax*nmax
         ddpdvr=dpdvr+nmax*nmax
         eigc=ddpdvr+nmax*nmax
         tham=eigc+nmax
         if(.not.tcalc) then
            f=tham+nmax*nmax
            rho=f+nmax*nmax*(maxit+1)
            eig=rho+nmax*nmax
            psi=eig+nmax
            err=psi+nmax
            bb=err+nmax*nmax*(maxit+1)
            bbtmp=bb+(maxit+1)*(maxit+1)
            sol=bbtmp+(maxit+1)*(maxit+1)
            ipvt=wpadti(sol+maxit+1)
         else
            matrix=tham+nmax*nmax
            rhs=matrix+nmax*nmax
            xold=rhs+nmax*nrhs
            err=xold+nmax*(maxit+1)
            xcur=err+nmax*(maxit+1)
            errcur=xcur+nmax
            bb=errcur+nmax
            bbtmp=bb+(maxit+1)*(maxit+1)
            sol=bbtmp+(maxit+1)*(maxit+1)
            ipvt=wpadti(sol+maxit+1)
         endif
         work(1)=iadtwp(ipvt+maxit+1)
         work(2)=work(1)+2*maxd*maxd
         words=wpadti(work(2)+2*maxd*maxd)                       
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'diis',0)
         endif
   10 continue
c
c     in order to proceed it is necessary to be able to calculate integrals
c     between the orthogonal polynomials to be generated on the desired
c     integration interval.  the calculation of the points and weights is
c     done on the interval ( -1.,1. ) and then converted to ( left, right ).
c     the fixed points on the standard interval are the two end points.      
c
      temp(1)=-1.d0
      temp(2)=1.d0
      if(nfix.eq.1) then
         temp(1)=1.d0
      endif
      call gaussq('legendre',npt,0.d0,0.d0,nfix,temp,z(work(1)),
     1             z(q),z(wt))
c     
c           convert to points and weights on (left,right) if needed
c
      call convt(z(q),z(wt),endpts(1),endpts(2),norm0,npt,doall)
c
c           calculate the recursion coefficients
c
      call cpoly(z(px),z(q),z(wt),z(a),z(b),endpts(1),endpts(2),
     1           z(work(1)),nmax,npt,pleft,pright)
      if(doall) then
         title='a coefficients'
         call prntrm(title,z(a+1),nmax,1,nmax,1,iout)
         title='b coefficients'
         call prntrm(title,z(b+1),nmax-1,1,nmax-1,1,iout)
      endif   
c
c           calculate the polynomials and their first and 
c                        second derivatives.
c
      call gpoly(z(px),z(dpx),z(ddpx),z(q),z(a),z(b),endpts(1),
     1           endpts(2),pleft,pright,nmax,npt,doall)
      call diagx(z(px),z(dpx),z(ddpx),z(a+1),z(b+1),z(pnx),
     1           z(dpnx),z(ddpnx),z(eigc),z(work(1)),z(tham),
     2           nmax,npt,doall)
      if(doall) then
         title='polynomials at original points'
         call prntrm(title,z(px),npt,nmax,npt,nmax,iout)
         title='first derivative of polynomials at original points'
         call prntrm(title,z(dpx),npt,nmax,npt,nmax,iout)
         title='second derivative of polynomials at original points'
         call prntrm(title,z(ddpx),npt,nmax,npt,nmax,iout)
         title='dvr polynomials at original points'
         call prntrm(title,z(pnx),npt,nmax,npt,nmax,iout)
         title='first derivative of dvr polynomials at original points'
         call prntrm(title,z(dpnx),npt,nmax,npt,nmax,iout)
         title='second derivative of dvr polynomials at original points'
         call prntrm(title,z(ddpnx),npt,nmax,npt,nmax,iout)
         call chk(z(px),z(wt),z(work(1)),z(work(2)),nmax,npt)
         call chk(z(pnx),z(wt),z(work(1)),z(work(2)),nmax,npt)
         call xmtrx(z(px),z(pnx),z(q),z(wt),z(work(1)),nmax,npt)
      endif   
      call newply(z(px),z(dpx),z(ddpx),z(a),z(b),z(pdvr),z(dpdvr),
     1            z(ddpdvr),z(eigc),z(tham),endpts(1),endpts(2),
     2            pleft,pright,nmax)
      if(doall) then
         title='polynomials at dvr points'
         call prntrm(title,z(px),nmax,nmax,nmax,nmax,iout)
         title='first derivative of polynomials at dvr points'
         call prntrm(title,z(dpx),nmax,nmax,nmax,nmax,iout)
         title='second derivative of polynomials at dvr points'
         call prntrm(title,z(ddpx),nmax,nmax,nmax,nmax,iout)
         title='dvr polynomials at dvr points'
         call prntrm(title,z(pdvr),nmax,nmax,nmax,nmax,iout)
         title='first derivative of dvr polynomials at dvr points'
         call prntrm(title,z(dpdvr),nmax,nmax,nmax,nmax,iout)
         title='second derivative of dvr polynomials at dvr points'
         call prntrm(title,z(ddpdvr),nmax,nmax,nmax,nmax,iout)
      endif
      call ham0(z(pnx),z(dpnx),z(ddpnx),z(wt),z(eigc),z(tham),
     1          omega,nmax,npt,pottyp,tfermi,doall)
      if(.not.tcalc) then
         call nlschr(z(f),z(eig),z(rho),z(psi),z(err),z(pdvr),z(eigc),
     1               z(tham),z(bb),z(bbtmp),z(sol),z(work(1)),
     2               z(work(2)),suminfo,ia(ipvt),natmin,natmax,hbar,
     3               amass,omega,scatl,todiis,switch,cnverg,nmax,
     4               maxit,nstep,number,trunc,prnt,drctv)
      else
         call mkrhs(z(rhs),nmax,nrhs)
         do 20 ien=1,nen
            energy(ien)=energy(ien)*omega
            call linsys(z(tham),z(rhs),energy(ien),z(matrix),z(xold),
     1                  z(err),z(xcur),z(errcur),z(bb),z(bbtmp),z(sol),
     2                  ia(ipvt),todiis,switch,cnverg,nmax,nrhs,maxit,
     3                   trunc,prnt,drctv,type)
 20      continue   
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'extrapolating on NLSE')
 2    format(/,5x,'extrapolating on linear system',/,5x,
     1            'number of right hand sides = ',i3,/,5x,
     2            'number of energies         = ',i3)      
 3    format(/,5x,'order of radial polynomials       = ',i3,/,5x,
     1            'number of integration points      = ',i3,/,5x,
     2            'left boundary                     = ',e15.8,/,5x,
     3            'right boundary                    = ',e15.8,/,5x,
     4            'derivative at left boundary       = ',e15.8,/,5x,
     5            'derivative at right boundary      = ',e15.8,/,5x,
     6            'order-of-left-polynomials         = ',i2,/,5x,
     7            'order-of-right-polynomials        = ',i2,/,5x,
     8            'number of right hand sides        = ',i2,/,5x,
     9            'maximum number of diis iterations = ',i3,/,5x,
     x            'number of energies                = ',i2)
 4    format(/,5x,'number fixed radial end points = ',i1,/,5x,     
     1            'end point                      = ',e15.8)
 5    format(/,5x,'number fixed radial end points = ',i1,/,5x,     
     1            'left end point                 = ',e15.8,/,5x,
     2            'right end point                = ',e15.8)
 6    format(/,15x,'diis procedure is turned on')
 7    format(/,15x,'diis procedure is turned off')
      end

