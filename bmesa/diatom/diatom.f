*deck diatom.f 
c***begin prologue     diatom
c***date written       951230   (yymmdd)
c***revision date               (yymmdd)
c***keywords           diatom
c***author             schneider, b. i.(nsf)
c***source             diatom
c***purpose            solve two-dimensional (r,theta) problem
c***                   in orthogonal polynomials
c***                   hamiltonian is;
c***          
c***                  -1./(2*mu) [ d**2/dr**2 
c                                      + 
c                      1./r**2 { (1.-z*z) * d**2/dz**2 -2*z*d/dz 
c                                   - m**2/(1.-z**2) } ] Psi
c                                       = E Psi
c***                                          
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       diatom
      program diatom
c
      implicit integer (a-z)
      parameter ( maxl=20 )
      character*4096 ops
      character*8 cpass
      character*80 title
      character*800 card
      character*32 chrkey, pottyp, diatyp, kedir, hdir
      character*128 fillam
      character*3 prnopt
      logical posinp, logkey, prncof, prnply, prnwpt, check
      logical prnke, prnh, coord, rbcond
      common z(1)
      dimension ia(1), temp(2), prnopt(10)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, half, fpkey, nrmr, nrmth, rmax, temp
      data half /.5d0/
c
      call drum
      write(iout,*)
      write(iout,*) '          diatom diagonalization'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prncof=logkey(ops,'print=m6255=polynomial-coefficients',
     1              .false.,' ')
      prnply=logkey(ops,'print=m6255=polynomials',.false.,' ')
      prnwpt=logkey(ops,'print=m6255=points/weights',.false.,' ')
      prnke=logkey(ops,'print=m6255=kinetic-energy',.false.,' ')
      prnh=logkey(ops,'print=m6255=hamiltonian',.false.,' ')
      prnopt(1)='off'
      if(prncof) then
         prnopt(1)='on'
      endif
      prnopt(2)='off'
      if(prnply) then
         prnopt(2)='on'
      endif
      prnopt(3)='off'
      if(prnwpt) then
         prnopt(3)='on'
      endif
      prnopt(4)='off'
      if(prnke) then
         prnopt(4)='on'
      endif
      prnopt(5)='off'
      if(prnh) then
         prnopt(1)='on'
      endif                                             
      coord=logkey(ops,'to-x-representation',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      diatyp=chrkey(ops,'diagonalization-procedure','direct',' ')
      kedir=chrkey(ops,'diagonalize-kinetic-energy','no',' ')
      hdir=chrkey(ops,'ci','eigenvalues',' ')
      pottyp=chrkey(card,'potential','none',' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$diatom',cpass) ) then
           call cardin(card)
           pottyp=chrkey(card,'potential','none',' ')
           mval=intkey(card,'m-value',0,' ')
           nr=intkey(card,'maximum-order-of-r-polynomials',10,' ')
           nth=intkey(card,'maximum-order-of-theta-polynomials',10,' ')
           nptr=intkey(card,'maximum-number-of-r-points',nr,' ')
           nptth=intkey(card,'maximum-number-of-theta-points',
     1                  nth,' ')
           rbcond=logkey(card,'zero-derivative',.false.,' ')
           rmax=fpkey(card,'maximum-r-value',10.d0,' ')
      endif
      write(iout,1) nr, nth, nptr, nptth, rmax, pottyp, 
     1              diatyp, hdir, mval 
      if(rbcond) then
         write(iout,2)
      else
         write(iout,3)
      endif
      write(iout,4) ( prnopt(i), i=1,5 )
      if(coord) then  
         write(iout,5)
      endif
      if(check) then
         write(iout,6)
      endif
      nmax=max(nptr,nptth,nr,nth)
      ioff=1
      do 10 i=1,2
         r=ioff
         wtr=r+nptr
         th=wtr+nptr
         wtth=th+nptth
         ar=wtth+nptth
         br=ar+nr+1
         ath=br+nr+1
         bth=ath+nth+1
         pr=bth+nth+1
         dpr=pr+nptr*nr
         ddpr=dpr+nptr*nr
         pth=ddpr+nptr*nr
         dpth=pth+nptth*nth
         ddpth=dpth+nptth*nth
         scr=ddpth++nptth*nth
         dum=scr+nmax*nmax
         wds0=dum+nmax*nmax
         if(coord) then
            pnr=wds0
            dpnr=pnr+nptr*nr
            ddpnr=dpnr+nptr*nr
            pnth=ddpnr+nptr*nr
            dpnth=pnth+nptth*nth
            ddpnth=dpnth+nptth*nth  
            eigxth=ddpnth+nptth*nth
            eigxr=eigxth+nth
            eigkth=eigxr+nr
            eigkr=eigkth+nth
            tr=eigkr+nr
            tth=tr+nr*nr
            work=tth+nth*nth
            dum=work+nmax*nmax
            wds0=dum+nmax*nmax
         endif
         if(diatyp.eq.'direct') then
            ham=wds0
            eig=ham+nr*nth*nr*nth
            workd=eig+nr*nth
            dumd=workd+nr*nth
            wds0=dumd+nr*nth            
         else
c            
         endif
         words=wpadti(wds0)
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
      temp(1)=-1.d0
      temp(2)=1.d0        
      call gaussq('legendre',nptth,0.d0,0.d0,0,temp,z(scr),
     1             z(th),z(wtth))
      call convt(z(th),z(wtth),-1.d0,1.d0,nrmth,nptth,'legendre',
     1           prnwpt)
      call gaussq('legendre',nptr,0.d0,0.d0,2,temp,z(scr),
     1             z(r),z(wtr))
      call convt(z(r),z(wtr),0.d0,rmax,nrmr,nptr,'legendre',
     1           prnwpt)
c
c           calculate the recursion coefficients
c
      if(mod(mval,2).eq.0) then
         call cpoly(z(pth),z(th),z(wtth),z(ath),z(bth),
     1              -1.d0,1.d0,z(scr),nth,nptth,0,0)
      else
         call lpoly(z(pth),z(th),z(wtth),z(ath),z(bth),
     1              -1.d0,1.d0,z(scr),nth,nptth,half,half)
      endif     
      if(rbcond) then
         call cpoly(z(pr),z(r),z(wtr),z(ar),z(br),0.d0,rmax,
     1              z(scr),nr,nptr,1,0)
      else
         call cpoly(z(pr),z(r),z(wtr),z(ar),z(br),0.d0,rmax,
     1              z(scr),nr,nptr,1,1)
      endif
      if (prncof) then
          title='theta a coefficients'
          call prntrm(title,z(ath+1),nth,1,nth,1,iout)
          title='theta b coefficients'
          call prntrm(title,z(bth+1),nth-1,1,nth-1,1,iout)
          title='r a coefficients'
          call prntrm(title,z(ar+1),nr,1,nr,1,iout)
          title='r b coefficients'
          call prntrm(title,z(br+1),nr-1,1,nr-1,1,iout)
      endif
c
c           calculate the polynomials and their first and 
c                        second derivatives.
c
      if(mod(mval,2).eq.0) then 
         call gpoly(z(pth),z(dpth),z(ddpth),z(th),
     1              z(ath),z(bth),-1.d0,1.d0,0,0,nth,
     2              nptth,.false.)
      else              
         call lgpoly(z(pth),z(dpth),z(ddpth),z(th),
     1               z(ath),z(bth),-1.d0,1.d0,half,half,nth,
     2               nptth,.false.)
      endif
      if(rbcond) then
         call gpoly(z(pr),z(dpr),z(ddpr),z(r),
     1              z(ar),z(br),0.d0,rmax,1,0,nr,nptr,.false.)
      else
         call gpoly(z(pr),z(dpr),z(ddpr),z(r),
     1              z(ar),z(br),0.d0,rmax,1,1,nr,nptr,.false.)
      endif
      if (prnply) then
          title='theta polynomials'
          call prntrm(title,z(pth),nptth,nth,nptth,nth,iout)
          title='first derivative of theta polynomials'
          call prntrm(title,z(dpth),nptth,nth,nptth,nth,iout)
          title='second derivative of theta polynomials'
          call prntrm(title,z(ddpth),nptth,nth,nptth,nth,iout)
          title='r polynomials'
          call prntrm(title,z(pr),nptr,nr,nptr,nr,iout)
          title='first derivative of r polynomials'
          call prntrm(title,z(dpr),nptr,nr,nptr,nr,iout)
          title='second derivative of r polynomials'
          call prntrm(title,z(ddpr),nptr,nr,nptr,nr,iout)
      endif
      if(check) then
         call chk(z(pth),z(wtth),z(scr),z(dum),nth,nptth)
         call chk(z(pr),z(wtr),z(scr),z(dum),nr,nptr)
      endif   
      if(coord) then
         call diagx(z(pth),z(dpth),z(ddpth),z(ath+1),z(bth+1),
     1              z(pnth),z(dpnth),z(ddpnth),z(eigxth),
     2              z(work),z(dum),nth,nptth,prnke)
         call diagx(z(pr),z(dpr),z(ddpr),z(ar+1),z(br+1),
     1              z(pnr),z(dpnr),z(ddpnr),z(eigxr),z(work),z(dum),
     2              nr,nptr,prnke)
         call kethe(z(pnth),z(dpnth),z(ddpnth),z(th),z(wtth),
     1              z(tth),z(eigkth),z(work),z(dum),mval,nth,nptth,
     2              kedir,prnke)
         call ker(z(pnr),z(dpnr),z(ddpnr),z(r),z(wtr),z(tr),z(eigkr),
     1            z(work),z(dum),nr,nptr,kedir,prnke)
         if(diatyp.eq.'direct') then
            call ddiag(z(tr),z(tth),z(eigxr),z(eigxth),z(ham),
     1                 z(eig),z(workd),z(dumd),1,nth,nr,nth,nr*nth,
     2                 pottyp,hdir,prnh)
         elseif(diatyp.eq.'iterative') then
            call itdiag
         endif
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'maximum number of r polynomials     = ',i5,/,5x,
     1            'maximum number of theta polynomials = ',i5,/,5x,
     2            'maximum number of r points          = ',i5,/,5x,
     3            'maximum number of theta points      = ',i5,/,5x,
     4            'maximum r value                     = ',e15.8,
     5                                                    /,5x,
     6            'potential type                      = 'a16,/,5x,
     7            'diagonalization procedure           = ',a16,/5x,
     8            'diagonalization information         = ',a16,/5x,
     7            'm value                             = ',i2)
 2    format(/,15x,'zero derivative at rmax for boundary condition')
 3    format(/,15x,'zero function at rmax for boundary condition')
 4    format(/,5x,'print options:'/,25x,'polynomial coefficients = ',a3,
     1                            /,25x,'polynomials             = ',a3,
     2                            /,25x,'points/weights          = ',a3,
     3                            /,25x,'kinetic energy          = ',a3,
     4                            /,25x,'hamiltonian             = ',a3)
 5    format(/,15x,'using co-ordinate representation')
 6    format(/,15x,'check orthonormality of polynomials')
      end





