*deck m6031
c***begin prologue     m6031
c***date written       940101   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6031 ,link 6031
c***author             schneider, barry (nsf)
c***source             m6031
c***purpose            solve one-dimensional scattering problem using
c***                   r-matrix method with analytic basis sets.
c***description        the r-matrix method is used for some simple model
c***                   one-dimensional problems to examine resonances and
c***                   to test convergence.
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6031
      program rmone
      implicit integer (a-z)
      parameter ( nbasis=200, nen=1000 )
      real*8 z, rbox, fpkey, convg, v0, ebeg, deltae
      real*8 rmta, drmta, rmte, drmte, energy
      common a(1)
      dimension z(1)
      dimension power(nbasis), energy(nen)
      equivalence (z,a)
      common /memory/ ioff
      common /io / inp, iout
      logical logkey, exact, calrmt, nospol, noppol
      character *4096 ops
      character *1600 card
      character *16 chrkey, cpass, type
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      exact=logkey(ops,'use-exact-r-matrix',.false.,' ')
      type=chrkey(ops,'model-potential','square-well',' ')
      calrmt=logkey(ops,'calculate-r-matrix',.false.,' ')
      nospol=logkey(ops,'no-s-wave-poles',.false.,' ')
      noppol=logkey(ops,'no-p-wave-poles',.false.,' ')
      call posinp('$rmat',cpass)
      call cardin(card)
      rbox=fpkey(card,'r-matrix-box-size',5.d0,' ')
      nprim=intkey(card,'number-of-primitives',5,' ')
      call intarr(card,'powers-of-polynomials',power,nprim,' ')
      pmax=0
      do 10 i=1,nprim
         if (power(i).eq.0) then
             call lnkerr('quit. a polynomial of zero power entered')
         endif
         pmax=max(pmax,power(i))
 10   continue
      pmax=pmax+pmax+2       
      write(iout,1) nprim,rbox
      call posinp('$energy',cpass)
      call cardin(card)
      noene=intkey(card,'number-of-energies',1,' ')
      ebeg=fpkey(card,'first-energy',0.d0,' ')
      deltae=fpkey(card,'energy-step',.1d0,' ')
      niter=intkey(card,'number-of-iterations',50,' ')
      convg=fpkey(card,'convergence-criterion',1.d-08,' ')
      write(iout,2) noene, niter, convg
      ioff=1
      do 20 i=1,2
         bfbox=ioff
         dbfbox=bfbox+nprim
         hfbox=dbfbox+nprim
         dhfbox=hfbox+nprim
         vec=dhfbox+nprim
         smat=vec+nprim*nprim
         eig=smat+nprim*nprim
         enrm=eig+nprim
         scr=enrm+2*nprim
         ekp=scr+max(nprim*nprim,pmax)
         ekin=ekp+nprim*nprim
         vp=ekin+nprim*nprim
         v=vp+nprim*nprim
         ham=v+nprim*nprim
         words=wpadti(ham+nprim*nprim)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
             call getscm(words,z,ngot,'m6031',0)
             write (iout,*) 'get ',words,' words of core'      
         endif
 20   continue
c
c     orthonormalize basis
c
      call orthbs(z(smat),z(eig),z(vec),z(scr),rbox,power,nprim,ncon)          
c     get basis set values at r-matrix boundary
      call bsabox(z(bfbox),z(dbfbox),z(vec),z(smat),z(scr),rbox,
     1            power,nprim,ncon)
c     kinetic energy matrix elements for ground state
      lval=0
      call kinmat(z(ekp),z(ekin),z(vec),z(smat),z(scr),lval,rbox,
     1            power,nprim,ncon)
c     potential matrix elements
      call vmat(z(vp),z(v),z(vec),z(smat),z(scr),rbox,power,v0,
     1          nprim,ncon,pmax,type)
c     hamiltonian formation
      call vadd(z(ham),z(ekin),z(v),ncon*ncon)
c     diagonalize
      call diag(z(ham),z(bfbox),z(dbfbox),z(hfbox),z(dhfbox),z(eig),
     1          z(scr),ncon,'s-wave')
      if (.not.nospol) then
           call poles(z(eig),z(hfbox),z(enrm),v0,rbox,convg,lval,
     1                niter,ncon,exact)
      endif
      call emake(energy,ebeg,deltae,noene)
      if (calrmt) then
          call rmat(rmta,drmta,rmte,drmte,z(hfbox),z(eig),rbox,
     1              energy,v0,lval,ncon,noene,.true.)
      else
          call match(z(hfbox),z(eig),rbox,energy,lval,v0,ncon,noene)
      endif
c     kinetic energy matrix elements for scattering state
      lval=1
      call kinmat(z(ekp),z(ekin),z(vec),z(smat),z(scr),lval,rbox,
     1            power,nprim,ncon)
c     hamiltonian formation
      call vadd(z(ham),z(ekin),z(v),ncon*ncon)
c     diagonalize
      call diag(z(ham),z(bfbox),z(dbfbox),z(hfbox),z(dhfbox),z(eig),
     1          z(scr),ncon,'p-wave')
      if (.not.noppol) then
           call poles(z(eig),z(hfbox),z(enrm),v0,rbox,convg,lval,
     1                niter,ncon,exact)
      endif
      if (calrmt) then
          call rmat(rmta,drmta,rmte,drmte,z(hfbox),z(eig),rbox,
     1              energy,v0,lval,ncon,noene,.true.)
      else
          call match(z(hfbox),z(eig),rbox,energy,lval,v0,ncon,noene)
      endif
      call iosys ('write integer maxsiz to rwf',1,
     1             canget,0,' ')
      call chainx(0)
 1    format(/,10x,'R-matrix Calculation',//,1x,'no. primitives= ',i4,
     1          5x,'r-matrix box size = ',f10.5)
    2 format(/,1x,'number of energies   = ',i4,
     1         1x,'number of iterations = ',i4,/,1x,
     2            'convergence criterion = ',f10.5)
      stop
      end


