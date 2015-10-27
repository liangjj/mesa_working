*deck hypsph.f 
c***begin prologue     hypsph
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hyperspherical functions
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            solve the 2D hyperspherical equation
c***                   using dvr in angle and radius.
c***                   
c***                         2       2                   2       2
c***              - 1 / 2 [ d G/d(rho) + 1. / (rho*rho)(d G/d(phi) ]
c***                            -(1./8.*rho*rho) G + V G = E G
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       hypsph
      program hypsph
c
      implicit integer (a-z)
      parameter ( maxreg=100, nume=100, wdim=5 ) 
      character*4096 ops
      character*2 itoc
      character*8 dtype
      character*80 cpass, phrse, qtyp, drctv, prntyp, kwd
      character*800 card
      character*24 precon, sym, typpot
      character*16 eunit, dertyp, coord
      character*128 chrkey, filkohn, filham, units
      logical fitprd, ondisk, pack1, cgrid, hamd
      logical dollar, logkey, itsolv, itdiag, dvdprt, prn
c
c      pointers for lobatto routine
c
      integer*8 pgrid(3), pham(3), phamil(3)
      integer*8 pjunk(maxreg,2)
c
c     pointers to grids and potentials for each dimension
c
      integer*8 px(3), pv(3)
      real*8 dscale, pi, pi2, energy, edge, fpkey 
      real*8 thresh, cnverg, eps
      real*8 a, b
      dimension npt(maxreg), n(3), nphy(3), a(2,2), b(2,2)
      dimension coord(3), edge(maxreg+1), kwd(2)
c
c      words returned from memory calls
c
      dimension nword(wdim,3)
c
c
      dimension energy(nume)
      dimension junk(maxreg,2), prntyp(13), prn(20)
      dimension lenbuf(3), nonzro(3)
      common/io/inp, iout      
      data dscale/-.5d0/
      data pi  / 3.141592653589793238462643d+00 /
      data pi2 / 1.5707963267948966192313215d+00 /
      data prntyp / 'adiabatic-potentials','adiabatic-wavefunctions',
     1              'radial-wavefunctions','sector-points',
     2              'sector-polynomials','sector-matrix',
     3              'global-points','global-polynomials','potential',
     4              'regional-matrix-elements','hamiltonian',
     5              'eigenvalues','eigenvectors' /
c
      call drum
      write(iout,*)
      write(iout,1)
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c
      do 10 i=1,3
         prn(i)=logkey(ops,'print=hyperspherical='//prntyp(i),
     1                 .false.,' ')
 10   continue   
      prn(4)=logkey(ops,'print=hyperspherical=all',.false.,' ')     
      if(prn(4)) then
         call setprn(prn(1),3)
      endif
      do 20 i=4,13
         prn(i)=logkey(ops,'print=sector='//prntyp(i),.false.,' ')
 20   continue   
      prn(14)=logkey(ops,'print=sector=all',.false.,' ')
      if(prn(14)) then
         call setprn(prn(4),10)
      endif
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      dim=intkey(ops,'number-of-dimensions',2,' ')
      pack1=logkey(ops,'pack',.false.,' ')
      drctv=chrkey(ops,'directive','form',' ')
      if(sym.eq.'symmetric') then
         pack1=.false.
      endif
c
c     calculate one body operator for each dimension.
c     the assumed form coming from lobatto is
c
c        h =  dscale * <psi(i)| T + L |psi(j)> + v(i)*delta(i,j)
c                     where
c                2    2 
c           T = d / dq
c           L = - [ delta(q-qright) * d/dq - delta(q-qleft) * d/dq ]
c           v = one=body potential
c      
      call rdlabl(ops,coord,2)
      typpot='none'
      nphy(dim+1)=1
      do 30 i=1,2
         len=length(coord(i))
         kwd(1)='$h0('//coord(i)(1:len)//')'
         kwd(2)='$v0('//coord(i)(1:len)//')'
         if ( dollar(kwd(1),card,cpass,inp) ) then
              nreg=intkey(card,'number-of-regions',1,' ')
              call intarr(card,'number-of-points-per-region',
     1                    npt,nreg,' ')
              qtyp=chrkey(card,'quadrature-type','legendre',' ')
              call fparr(card,'region-boundaries',edge,nreg+1,' ')
              bcleft=intkey(card,'left-boundary-condition',0,' ')
              bcright=intkey(card,'right-boundary-condition',0,' ')
              if(coord(i).eq.'hyperangle') then
                 bcleft=0
                 bcright=0
                 call fparr(card,'angular-boundaries',edge,nreg+1,' ')
                 do 40 j=1,nreg+1
                    edge(j)=edge(j)*pi2
 40              continue   
              endif
              write(iout,2) nreg, (npt(j), j=1,nreg)
              write(iout,3) (edge(j),j=1,nreg+1)
              write(iout,4) bcleft, bcright
         endif
         call lobatto(pgrid(i),pham(i),phamil(i),edge,dscale,
     1                kwd(2),qtyp,typpot,bcleft,bcright,
     2                n(i),npt,nreg,nword(1,i),prn(4),pjunk,junk)
         start=1
         nphy(i)=n(i)
         if(bcleft.eq.0) then
            start=2
            nphy(i)=nphy(i)-1
         endif
         if(bcright.eq.0) then
            nphy(i)=nphy(i)-1
         endif
         write(iout,5) nphy(i)
         nphy(dim+1)=nphy(dim+1)*nphy(i)
         call filxv1(pgrid(i),pham(i),px(i),pv(i),start,
     1               n(i),nphy(i),nword(4,i))
         call memory(-nword(1,i),pgrid(i),idum,'grid',idum)
         call memory(-nword(2,i),pham(i),idum,'ham',idum)
 30   continue   
c
c                
c
      kwd(2)='$v(hyperangle)'
      call vparm(typpot,a,b,card,cpass,kwd(2),2)
c
c     get the channel functions
c
      call hypfn(phamil(1),px,pv(3),typpot,a,b,nword(4,3),nphy,prn(9))
      call lnkerr('quit')
c
c      are we using iterative or direct solution methods.
c      if iterative, read in data
c
      if(itdiag) then
        call dvddat(card,cpass,nphy(dim+1),nroots,ntrials,nattim,cnverg,
     1              thresh,niter,nvec,lenbuf(dim+1),dvdprt,cgrid,hamd,
     2              n0,filham)
        write(iout,6) nroots, nattim, thresh, cnverg, niter, 
     1                 lenbuf(dim+1)
      elseif(itsolv) then
        call lindat(card,cpass,nphy(dim+1),cnverg,thresh,eps,precon,
     1              nblck,lenbuf(dim+1),dvdprt,filham)
        write(iout,7) thresh, cnverg, precon, nblck, lenbuf(4)
      else
        nroots=intkey(ops,'number-of-roots',nphy(dim+1),' ')
        nroots=min(nroots,nphy(dim+1))
        lenbuf(dim+1)=intkey(ops,'hamiltonian-buffer',
     1                        min(1000000,nphy(dim+1)*nphy(dim+1)),' ')
        write(iout,8) nphy(dim+1), nroots, lenbuf(dim+1)       
        call iosys ('read character "hamiltonian filename" from rwf',
     1               -1,0,0,filham)
        call iosys('open ham as new',0,0,0,filham)
      endif
c
c     set up the hamiltonian.  the method is very different if
c     a direct construction of the operator is used.
c
      if(.not.itdiag.or..not.itsolv) then
         write(iout,9)
         call drham(phamil,pv(3),px,dscale,nphy,nword(4,1))
      else
         write(iout,11)
      endif
c
c      nen=intkey(ops,'number-of-energies',1,' ')
c      eunit=chrkey(ops,'units','energy',' ')
c      call fparr(ops,'energy',energy,nen,' ')
c      if(eunit.eq.'wave-vector') then
c         do 50 ene=1,nen
c            energy(ene)=energy(ene)*energy(ene)*.5d0
c 50      continue
c      endif   
      do 60 i=1,dim
         call memory(-nword(1,i),pgrid(i),idum,'grid',idum)
         call memory(-nword(2,i),pham(i),idum,'ham',idum)
         call memory(-nword(3,i),phamil(i),idum,'diag',idum)
 60   continue
      call chainx(0)               
      stop
 1    format(/,20x,'***** Hyperspherical Functions *****')      
 2    format(/,1x,'number of regions = ',i4,/,1x,
     1            'number of points per region = ',/,(50x,5(i4,1x)))
 3    format(/,1x,'edges = ',/,(9x,5(e15.8,1x)))
 4    format(/,1x,'left boundary condition  = ',i1,/,1x,
     1            'right boundary condition = ',i1)
 5    format(/,1x,'number of functions = ',i4)
 6    format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'number of roots at a time          = ',i3,/,5x
     3             'overlap tolerance                  = ',e15.8,/,5x,
     4             'convergence criterion              = ',e15.8,/,5x,
     5             'maximum number of iterations       = ',i6,/,5x,
     6             'hamiltonian buffer length          = ',i8 )           
 7    format(/,15x,'iterative linear system information',/,/,5x,
     1             'overlap tolerance                  = ',e15.8,/,5x,
     2             'convergence criterion              = ',e15.8,/,5x,
     3             'preconditioning                    = ',a24,/,5x,
     4             'block size                         = ',i5,/,5x,
     5             'hamiltonian buffer length          = ',i8 )           
 8    format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix            = ',i6,/,5x,
     2             'number of roots           = ',i6,/,5x,
     3             'hamiltonian buffer length = ',i8 )           
 9    format(/,15x,'hamiltonian explicitly constructed')
 11   format(/,15x,'hamiltonian prepared for iterative solve')
      end
















