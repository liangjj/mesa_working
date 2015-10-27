*deck hypfit.f 
c***begin prologue     hypfit
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hyperspherical functions
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            using dvr functions fit the
c***                   hyperspherical functions.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       hypfit
      program hypfit
c
      implicit integer (a-z)
      parameter ( nume=100, maxreg=50, wdim=5 ) 
      character*4096 ops
      character*1 itoc
      character*2 coord
      character*16 eunit, dertyp
      character*24 typpot, type, key
      character*128 chrkey
      character*80 qtyp
      logical logkey, fitprd
      logical prn(20)
      integer*8 pgrid(2), pham(2), phamil(2), pjunk(maxreg,2)
      real*8 energy, edge, dscale 
      dimension energy(nume)
      dimension nword(wdim,3), n(2), nphy(2)
      dimension npt(maxreg), junk(maxreg,2)
      dimension edge(maxreg+1)
      common/io/inp, iout      
      data dscale/-.5d0/
c
      call drum
      write(iout,*)
      write(iout,1)
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c
      prn(1)=logkey(ops,'print=m6299=hypersphericals',.false.,' ')
      prn(2)=logkey(ops,'print=m6299=variables',.false.,' ')
      prn(3)=logkey(ops,'print=m6299=derivatives',.false.,' ')
      prn(4)=logkey(ops,'print=m6299=radial-grid',.false.,' ')
      prn(5)=logkey(ops,'print=m6299=adiabatic-potential',.false.,' ')
      prn(6)=logkey(ops,'print=m6299=adiabatic-wavefunction',
     1              .false.,' ')
      prn(7)=logkey(ops,'print=m6299=all',.false.,' ')
      if(prn(7)) then
         call setprn(prn(1),6)
      endif
      prn(7)=logkey(ops,'print=sector-points',.false.,' ')
      prn(8)=logkey(ops,'print=sector-polynomials',.false.,' ')
      prn(9)=logkey(ops,'print=sector-matrix',.false.,' ')
      prn(10)=logkey(ops,'print=global-points',.false.,' ')
      prn(11)=logkey(ops,'print=global-polynomials',.false.,' ')
      prn(12)=logkey(ops,'print=potential',.false.,' ')
      prn(13)=logkey(ops,'print=regional-matrix-elements',
     1               .false.,' ')
      prn(14)=logkey(ops,'print=hamiltonian',.false.,' ')
      prn(15)=logkey(ops,'print=eigenvalues',.false.,' ')
      prn(16)=logkey(ops,'print=eigenvectors',.false.,' ')
      prn(17)=logkey(ops,'print=amplitudes',.false.,' ')      
      prn(18)=logkey(ops,'print=lobatto=all',.false.,' ')      
      if(prn(18)) then
         call setprn(prn(7),11)
      endif
      dertyp=chrkey(ops,'type-derivative','hyperspherical',' ')
      mmax=intkey(ops,'maximum-m',0,' ')
      type=chrkey(ops,'function-type','bessel',' ')
      fitprd=logkey(ops,'fit-product',.false.,' ')
c
c     read in needed information
c
      nen=intkey(ops,'number-of-energies',1,' ')
      eunit=chrkey(ops,'units','energy',' ')
      nreg=intkey(ops,'number-of-regions',1,' ')
      call intarr(ops,'number-of-points-per-region',npt,nreg,' ')
      call fparr(ops,'region-boundaries',edge,nreg+1,' ')
      qtyp=chrkey(ops,'quadrature-type','legendre',' ')
      bcleft=intkey(ops,'left-boundary-condition',0,' ')
      bcright=intkey(ops,'right-boundary-condition',0,' ')
      call fparr(ops,'energy',energy,nen,' ')
      do 10 i=1,2
         coord='r'//itoc(i)
         key='$v0('//coord//')'
         call lobatto(pgrid(i),pham(i),phamil(i),edge,
     1                dscale,key,qtyp,typpot,bcleft,bcright,
     2                n(i),npt,nreg,nword(1,i),prn(7),pjunk,junk)
         call memory(-nword(2,i),pham(i),idum,'ham',idum)
         call memory(-nword(3,i),phamil(i),idum,'hamil',idum)
 10   continue   
      if(eunit.eq.'wave-vector') then
         do 20 ene=1,nen
            energy(ene)=energy(ene)*energy(ene)*.5d0
 20      continue
      endif   
      call drvhyp(pgrid,energy,mmax,n,nen,type,
     1            dertyp,fitprd,prn)
      call chainx(0)               
      stop
 1    format(/,20x,'***** Hyperspherical Functions *****')      
      end








