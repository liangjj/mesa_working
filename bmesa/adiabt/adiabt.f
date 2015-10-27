*deck adiabt.f 
c***begin prologue     adiabt
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           adiabatic hyperspherical functions
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            dvr expansion of hyperspherical functions.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       adiabt
      program adiabt
c
      implicit integer (a-z)
      parameter ( maxreg=100 , maxdim=6 ) 
      character*4096 ops
      character*80 cpass, qtyp, chrkey
      character*800 card
      character*24 typpot
      character*16 coord, key
      logical dollar, logkey
      logical prn(20)
      real*8 glbg, ham, hamil, edge, fparr 
      real*8 pi, pi2, del
      integer*8 pglbg(maxdim), pham(maxdim), phamil(maxdim)
      dimension pjunk(maxreg,2)
      dimension npt(maxreg), n(4), start(4), end(4), typpot(maxdim)
      dimension edge(maxreg+1), nword(3,4), nphy(4)
      dimension first(4), last(4), junk(maxreg,2), coord(maxdim)
      common/io/inp, iout      
      data pi  / 3.141592653589793238462643d+00 /
      data pi2 / 1.5707963267948966192313215d+00 /
c
      call drum
      write(iout,*)
      write(iout,1)
      call iosys ('read character options from rwf',-1,0,0,ops)
c
      prn(1)=logkey(ops,'print=m7100=sector-points',.false.,' ')
      prn(2)=logkey(ops,'print=m7100=sector-polynomials',.false.,' ')
      prn(3)=logkey(ops,'print=m7100=sector-matrix',.false.,' ')
      prn(4)=logkey(ops,'print=m7100=global-points',.false.,' ')
      prn(5)=logkey(ops,'print=m7100=global-polynomials',.false.,' ')
      prn(6)=logkey(ops,'print=m7100=potential',.false.,' ')
      prn(7)=logkey(ops,'print=m7100=regional-matrix.elements',
     1              .false.,' ')
      prn(8)=logkey(ops,'print=m7100=hamiltonian',.false.,' ')
      prn(9)=logkey(ops,'print=m7100=eigenvalues',.false.,' ')
      prn(10)=logkey(ops,'print=m7100=eigenvectors',.false.,' ')
      prn(11)=logkey(ops,'print=m7100=amplitudes',.false.,' ')      
      prn(12)=logkey(ops,'print=m7100=radial-grid',.false.,' ')
      prn(13)=logkey(ops,'print=m7100=adiabatic-potential',.false.,' ')
      prn(14)=logkey(ops,'print=m7100=adiabatic-wavefunction',
     1               .false.,' ')
      dim=intkey(ops,'number-of-dimensions',1,' ')
c
c     read in needed information
c
c
c     read in the basis set information
c
      call rdlabl(ops,coord,dim)
      do 10 i=1,dim
         qtyp='legendre'
         len=length(coord(i))
         if ( dollar('$dimension('//coord(i)(1:len)//')',card,
     1                 cpass,inp) ) then
              write(iout,2) coord(i) 
              if(coord(i).eq.'hyperangle') then
                 nreg=intkey(card,'number-of-regions',1,' ')
                 call intarr(card,'number-of-points-per-region',
     1                       npt,nreg,' ')
                 bcleft=0
                 bcright=0
                 call fparr(card,'angular-boundaries',edge,nreg+1,' ')
                 do 20 j=1,nreg+1
                    edge(j)=edge(j)*pi2
 20              continue   
              else
                 nreg=intkey(card,'number-of-regions',1,' ')
                 call intarr(card,'number-of-points-per-region',
     1                       npt,nreg,' ')
                 call fparr(card,'region-boundaries',edge,nreg+1,' ')
                 bcleft=intkey(card,'left-boundary-condition',0,' ')
                 bcright=intkey(card,'right-boundary-condition',0,' ')
              endif
              write(iout,3) nreg, (npt(j), j=1,nreg)
              write(iout,4) (edge(j),j=1,nreg+1)
              write(iout,5) bcleft, bcright
              key='$v0('//coord(i)(1:len)//')'
              write(iout,*) key
              call lobatto(pglbg(i),pham(i),phamil(i),edge,
     1                     key,qtyp,typpot(i),bcleft,
     2                     bcright,n(i),npt,nreg,nword(1,i),prn,
     3                     pjunk,junk)
              nphy(i)=n(i)
              start(i)=1
              end(i)=n(i)
              if(bcleft.eq.0) then
                 start(i)=2
                 nphy(i)=nphy(i)-1
              endif
              if(bcright.eq.0) then
                 nphy(i)=nphy(i)-1
                 end(i)=end(i)-1
              endif
         endif
c
c                
 10   continue 
      if(coord(3).eq.'hyperangle') then
c
c        calculate adiabatic potential and the
c              hyperspherical harmonics
c     
         call ylamda(pglbg(3),phamil(3),nphy(3),start(3),prn(12))
      else
         call lnkerr('error in coordinate specification')
      endif
      do 30 i=1,dim
         call memory(-nword(1,i),pglbg(i),idum,'grid',idum)
         call memory(-nword(2,i),pham(i),idum,'ham',idum)
         call memory(-nword(3,i),phamil(i),idum,'hamil',idum)
 30   continue   
      call chainx(0)               
      stop
 1    format(/,20x,'***** Adiabatic Hyperspherical Functions *****')      
 2    format(/,1x,'points, weights and functions for ',a8,' coordinate')
 3    format(/,1x,'number of regions = ',i4,/,1x,
     1            'number of points per region = ',/,(30x,5(i4,1x)))
 4    format(/,1x,'edges = ',/,(9x,5(e15.8,1x)))
 5    format(/,1x,'left boundary condition  = ',i1,/,1x,
     1            'right boundary condition = ',i1)
      end
