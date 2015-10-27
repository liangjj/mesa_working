*deck basis.f 
c***begin prologue     basis
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             mgtdvr
c***purpose            driver for dvr basis sets
c
c***description           this routine returns a pointer pham
c***                      which points to the first word address
c**                       of a block of memory containing the following:
c***                      hamiltonian(nphy*nphy)
c***                      potential(nphy)
c***                      kinetic energy matrix(nphy*nphy)
c***                      primitive surface values(2)
c***                      first point in the mesh(1)
c***                      coordinates(nphy) 
c***                      weights(nphy) 
c***                      dvr functions(nphy*nphy) 
c***                      dvr first derivatives(nphy*nphy)
c***                      dvr second derivatives(nphy*nphy)    
c***
c***                   if a diagonalization is performed
c***                      eigenvectors(nphy*nphy) of hamiltonian
c***                      surface values(2*nphy) of hamiltonian
c***                      eigenvalues(nphy*nphy) of hamiltonian
c***                      eigenvectors(nphy*nphy) of kinetic energy matrix 
c***                      surface values(2*nphy) of kinetc energy matrix
c***                      eigenvalues(nphy)  of kinetc energy matrix
c***
c***                   for the time variable the order is;
c***                      hamiltonian(nphy*nphy)  
c***                      first point(1)
c***                      coordinates(nphy)
c***                      weights(nphy)
c***                      dvr functions(nphy*nphy)
c***                      dvr first derivative(nphy*nphy)
c***                      dvr second derivative(nphy*nphy)
c***references       
c
c***routines called    lobatto(dvrlib) 
c***end prologue       basis
      subroutine basis(pham,coord,n,nphy)
c
      implicit integer (a-z)
      parameter ( maxreg=500 )
      character*(*) coord
      character*2 itoc, ichr
      character*80 units, typpot, key
      character*80 chrkey, typwt, typarg, typint, refwt, cpass, prnkey
      character*1600 card
      logical dollar, logkey
      logical prn, reuse, nodiag
      real*8 edge, endpts, fpkey 
      real*8 hbar, massau, lenau, timau, mass, pmass, lenscl, escal
      real*8 massn2p, pi, zero, half
#ifdef DECPOINTER
      integer*8 pham
      integer*8 pjunk
#else
#ifdef SGIPOINTER
      integer*4 pham
      integer*4 pjunk
#endif SGIPOINTER
#endif DECPOINTER
      dimension n(*), nphy(*)
      dimension npt(maxreg), pjunk(maxreg,2), junk(maxreg,2)
      dimension edge(maxreg+1), nrq(maxreg), endpts(2)
      dimension prn(12), prnkey(12)
      dimension drop(2), fix(2)
      common/io/inp, iout      
      data prnkey / 'sector-points','sector-polynomials',
     1              'sector-matrix','global-points',
     2              'global-polynomials','potential',
     3              'global-matrix-elements','hamiltonian',
     4              'eigenvalues','eigenvectors','amplitudes','all' /
      data typwt / 'legendre' / 
      data typarg / 'linear' /
      data typint / 'finite' /
      data refwt / 'legendre' /
      data fix / .true., .true. /
      data zero, half / 0.d0, .5d0 /
      data pi  / 3.141592653589793238462643d+00 /
c     hbar in joule-sec
      data hbar/1.054571596d-34/
      data massau, lenau, timau, pmass / 9.10938188d-31,
     1     5.291772083d-11,
     2     2.418884326d-17,
     3     1.67262158d-27 /
      data massn2p / 1.00137841887d0 /                                                                                                                                                     
c
      write(iout,1) coord
      len=length(coord)
c     read the input by stopping at the keyword in the input stream.

      if ( dollar('$h0('//coord(1:len)//')',card,cpass,inp) ) then

c          set the print variables for printing.

           do 10 i=1,11
              prnkey(i)='print='//prnkey(i)
 10        continue
           prnkey(12)=chrkey(card,'print',prnkey(12),' ')
           if(prnkey(12).eq.'all') then
              call setprn(prn,11)
           else
              call setlog(prn,prnkey,card,11)
           endif
           if(coord(1:1).eq.'t') then

c        this is a time calculation

              write(iout,2)
              call fparr(card,'region-boundaries',edge,2,' ')
              typwt=chrkey(card,'weight-type','legendre',' ')
              endpts(1)=edge(1)
              endpts(2)=edge(2)
              n(1)=intkey(card,'polynomial-order',1,' ')
              reuse=logkey(card,'reuse',.false.,' ')
              write(iout,3) n(1)
              write(iout,4) (edge(i),i=1,2)

c             number of quadrature points is one greater
c             than polynomial order.

              npt(1)=n(1)+1
              call tlobato(pham,edge,typwt,npt(1),ndum,prn)
           else

c             this is a space calculation

              write(iout,5)
              mass=fpkey(card,'mass',massau,' ')
              units=chrkey(card,'units','atomic-units',' ')
              if(units.eq.'atomic-units') then
                 mass=mass/massau
              endif
              key='$v0('//coord(1:len)//')'
              nodiag=logkey(card,'do-not-diagonalize',.false.,' ')
              nreg=intkey(card,'number-of-regions',1,' ')
              write(iout,6) nreg
              call fparr(card,'region-boundaries',edge,nreg+1,' ')
              write(iout,7) (edge(i),i=1,nreg+1)
              call intarr(card,'polynomial-order-per-region',
     1                    n,nreg,' ')
              write(iout,8) (n(i), i=1,nreg)
              do 20 i=1,nreg
                 npt(i)=n(i)+1
 20           continue
              do 30 i=1,nreg
                 nrq(i)=npt(i)
 30           continue   
              typwt=chrkey(card,'weight-type','legendre',' ')
              call intarr(card,'number-of-reference-quadrature-'//
     1                    'points-per-region',nrq,nreg,' ')
              write(iout,9) (nrq(i), i=1,nreg)
              nfix=intkey(card,'number-of-fixed-points',2,' ')
              if(nfix.ne.0) then
                 fix(1)=logkey(card,'left-fixed-point',.true.,' ')
                 fix(2)=logkey(card,'right-fixed-point',.true.,' ')
                 drop(1)=logkey(card,'drop-first-function',.false.,' ')
                 drop(2)=logkey(card,'drop-last-function',.false.,' ')
              endif
              bcl=1
              bcr=1
              if(drop(1)) then
                 bcl=0
              endif
              if(drop(2)) then
                 bcr=0
              endif
              call lobatto(pham,mass,edge,coord(1:len),
     1                     key,typwt,typpot,bcl,bcr,n(i),
     2                     npt,nrq,nreg,nodiag,nword,prn,pjunk,junk)
           endif
      endif
      return
 1    format(/,20x,'coordinate = ',a24)
 2    format(/,20x,'generating dvr representation in time')
 3    format(/,1x,'polynomial order in time = ',i4)
 4    format(/,1x,'time edges = ',/,5(e15.5,1x))
 5    format(/,20x,'generating dvr representation in space')
 6    format(/,1x,'number of regions = ', i3)
 7    format(/,1x,'space edges = ',/,5(e15.5,1x))
 8    format(/,1x,'polynomial order in each region = ',/,5(i4,1x))
 9    format(/,1x,'number of reference quadrature points in each'
     1            ' region = ',/,5(i4,1x))
 11   format(/,1x,'basis size = ',i4)
      end
