c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Driver for Split Operator Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck Program Drivso.f 
***begin prologue     drivso
***date written       960718   (yymmdd)
***revision date               (yymmdd)
***keywords           time, dvr, orthogonal polynomial
***                   
***author             schneider, b. i.(nsf)
***source             drivso
***purpose            driver for solution of time dependent
***                   schroedinger equation split operator
***description
***references
***routines called    iosys, util and mdutil
***end prologue       drivso

c The wavefunction is propagated for small times using the formula
c\begin{eqnarray}
c    \Psi(x,y,z,t + {\delta}t) =  \nonumber \\ 
c exp( - i V(x,y,z,t) {\delta}t/2 ) exp( -i T(x,y,z) {\delta}t )
c exp( - i V(x,y,z,t) {\delta}t/2 ) \Psi(x,y,z,t) 
c\end{eqnarray}
c where, $V(x,y,z,t)$, is the potential, which is local in 
c coordinate space and $T(x,y,z)$ is the
c multidimensional kinetic energy operator or more generally, a sum of terms in 
c each coordinate.
c\begin{equation}
c         T(x,y,z) = T(x) + T(y) + T(z)
c\end{equation}
c Since $T(x,y,z)$ is separable, the propagation can be done, one coordinate at a time.
c The consequence is that the propagation becomes an order N process.  The prefactor
c is not as small as for an FFT propagation but depending on the representation of
c $T(x,y,z)$, can be comparable.  The issue then becomes the value of N for the
c various methods.

      program drivso
      implicit integer (a-z)
      parameter ( maxreg=5000 )
      character*4096 ops
      character*2 itoc
      character*8 prtflg
      character*8 type
      character*80 cpass, chrkey, qtyp, prnkey, vtyp, units
      character*24 coord, key, typpot, system
      character*800 card
      character*128 filbec, filham
      character*3 trapon, typke
      logical prn, dollar, logkey
      logical itdiag, prdvd
      logical nlse, reuse, genpts, pnch, proj
      real*8 fpkey 
      real*8 energy, edge, scale, omega, width, shift, gamma 
      integer*8 pham
      dimension pham(4)
      dimension nwham(4), nwpot(4), n(4), nphy(4)
      dimension npts(maxreg), edge(maxreg+1)
      dimension coord(4), key(4), prn(40), prdvd(15)
      dimension type(2), typpot(4), vtyp(2)
      common/io/inp, iout      
      common/punch/pun(3)
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)

c         set spatial dimensionality of problem and program options 

      spdim=intkey(ops,'number-of-space-dimensions',1,' ')
      ntreg=intkey(ops,'number-of-time-regions',1,' ') 
      genpts=logkey(ops,'automate-points',.false.,' ')
      nlse=logkey(ops,'non-linear-equations',.false.,' ')
      system=chrkey(ops,'coordinate-system','cartesian',' ')
      typke=chrkey(ops,'kinetic-energy-type','dvr',' ')
      pnch=logkey(ops,'m8001=punch',.false.,' ')
      spac=logkey(ops,'no-spatial-hamiltonian',.false.,' ')
      if(pnch) then
         pun(1)=97
	 pun(2)=98
	 pun(3)=99
         open (unit=pun(1),file='real',access='sequential',
     1         form='formatted',err=1000,status='unknown')
         open (unit=pun(2),file='imaginary',access='sequential',
     1         form='formatted',err=1000,status='unknown')
         open (unit=pun(3),file='absolute',access='sequential',
     1         form='formatted',err=1000,status='unknown')
      endif    	 
      proj=logkey(ops,'m8001=projections',.false.,' ')
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      if(nlse) then
         itsolv=.true.
      endif
      type(1)='new'
      type(2)=chrkey(ops,'open-ham','unknown',' ')
      write(iout,1) spdim
      write(iout,2)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as '//type(1),0,0,0,filbec)
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
      write(iout,3) spac, ntreg

c     get all of the one-dimensional matrices needed to construct
c     the spatial part of the hamiltonian and associated quantities.

      do 30 i=1,spdim
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         do 40 j=1,20
            prn(j)=.false.
 40      continue   
         if(typke.eq.'dvr') then
            call basis(ops,cpass,card,i,pham(i),edge,ngrds,coord(i),
     1                 key(i),typpot(i),n(i),nphy(i),nwham(i),
     2                 reuse,prn)
         elseif(typke.eq.'fd') then
            call fdiff(ops,cpass,card,i,pham(i),edge,coord(i),key(i),
     1                 typpot(i),n(i),nphy(i),nwham(i),reuse,prn) 
         else
            call lnkerr('quit. bad ke type') 
         endif
 30   continue
      call lnkerr('quit')

c     get potential parameters

      call vcple(units,vtyp,scale,omega,width,shift,gamma,prn)
      if( dollar('$time',card,cpass,inp) ) then
          ntreg=intkey(card,'number-of-time-regions',1,' ')
          call fparr(card,'time-points',edge,ntreg+1,' ')
      endif
      call so3d(pham,edge,key,typke,typpot,spac,system,card,
     1          units,vtyp,scale,omega,width,shift,
     2          gamma,spdim,nphy,ntreg,prn(12),proj,pnch)
c      call pltwfn(y(driver),y(psi0),ham1(q(1)),ham2(q(2)),
c     1            ham3(q(3)),ham4(q(dim)),ham1(pq(1)),
c     2            ham2(pq(2)),ham3(pq(3)),ham4(pq(dim)),
c     3            ham1(eigv0(1)),ham2(eigv0(2)),ham3(eigv0(3)),
c     4            ham1(eig0(1)),ham2(eig0(2)),ham3(eig0(3)),
c     5            t0,energy,spdim,ntot,n3d,nc,nphy(1),
c     6            nphy(dim),i0stat,proj,pnch)
      if(pnch) then
         close (unit=pun(1),err=1000)
         close (unit=pun(2),err=1000)
         close (unit=pun(3),err=1000)
      endif       
      call chainx(0)
 1000 call lnkerr('error in file handling')                     
      stop
 1    format(/,20x,'time-dependent basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
 2    format(/,15x,'calculation = solve time-dependent schrodinger'
     1             ' equation')
 3    format(/,5x,'time-dependent data',/,5x,
     1            'no spatial hamiltonian   = ',l1,/,5x,
     2            'number of time intervals = ',i3)
      end









