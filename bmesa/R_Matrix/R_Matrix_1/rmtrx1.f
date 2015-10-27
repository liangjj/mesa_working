c*deck rmtrx1.f 
c***begin prologue     rmtrx1
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           r-matrix
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            stage one of r-matrix code.  compute
c***                   energy independent quantities. 
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       rmtrx1
      program rmtrx1
c
      implicit integer (a-z)
      parameter ( maxreg=100, wdim=10 ) 
      character*4096 ops
      character*2 itoc, ig
      character*3 oldnew
      character*80 cpass, chrkey, qtyp, drctv, kwd
      character*800 card
      character*24 typpot, sym, precon, coord
      character*128 filham
      character*80 prnkey
      logical dollar, logkey, prn, dvdprt, cgrid, hamd
      logical pack1, ondisk
      logical itsolv, itdiag, onoff, only, h0hov
      real*8 fpkey, edge, energy, rbox 
      real*8 cnverg, thresh, eps, eleft, eright
      integer*8 pham, pjunk
      real*8 mass
      character*4 parity
      logical nodiag
      dimension edge(maxreg+1), npt(maxreg), pjunk(maxreg,2)
      dimension junk(maxreg,2)
      dimension pham(4,2)
      dimension n(4,2), nphy(4,2)
      dimension prn(20), nword(wdim,4,2), lenbuf(4,2), 
     1          nonzro(4,2), kwd(2)
      dimension onoff(3,2), coord(3), prnkey(16)
      common/io/inp, iout      
      data pack1 / .false. /
      data prnkey / 'sector-points','sector-polynomials',
     1              'sector-matrix','global-points',
     2              'global-polynomials','potential',
     3              'regional-matrix-elements','one-body-hamiltonian',
     4              'hamiltonian','eigenvalues',
     5              'eigenvectors','diagonalization',
     6              'hamio','h0-to-h-overlaps',
     7              'r-matrix-amplitudes','eigenvalue' /
      data mass, parity, angmom, nodiag / 1.d0,'none',0,.false. /
      call drum
      write(iout,*)
      write(iout,1)
      call iosys ('read character options from rwf',-1,0,0,ops)
c
c
c
c     read in needed information
c
      dim=intkey(ops,'number-of-dimensions',1,' ')
      sym=chrkey(ops,'symmetry','unsymmetric',' ')
      only=logkey(ops,'only-one-body',.false.,' ')
      drctv=chrkey(ops,'directive','diagonalize',' ')
      h0hov=logkey(ops,'overlaps',.false.,' ')
      pack1=logkey(ops,'pack',pack1,' ')
      write(iout,2) dim, sym, only, drctv, h0hov 
c
      do 10 i=1,8
         prn(i)=logkey(ops,'print=sector='//prnkey(i),.false.,' ')
 10   continue   
      prn(9)=logkey(ops,'print=sector=all',.false.,' ')
      if(prn(9)) then
         call setprn(prn(1),8)
      endif
      do 20 i=9,15      
         prn(i)=logkey(ops,'print=r-matrix-1='//prnkey(i),.false.,' ')
 20   continue   
      prn(16)=logkey(ops,'print=r-matrix-1=all',.false.,' ')
      if(prn(16)) then
         call setprn(prn(9),7)
      endif
      prn(16)=logkey(ops,'print=r-matrix-1='//prnkey(16),.false.,' ')
      call iosys ('read character "hamiltonian filename" from rwf',
     1            -1,0,0,filham)
      oldnew=chrkey(ops,'open-hamiltonian-as','new',' ')
      call iosys('open ham as '//oldnew(1:3),0,0,0,filham)
      nphy(dim+1,1)=1
      nphy(dim+1,2)=1
      do 30 i=1,dim
         onoff(i,1)=.false.
         onoff(i,2)=.false.
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         len=length(coord(i))
         do 40 j=1,2
            if(j.eq.1) then
               kwd(1)='$h('//coord(i)(1:len)//')'
               kwd(2)='$v('//coord(i)(1:len)//')'
            else
               kwd(1)='$h0('//coord(i)(1:len)//')'
               kwd(2)='$v0('//coord(i)(1:len)//')'
            endif
            write(iout,3) kwd
            if ( dollar(kwd,card,cpass,inp) ) then
                 onoff(i,j)=.true.
                 nreg=intkey(card,'number-of-regions',1,' ')
                 call intarr(card,'number-of-points-per-region',
     1                       npt,nreg,' ')
                 call fparr(card,'region-boundaries',edge,nreg+1,' ')
                 qtyp=chrkey(card,'quadrature-type','legendre',' ')
                 bcleft=intkey(card,'left-boundary-condition',0,' ')
                 bcright=intkey(card,'right-boundary-condition',0,' ')
                 write(iout,4) nreg
                 write(iout,5) (npt(k), k=1,nreg)
                 write(iout,6) (edge(k),k=1,nreg+1)
                 write(iout,7) bcleft, bcright
c
                 call lobato(pham(i,j),mass,edge,coord,parity,qtyp,
     1                       kwd(2),typpot,bcleft,bcright,angmom,
     2                       nphy(i,j),npt,npt,nreg,nodiag,
     3                       nword(1,i,j),prn,pjunk,junk)
                 nphy(dim+1,j)=nphy(dim+1,j)*nphy(i,j)
            endif
 40      continue   
c----------------------------------------------------------------------c
c                the pham pointer points to;
c
c      h=1                      hamiltonian
c      vphy=h+nphy(i)*nphy(i)   potential
c      h0=vphy+nphy             unperturbed hamiltonian h0
c      srf=h0+nphy*nphy         primitive surface values
c      pt0=srf+ 2               first mesh point
c      q1=pt0+1                 rest of mesh
c      qwt=q1+nphy              weights
c      pq=qwt+nphy              polynomials
c      dpq=pq+nphy*nphy         derivative of polynomials
c      ddpq=dpq+nphy*nphy       second derivative of polynomials
c      eigv=ddpq+nphy*nphy      eigenvectors of h
c      rgama=eigv+nphy*nphy     surface values of eigenvectors
c      eig=rgama+2*nphy         eigenvalues of h
c      eigv0=eig+nphy           eigenvectors of h0      
c      rgama0=eigv0+nphy*nphy   surface values of eigenvectors of h0
c      eig0=rgama0+2*nphy       eigenvalues of h0
c----------------------------------------------------------------------c
 30   continue
c
c     
c
c      are we using iterative or direct solution methods.
c      if iterative, read in data
c
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      if(itdiag) then
        call dvddat(card,cpass,nphy(dim+1,1),nroots,ntrials,
     1              nattim,cnverg,thresh,niter,nvec,lenbuf(dim+1,1),
     2              dvdprt,cgrid,hamd,n0,filham)
        write(iout,9) nroots, nattim, thresh, cnverg, niter, 
     1                 lenbuf(dim+1,1)
      elseif(itsolv) then
        call lindat(card,cpass,nphy(dim+1,1),cnverg,thresh,eps,precon,
     1              nblck,lenbuf(dim+1,1),dvdprt,filham)
        write(iout,11) thresh, cnverg, precon, nblck, lenbuf(dim+1,1)
      else
        nroots=intkey(ops,'number-of-roots',nphy(dim+1,1),' ')
        nroots=min(nroots,nphy(dim+1,1))
        lenbuf(dim+1,1)=intkey(ops,'hamiltonian-buffer',
     1                         min(1000000,
     2                         nphy(dim+1,1)*nphy(dim+1,1)),' ')
        write(iout,12) nphy(dim+1,1), nroots, lenbuf(dim+1,1)       
      endif
c
      call gtrbox(pham,rbox,nphy)
      call iosys('write real "r-matrix box" to ham',1,rbox,0,' ')
c
c     diagonalize the full hamiltonian
c
      if(only) then
         call chainx(0)               
         stop
      endif
      call drvh(pham,nphy,lenbuf,nonzro,nword,pack1,
     1          drctv,'h',ondisk,sym,dim,wdim,prn(9))
      if(sym.eq.'symmetric') then
         nphy(dim+1,1)=nphy(1,1)*(nphy(1,1)+1)/2
      endif
      call hamio(pham(dim+1,1),'h','output',.true.,
     1           sym,nphy(dim+1,1),ndum,prn(13))
      if(h0hov) then
         call h2h0(pham,sym,nphy,dim,prn(14))
      endif
c
c     calculate projections of r-matrix states on to surface functions.
c
      kwd(1)=chrkey(ops,'type-calculation','r-matrix',' ')
      if(kwd(1).eq.'r-matrix' ) then
         call rgamma(pham,sym,nphy(1,1),dim,prn(15))
      elseif( dollar('$'//kwd(1)(1:10),card,cpass,inp) ) then
              niter=intkey(card,'maximum-number-of-iterations',1,' ')
              cnverg=fpkey(card,'convergence-criterion',1.d-08,' ')
              mrts=intkey(card,'number-of-roots',50,' ')
              step=intkey(card,'number-of-steps',1000,' ')
              eleft=fpkey(card,'lowest-energy',-5.d0,' ')
              eright=fpkey(card,'highest-energy',5.d0,' ')
              call bnd1d(pham,eleft,eright,rbox,cnverg,
     1                   nphy(1,1),mrts,step,niter,prn(16))
      else
              call lnkerr('error in input keyword')
      endif
      do 50 i=1,dim
         if(onoff(i,1)) then
            call getmem(-nword(1,i,1),pham(i,1),idum,coord,idum)
         endif
         if(onoff(i,2)) then 
            call getmem(-nword(1,i,2),pham(i,2),idum,coord,idum)
         endif
 50   continue   
      call chainx(0)               
      stop
 1    format(/,20x,'***** R-Matrix Code using DVR Basis Sets *****',
     1       /,5x,'***** Stage One:Hamiltonian Formation,'
     2             'Diagonalization and R-Matrix Amplitudes *****')      
 2    format(/20x,'basic input information',/,5x,
     1       'dimensions                    = ',i2,/,5x,
     2       'symmetry                      = ',a12,/,5x,
     3       'one body calculation          = ',l1,/,5x,
     4       'diagonalize                   = ',a12,/,5x,
     5       'compute overlaps              = ',l1)
 3    format(/,5x,'data keyword = ',a24,2x,'potential keyword = ',a24)
 4    format(/,1x,'number of regions = ', i3)
 5    format(/,1x,'number of points per region = ',/,5(i4,1x))
 6    format(/,1x,'region edges = ',/,5(e15.5,1x))
 7    format(/,1x,'left boundary condition  = ',i1,/,
     1         1x,'right boundary condition = ',i1)
 8    format(/,1x,'number of functions = ',i4)
 9    format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'number of roots at a time          = ',i3,/,5x
     3             'overlap tolerance                  = ',e15.8,/,5x,
     4             'convergence criterion              = ',e15.8,/,5x,
     5             'maximum number of iterations       = ',i6,/,5x,
     6             'hamiltonian buffer length          = ',i8 )           
 11   format(/,15x,'iterative linear system information',/,/,5x,
     1             'overlap tolerance                  = ',e15.8,/,5x,
     2             'convergence criterion              = ',e15.8,/,5x,
     3             'preconditioning                    = ',a24,/,5x,
     4             'block size                         = ',i5,/,5x,
     5             'hamiltonian buffer length          = ',i8 )           
 12   format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix            = ',i6,/,5x,
     2             'number of roots           = ',i6,/,5x,
     3             'hamiltonian buffer length = ',i8 )           
 13   format(/,15x,'hamiltonian created, diagonalized and '/,15x,
     1             'eigenvalues and eigenvectors put on disk. quit')
      end
