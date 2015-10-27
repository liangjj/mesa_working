*deck rmtrx2.f 
c***begin prologue     rmtrx2
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           r-matrix
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            1. calculate r-matrices and scattering
c***                      quantities.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       rmtrx2
      program rmtrx2
c
      implicit integer (a-z)
      parameter ( nume=100, wdim=10 ) 
      character*4096 ops
      character*2 itoc, ig
      character*80 cpass, chrkey, prnkey
      character*800 card
      character*24 sym
      character*128 filham
      logical dollar, logkey, prn, smatrix
      integer*8 pham
      real*8 fpkey, energy, rbox, dscale 
      dimension ngot(wdim,4,2), energy(nume)
      dimension prn(6), prnkey(6)
      common/io/inp, iout      
      data dscale / -.5d0 /
      data pack1 / .false. /
      data prnkey / 'r-matrix','solutions','k-matrix',
     1              's-matrix','hamio','matching' /
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
      if(dim.eq.1) then
         sym='symmetric'
      endif
      smatrix=logkey(ops,'s-matrix=on',.false.,' ')
      write(iout,2) dim, sym, smatrix 
c
      do 10 i=1,6
         prn(i)=logkey(ops,'print=r-matrix-2='//prnkey(i),.false.,' ')
 10   continue   
      prn(7)=logkey(ops,'print=r-matrix-2=all',.false.,' ')
      if(prn(7)) then
         call setprn(prn,6)
      endif
      call iosys ('read character "hamiltonian filename" from rwf',
     1            -1,0,0,filham)
      call iosys('open ham as old',0,0,0,filham)
      call hamio(pham,'h','input',.false.,sym,n,nc,prn(5))
c

      call iosys('read real "r-matrix box" from ham',1,rbox,0,' ')
c
      if ( dollar('$scattering',card,cpass,inp) ) then
           call setsct(pham,energy,sym,n,nc,nen,card,dim)
           if(dim.eq.1) then
              call scat1d(pham,energy,rbox,n,nen,smatrix,prn)
           elseif(dim.eq.2) then
c              call scat2d(phamil(1,1),energy,rbox,sym,nphy(1,1),
c     1                    nen,smatrix,dim,prn)
           else
              call lnkerr('dimension error')
           endif
      else
           call lnkerr('error in input keyword')
      endif
      call chainx(0)               
      stop
 1    format(/,20x,'***** R-Matrix Code using DVR Basis Sets *****',
     1       /,20x,'*****       Stage Two: Scattering        *****')
 2    format(/20x,'basic input information',/,5x,
     1       'dimensions                    = ',i2,/,5x,
     2       'symmetry                      = ',a12,/,5x,
     3       'calculate s-matrix            = ',l1)
      end






