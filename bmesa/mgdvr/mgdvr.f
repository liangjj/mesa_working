*deck mgdvr.f 
c***begin prologue     mgdvr
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           pde, dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             mgdvr
c***purpose            driver for solution of time independent
c***                   schroedinger equation in dvr basis.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       mgdvr
      program mgdvr
c
      implicit integer (a-z)
      parameter ( maxreg=100, maxgrd=30 ) 
      character*4096 ops
      character*2 itoc
      character*80 cpass, chrkey, phrse
      character*24 coord, precon, sym
      character*800 card
      character*8 region
      character*128 filkohn, filham
      logical dollar, logkey
      logical itsolv, itdiag
      logical prnt, dvdprt, pnch
      logical cgrid, hamd, addv, pack1
      real*8 edge, fpkey 
      real*8 thresh, cnverg, eps
      integer*8 qpnt(3), nrmpnt(3), hampnt(4), vpnt(4), indpnt
      dimension nreg(3), coord(3), bcleft(3), bcright(3)
      dimension npt(maxgrd), nfun(maxreg,3), n(3), grdno(maxreg,3)
      dimension edge(maxreg+1,3)
      dimension prnt(30), dvdprt(30)
      dimension ngot(10,4), lenbuf(4), nonzro(4)
      common/io/inp, iout      
      common/punch/pun
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      addv=logkey(ops,'add-one-body-potential',.false.,' ')
      sym=chrkey(ops,'symmetry','unsymmetric',' ')
      pack1=logkey(ops,'pack',.false.,' ')
      prnt(1)=logkey(ops,'print=m6296=input-data',.false.,' ')
      prnt(2)=logkey(ops,'print=m6296=perturbation',.false.,' ')
      prnt(3)=logkey(ops,'print=m6296=hamiltonian',.false.,' ')
      prnt(4)=logkey(ops,'print=m6296=indices',.false.,' ')
      prnt(5)=logkey(ops,'print=m6296=buffer',.false.,' ')
      pnch=logkey(ops,'punch=on',.false.,' ')
      region=chrkey(ops,'region','internal',' ')
      call iosys ('read character "kohn filename" from rwf',-1,0,
     1             0,filkohn)
      call iosys ('open kohn as old',0,0,0,filkohn)
      call iosys('read integer "number of dimensions" from kohn',1,
     1            dim,0,' ')
      call iosys('write character region to kohn',0,0,0,region)
      write(iout,1) dim
c
c     read some basic grid information from file.  this is used for all
c     dimensions and regions. 
c
      call iosys('read integer "no. grids" from kohn',1,ngrid,0,' ')
      call iosys('read integer "no. pts per grid" from kohn',ngrid,
     1            npt,0,' ')
      if(prnt(1)) then
         write(iout,2) ngrid 
         write(iout,3) (npt(i),i=1,ngrid)
      endif
c
c     read in information for each dimension
c
      call iosys('write character symmetry to kohn',0,0,0,sym)
      n3d=1
      do 10 i=1,dim
         call iosys('read character "coordinate type for '//
     1              'dimension '//itoc(i)//'" from kohn',
     2               -1,0,0,coord(i))
         write(iout,4) coord(i)
         call iosys('read integer "number of regions for '//
     1               coord(i)//'" from kohn',1,nreg(i),0,' ')
         call iosys('read integer "grid number for regions '//
     1              'for '//coord(i)//'" from kohn',nreg(i),
     2               grdno(1,i),0,' ')
         call iosys('read real "region edges for '//coord(i)
     1              //'" from kohn',nreg(i)+1,edge(1,i),0,' ')
         call iosys('read integer "left boundary condition '//
     1              'for '//coord(i)//'" from kohn',1,bcleft(i),0,' ')
         call iosys('read integer "right boundary condition '//
     1              'for '//coord(i)//'" from kohn',1,bcright(i),0,' ')
         call iosys('read integer "number of functions per region '//
     1               'for '//coord(i)//'" from kohn',nreg(i),
     2                nfun(1,i),0,' ')     
         call iosys('read integer "size of global basis for '//
     1               coord(i)//'" from kohn',1,n(i),0,' ')
         n3d=n3d*n(i)
         if(prnt(1)) then
            write(iout,5) nreg(i)
            write(iout,6) (grdno(j,i),j=1,nreg(i))
            write(iout,7) (edge(j,i),j=1,nreg(i)+1)
            write(iout,8) bcleft(i)
            write(iout,9) bcright(i)
            write(iout,11) (nfun(j,i),j=1,nreg(i))
         endif
         phrse='read real "scaled points for '//coord(i)//'" from kohn'
         call arayin(phrse,qpnt(i),'real',n(i),ngot(1,i))
         if(prnt(1)) then
            cpass='scaled points'
            call preal(cpass,qpnt(i),n(i),1,iout)
         endif
         phrse='read real "normalization integrals for '//
     1          coord(i)//'" from kohn'
         call arayin(phrse,nrmpnt(i),'real',n(i),ngot(2,i))
         if(prnt(1)) then
            cpass='normalization integrals'
            call preal(cpass,nrmpnt(i),n(i),1,iout)
         endif
         call memory(-ngot(2,i),nrmpnt(i),idum,'dum',idum)
         phrse='read real "unperturbed hamiltonian for '//
     1          coord(i)//'" from kohn'
         call arayin(phrse,hampnt(i),'real',n(i)*n(i),ngot(3,i))
         if(prnt(1)) then
            cpass='kinetic energy'
            call preal(cpass,hampnt(i),n(i),n(i),iout)
         endif
         phrse='read real "one body potential for '//
     1          coord(i)//'" from kohn'
         call arayin(phrse,vpnt(i),'real',n(i),ngot(4,i))
         if(prnt(1)) then
            cpass='potential energy'
            call preal(cpass,vpnt(i),n(i),1,iout)
         endif
c
c        add kinetic and potential energy matrices
c
         if(addv) then
            call addpot(hampnt(i),vpnt(i),n(i))
         endif
 10   continue
c
c      are we using iterative or direct solution methods.
c      if iterative, read in data
c
      if(itdiag) then
        call dvddat(card,cpass,n3d,nroots,ntrials,nattim,cnverg,
     1              thresh,niter,nvec,lenbuf(4),dvdprt,cgrid,hamd,
     2              n0,filham)
        write(iout,12) nroots, nattim, thresh, cnverg, niter, 
     1                 lenbuf(4)
      elseif(itsolv) then
        call lindat(card,cpass,n3d,cnverg,thresh,eps,precon,nblck,
     1              lenbuf(4),dvdprt,filham)
        write(iout,13) thresh, cnverg, precon, nblck, lenbuf(4)
      else
        nroots=intkey(ops,'number-of-roots',n3d,' ')
        nroots=min(nroots,n3d)
        lenbuf(4)=intkey(ops,'hamiltonian-buffer',
     1                        min(1000000,n3d*n3d),' ')
        write(iout,14) n3d, nroots, lenbuf(4)       
        call iosys ('read character "hamiltonian filename" from rwf',
     1               -1,0,0,filham)
        call iosys('open ham as new',0,0,0,filham)
      endif
c
c     set up and diagonalize the matrix
c
      call hamtot(qpnt,hampnt,indpnt,vpnt,dim,n,n3d,nreg,nfun,maxreg,
     1            ngot,nonzro,lenbuf,addv,filham,prnt(2),
     2            pack1,sym,coord,region)
      size=n3d
      if(sym.eq.'symmetric') then
         size=n(1)*(n(1)+1)/2
      endif
      call iosys('write integer "size of full hamiltonian" to kohn',
     1     1,size,0,' ')            
      call chainx(0)               
      stop
 1    format(/,20x,'time-independent basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
 2    format(/,1x,'number of grids = ',i3)
 3    format(/,1x,'no. pts per grid = ',10(i4,1x)) 
 4    format(/,20x,'dimension = ',a24)
 5    format(/,1x,'number of regions = ', i3)
 6    format(/,1x,'grid no. per region = ',/,10(i4,1x))
 7    format(/,1x,'region edges = ',/,5(e15.5,1x))
 8    format(/,1x,'left boundary condition = ',i1)
 9    format(/,1x,'right boundary condition = ',i1)
 11   format(/,1x,'functions per region = ',/,5(i4,1x))
 12   format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'number of roots at a time          = ',i3,/,5x
     3             'overlap tolerance                  = ',e15.8,/,5x,
     4             'convergence criterion              = ',e15.8,/,5x,
     5             'maximum number of iterations       = ',i6,/,5x,
     6             'hamiltonian buffer length          = ',i8 )           
 13   format(/,15x,'iterative linear system information',/,/,5x,
     1             'overlap tolerance                  = ',e15.8,/,5x,
     2             'convergence criterion              = ',e15.8,/,5x,
     3             'preconditioning                    = ',a24,/,5x,
     4             'block size                         = ',i5,/,5x,
     5             'hamiltonian buffer length          = ',i8 )           
 14   format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix            = ',i3,/,5x,
     2             'number of roots           = ',i3,/,5x,
     3             'hamiltonian buffer length = ',i8 )           
      end
