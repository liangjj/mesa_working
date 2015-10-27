c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{MGDVR: Driver for Schroedinger Equation in DVR Representation}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

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
      parameter (  wdim=5 ) 
      character*4096 ops
      character*2 itoc
      character*80 cpass, chrkey, phrse, drctv
      character*24 coord, precon, sym
      character*800 card
      character*128 filkohn, filham
      character*8 dtype
      logical dollar, logkey
      logical itsolv, itdiag
      logical prn, dvdprt, pnch
      logical cgrid, hamd, pack1, ondisk
      real*8 thresh, cnverg
#ifdef DECPOINTER
      integer*8 pham, px, pv, phtot
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, px, pv, phtot
#endif SGIPOINTER
      dimension n(4), nphy(4), coord(4)
      dimension prn(30), dvdprt(30)
      dimension lenbuf(4), nonzro(4)
      dimension pham(4), px(4), pv(4), ngrds(4)
      dimension nword(10,wdim)
      common/io/inp, iout      
      common/punch/pun
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      itdiag=logkey(ops,'iterative-diagonalization',.false.,' ')
      sym=chrkey(ops,'symmetry','unsymmetric',' ')
      dim=intkey(ops,'number-of-space-dimensions',1,' ')
      pack1=logkey(ops,'pack',.false.,' ')
      dtype=chrkey(ops,'open-as','unknown',' ')
      drctv=chrkey(ops,'directive','diagonalize',' ')
      ondisk=logkey(ops,'buffer-on-disk',.false.,' ')
      if(sym.eq.'symmetric') then
         pack1=.false.
      endif
      pnch=logkey(ops,'punch=on',.false.,' ')
      call iosys ('read character "kohn filename" from rwf',-1,0,
     1             0,filkohn)
      call iosys ('open kohn as '//dtype,0,0,0,filkohn)
      write(iout,1) dim
c
c     read some basic grid information from file.  this is used for all
c     dimensions and regions. 
c
c
c     read in information for each dimension
c
      call iosys('write character symmetry to kohn',0,0,0,sym)
      nphy(dim+1)=1
      do 10 i=1,dim
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         call basis(pham(i),coord(i),nphy(i))
         nphy(dim+1)=nphy(dim+1)*nphy(i) 
         call filxv1(pham(i),px(i),pv(i),nphy(i),nword(1,i))
 10   continue
c
c      are we using iterative or direct solution methods.
c      if iterative, read in data
c
      lenbuf(dim+1)=intkey(ops,'hamiltonian-buffer',
     1                     min(1000000,nphy(dim+1)*nphy(dim+1)),' ')
      call iosys ('read character "hamiltonian filename" from rwf',
     1             -1,0,0,filham)
      call iosys('open ham as new',0,0,0,filham)
      if(itdiag) then
        call dvddat(card,cpass,cnverg,thresh,precon,nphy(dim+1),nroots,
     1              nattim,ntrials,maxit,maxvec,nblck,dvdprt)
        call prepit(pham,px,pv,dim,nphy(dim+1),nphy,nroots,
     1              ntrials,nattim,cnverg,thresh,precon,
     1              maxit,maxvec,nblck,dvdprt)
      elseif(itsolv) then
        call lindat(card,cpass,cnverg,thresh,precon,nphy(dim+1),
     1              maxit,maxvec,nblck,dvdprt)
      else
        nroots=intkey(ops,'number-of-roots',nphy(dim+1),' ')
        nroots=min(nroots,nphy(dim+1))
        write(iout,2) nphy(dim+1), nroots, lenbuf(dim+1)       
c
c       set up and diagonalize the matrix
c
        call drvh(pham,px,pv,phtot,nphy,lenbuf,nonzro,nword,pack1,
     1            drctv,'h',ondisk,sym,dim,wdim,prn)
        size=nphy(dim+1)
        if(sym.eq.'symmetric') then
           size=nphy(1)*(nphy(1)+1)/2
        endif
        call iosys('write integer "size of full hamiltonian" to kohn',
     1              1,size,0,' ')            
      endif
      call chainx(0)               
      stop
 1    format(/,20x,'time-independent basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
 2    format(/,15x,'direct diagonalization',/,/,5x,
     1             'size of matrix            = ',i6,/,5x,
     2             'number of roots           = ',i6,/,5x,
     3             'hamiltonian buffer length = ',i8 )           
      end

