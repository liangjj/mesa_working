c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Time Dependent DVR Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck Program Timdvr.f 
c***begin prologue     timdvr
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            driver for solution of time dependent
c***                   schroedinger equation using time based dvr
c***references
c***routines called    iosys, util and mdutil
c***end prologue       timdvr

c The time-dependent Schroedinger is,

c \begin{equation}
c \big( i \hbar { \partial \over \partial t } - 
c            H({\bf r},t) \big) \mid \Psi \> = 0  
c \end{equation}
c Define,
c \begin{equation}
c     \mid \Psi \> = \mid \Psi_{0} \>+ \mid \Phi \> 
c \end{equation}
c then,
c \begin{equation}
c    \big( i \hbar { \partial \over \partial t } 
c                - H({\bf r},t)\big) \mid \Phi \> = \mid U({\bf r},t) \>
c \end{equation}
c where,
c \begin{equation}
c  \mid U({\bf r},t) \> = \big( H({\bf r},t) 
c                              - i \hbar{ \partial \over \partial t }  \big) 
c                                 \mid \Psi_{0} \>
c \end{equation}
c Separate the Hamiltonian into an unperturbed operator, 
c $H_{0}({\bf r},t)$, and a perturbation, $V({\bf r},t)$, where,
c \begin{equation}
c \big( i \hbar { \partial \over \partial t } 
c               - H_{0}({\bf r},t) \big) \mid \Psi_{0} \> 
c                     = 0  
c \end{equation}
c then,
c \begin{equation}
c \mid U({\bf r},t) \> =  V({\bf r},t) \mid \Psi_{0} \>
c \end{equation}
c The perturbation may contain a term which depends on ${ {\mid \Psi \mid}^2 }$
c if we are dealing with the non-linear Schroedinger equation.

      program timdvr
c
      implicit integer (a-z)
      parameter ( maxreg=5000, mgr=10 )
      character*4096 ops
      character*2 itoc
      character*8 prtflg
      character*8 type
      character*80 cpass, chrkey, precon, qtyp, prnkey
      character*24 coord, key, typpot, system
      character*800 card
      character*128 filbec, filham
      character*3 trapon
      logical prn, dollar, logkey
      logical itsolv, itdiag
      logical spac, nlse, reuse, genpts, pnch, proj
      real*8 fpkey 
      real*8 thresh, cnverg, eps
      real*8 energy, edge 
#ifdef DECPOINTER
      integer*8 pham, phamil, pvt
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, phamil, pvt
#endif SGIPOINTER
      dimension phamil(mgr,4), pham(4), ngrds(4)
      dimension nwham(mgr,4), nwpot(mgr,4), ng(mgr,4), nphyg(mgr,4)
      dimension n(4), nphy(4), npts(maxreg), edge(maxreg+1)
      dimension coord(4), key(4), prn(40)
      dimension type(2), typpot(4)
      common/io/inp, iout      
      common/punch/pun(3)
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
c
c         set spatial dimensionality of problem and program options 
c
      spdim=intkey(ops,'number-of-space-dimensions',1,' ')
      ntreg=intkey(ops,'number-of-time-regions',1,' ') 
      genpts=logkey(ops,'automate-points',.false.,' ')
      nc=intkey(ops,'number-of-channels',1,' ') 
      nlse=logkey(ops,'non-linear-equations',.false.,' ')
      dim=spdim+1
      system=chrkey(ops,'coordinate-system','cartesian',' ')
      pnch=logkey(ops,'m8001=punch',.false.,' ')
      if(pnch) then
         pun(1)=97
	 pun(2)=98
	 pun(3)=99
         open (unit=pun(1),file='real',access='sequential',
     1         form='formatted',iostat=iostst,status='unknown')
         open (unit=pun(2),file='imaginary',access='sequential',
     1         form='formatted',iostat=iostat,status='unknown')
         open (unit=pun(3),file='absolute',access='sequential',
     1         form='formatted',iostat=iostat,status='unknown')
         if(iostat.ne.0) then
            call lnkerr('error in file handling')                     
         endif
      endif    	 
      proj=logkey(ops,'m8001=projections',.false.,' ')
      itsolv=logkey(ops,'iterative-linear-system-solve',.false.,' ')
      spac=logkey(ops,'no-spatial-hamiltonian',.false.,' ')
      if(nlse) then
         itsolv=.true.
      endif
      type(1)='new'
      type(2)=chrkey(ops,'open-ham','unknown',' ')
      write(iout,1) spdim
c
      write(iout,2)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as '//type(1),0,0,0,filbec)
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
      write(iout,3) spac, ntreg
c
c     get all of the one-dimensional matrices needed to construct
c     the spatial part of the hamiltonian and associated quantities.
c
      do 30 i=1,spdim
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         ngrds(i)=1
         do 40 j=1,20
            prn(j)=.false.
 40      continue  
         call basis(ops,cpass,card,i,phamil(1,i),edge,ngrds(i),coord(i),
     1              key(i),typpot(i),ng(1,i),nphyg(1,i),nwham(1,i),
     2              reuse,prn)
 30   continue      

c
c
c         Begin Time Calculation
c
c             set linear system procedures
c      
      if(itsolv) then
         if( dollar('$gmres',card,cpass,inp) ) then
             call lindat(card,cpass,cnverg,thresh,eps,precon,nblck,
     1                   prn(21),filham,type(2))
         endif
      endif
      do 50 i=1,spdim
         pham(i)=phamil(1,i)
         n(i)=ng(1,i)
         nphy(i)=nphyg(1,i)
 50   continue
      if(genpts) then
         if( dollar('$time',card,cpass,inp) ) then
            call intarr(card,'number-of-points-per-region',
     1                  npts,ntreg,' ')
            call fparr(card,'time-points',edge,ntreg+1,' ')
         else
            call timpts(edge,npts,ntreg)
         endif
      endif
      reuse=.false.
      do 60 tim=1,ntreg
c

         coord(dim)=itoc(tim)
         len=length(coord(dim))
         coord(dim)='t'//coord(dim)(1:len)
         len=length(coord(dim))
         key(dim)='$v0('//coord(dim)(1:len)//')'  
         ngrds(dim)=1
c
c        do the same for the time coordinate as the spatial basis
c
         do 70 i=1,20
            prn(i)=.false.
 70      continue   
         if(genpts) then
            n(dim)=npts(tim)
         endif
         call basis(ops,cpass,card,dim,pham(dim),edge(tim),ngrds(dim),
     1              coord(dim),key(dim),typpot(dim),n(dim),
     2              nphy(dim),nwham(1,dim),reuse,prn)
c
c        the procedures to be used depend on whether we are solving
c        a linear or non-linear set of equations.  for the linear case
c        one only needs to solve the set of equations once.  this may
c        be done by a direct of iterative linear system solve.
c        for the non-linear problem, iteration is intrisic and it
c        appears appropriate to use an iterative approach which updates
c        the solution "on the fly".
c
         call psit(pham,pvt,spdim,dim,nphy,nc,key,typpot,spac,
     1             tim,system,itsolv,cnverg,thresh,eps,precon,
     2             nblck,nwpot(1,dim),card,prn(12),hdiag,reuse,
     3             proj,pnch)
 60   continue
      if(pnch) then
         close (unit=pun(1),iostat=iostat)
         close (unit=pun(2),iostat=iostat)
         close (unit=pun(3),iostat=iostat)
         if(iostat.ne.0) then
            call lnkerr('error in file handling')
         endif
      endif       
      call chainx(0)
      stop
 1    format(/,20x,'time-dependent basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
 2    format(/,15x,'calculation = solve time-dependent schrodinger'
     1             ' equation')
 3    format(/,5x,'time-dependent data',/,5x,
     1            'no spatial hamiltonian   = ',l1,/,5x,
     2            'number of time intervals = ',i3)
      end









