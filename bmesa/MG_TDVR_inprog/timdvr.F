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
      implicit integer (a-z)
      character*4096 ops
      character*2 itoc
      character*8 prtflg
      character*80 cpass, chrkey, precon, title
      character*24 coord, system
      character*800 card
      character*128 filbec
      character*3 trapon
      logical dollar, logkey
      logical itsolv, itdiag
      logical spac, nlse, reuse, pnch, proj, automte
      real*8 fpkey 
      real*8 thresh, cnverg, eps
      real*8 energy, edge, delta 
      real*8 t_0, t_f
#ifdef DECPOINTER
      integer*8 pham
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham
#endif SGIPOINTER
      dimension pham(4), nphy(4), edge(2)
      dimension coord(4)
      common/io/inp, iout      
      common/punch/pun(4)
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)

c
c         set spatial dimensionality of problem and program options 
c

      spdim=intkey(ops,'number-of-space-dimensions',1,' ')
      ntreg=intkey(ops,'number-of-time-regions',1,' ') 
      nc=intkey(ops,'number-of-channels',1,' ') 
      nlse=logkey(ops,'non-linear-equations',.false.,' ')
      dim=spdim+1
      system=chrkey(ops,'coordinate-system','cartesian',' ')
      pnch=logkey(ops,'m8001=punch',.false.,' ')
      reuse=logkey(ops,'reuse-time-keyword',.false.,' ')
      automte=logkey(ops,'automate-time',.false.,' ')
      t_0=fpkey(ops,'initial-time',0.d0,' ')
      t_f=fpkey(ops,'final-time',0.d0,' ')
      npts = intkey(ops,'number-of-time-points-per-region',50,' ')
      delta=( t_f - t_0 )/ntreg
      if(pnch) then
         pun(1)=97
	 pun(2)=98
	 pun(3)=99
	 pun(4)=100
         title=chrkey(ops,'graph-1','graph-1',' ')
         open (unit=pun(1),file=title,access='sequential',
     1         form='formatted',iostat=iostst,status='unknown')
         title=chrkey(ops,'graph-2','graph-2',' ')
         open (unit=pun(2),file=title,access='sequential',
     1         form='formatted',iostat=iostat,status='unknown')
         title=chrkey(ops,'graph-3','graph-3',' ')
         open (unit=pun(3),file=title,access='sequential',
     1         form='formatted',iostat=iostat,status='unknown')
         title=chrkey(ops,'graph-4','graph-4',' ')
         open (unit=pun(4),file=title,access='sequential',
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
      write(iout,1) spdim
      write(iout,2)
      call iosys ('read character "bec filename" from rwf',
     1             -1,0,0,filbec)
      call iosys ('open bec as new',0,0,0,filbec)
      call iosys('rewind all on bec read-and-write',0,0,0,' ')
      write(iout,3) spac, ntreg

c
c     get all of the one-dimensional matrices needed to construct
c     the spatial part of the hamiltonian and associated quantities.
c

      do 10 i=1,spdim
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         call basis(pham(i),coord(i),nphy(i))
 10   continue      

c
c         Begin Time Calculation
c

      edge(1)=t_0
      edge(2)=edge(1)+delta
      reuse=.false.
      do 20 tim=1,ntreg
         coord(dim)=itoc(tim)
         len=length(coord(dim))
         coord(dim)='t'//coord(dim)(1:len)
         len=length(coord(dim))
         if(.not.automte) then

c
c        do the same for the time coordinate as the spatial basis
c

            call basis(pham(dim),coord(dim),nphy(dim))
         else
            call tlobato(pham(dim),edge,'legendre',
     1                   npts,nphy(dim),prn)
            edge(1)=edge(2)
            edge(2)=edge(1)+ delta
         endif
c
c        the procedures to be used depend on whether we are solving
c        a linear or non-linear set of equations.  for the linear case
c        one only needs to solve the set of equations once.  this may
c        be done by a direct of iterative linear system solve.
c        for the non-linear problem, iteration is intrisic and it
c        appears appropriate to use an iterative approach which updates
c        the solution "on the fly".
c

         call psit(pham,system,coord(dim)(1:len),tim,reuse,nphy,spdim,
     #             dim,nc,spac,itsolv,proj,pnch)
         if(automte) then
            reuse=.true.
         endif
 20   continue
      if(pnch) then
         close (unit=pun(1),iostat=iostat)
         close (unit=pun(2),iostat=iostat)
         close (unit=pun(3),iostat=iostat)
         close (unit=pun(4),iostat=iostat)
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









