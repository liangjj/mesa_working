*deck driver.f 
c***begin prologue     driver
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             driver
c***purpose            driver for orthogonal and DVR functions based on
c***                   generalized Gauss quadrature.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       driver
      program driver
c
      implicit integer (a-z)
      character*4096 ops
      character*2 itoc
      character*80 cpass, chrkey, typwt, refwt, testfn
      character*16 typint, parity
      character*24 units, cordsys, coord, tokeep
      logical dollar, logkey, prn, fix, diag, stndrd, drop
      real*8 ply, kmat, edge, alpha, beta, fpkey
      real*8 alf, bet, refalf, refbet, x, scr, ptfix
      real*8 hbar, massau, lenau, timau, mass, pmass, lenscl, escal
      real*8 massn2p, pi, zero, half 
      dimension prn(10), edge(2), fix(2)
      dimension ptfix(2), ngot(3)
      common/io/inp, iout
      data fix/.false.,.false./      
      data ptfix / -1.d0 , 1.d0 /
      data zero, half / 0.d0, .5d0 /
      data pi  / 3.141592653589793238462643d+00 /
c     hbar in joule-sec
      data hbar/1.054571596d-34/
      data massau, lenau, timau, pmass / 9.10938188d-31,
     1     5.291772083d-11,
     2     2.418884326d-17,
     3     1.67262158d-27 /
      data massn2p / 1.00137841887d0 /                
      data cordsys / 'cartesian' /
      pointer (pntply,ply(1))
      pointer (pntkmt,kmat(1))
      pointer (px,x(1))
      pointer (pscr,scr(1))
      pointer (pkplst,kplst(1))
      pointer (pdrlst,drlst(1))
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      prn(1)=logkey(ops,'print=m5000=reference-points-and-weights',
     1              .false.,' ')
      prn(2)=logkey(ops,'print=m5000=ratio-weight',.false.,' ')
      prn(3)=logkey(ops,'print=m5000=lanczos-coefficients',.false.,' ')
      prn(4)=logkey(ops,'print=m5000=points-and-weights',.false.,' ')
      prn(5)=logkey(ops,'print=m5000=orthogonal-polynomials',
     1             .false.,' ')
      prn(6)=logkey(ops,'print=m5000=coordinate-functions',.false.,' ')
      cordsys=chrkey(ops,'coordinate-system',cordsys,' ')
      coord=chrkey(ops,'dimension','x',' ')
      parity=chrkey(ops,'parity','none',' ')
      mass=fpkey(ops,'mass',massau,' ')
      units=chrkey(ops,'units','atomic-units',' ')
      if(units.eq.'atomic-units') then
         mass=mass/massau
      endif
      typint=chrkey(ops,'interval-type','finite',' ')
      typwt=chrkey(ops,'weight-type','legendre',' ')
      if(typwt.eq.'legendre'
     1        .or.
     2   typwt.eq.'hermite'
     3        .or.
     4   typwt.eq.'laguerre'
     5        .or.
     6   typwt.eq.'chebyshev-1'
     7        .or.
     8   typwt.eq.'chebyshev-2'
     9        .or.
     $   typwt.eq.'jacobi') then
         stndrd=.true.
      else
         stndrd=.false.
      endif
      nfix=intkey(ops,'number-of-fixed-points',0,' ')
      diag=logkey(ops,'diagonalize',.false.,' ')
      lensys=length(cordsys)
      lencrd=length(coord)
      write(iout,1) cordsys(1:lensys), coord(1:lencrd), typint, typwt

c     For general weight functions and finite intervals, either
c     [-1,1] or [0,1] works as a standard.  Other finite intervals
c     will convert the standard points and weights to those required.

      if(typint.eq.'finite') then
         call fparr(ops,'end-points',edge,2,' ')
         write(iout,2) (edge(i),i=1,2), nfix
         fix(1)=logkey(ops,'left-fixed-point',.false.,' ')
         fix(2)=logkey(ops,'right-fixed-point',.false.,' ')
         if(nfix.eq.1) then
            if(fix(2)) then
               ptfix(1) = ptfix(2)
            endif
         endif
         if(nfix.ne.0) then 
            write(iout,3) (ptfix(i),i=1,nfix)
         endif
      endif
      if(typwt.eq.'jacobi'.or.
     1   typwt.eq.'laguerre') then
         alpha=fpkey(ops,'alpha',0.d0,' ')
         beta=fpkey(ops,'beta',0.d0,' ')
         write(iout,4) alpha, beta
      elseif(typwt.eq.'rys') then
         alpha=fpkey(ops,'rys-alpha',0.d0,' ')      
         parity='even' 
         write(iout,5) alpha
      endif
      n=intkey(ops,'polynomial-n',10,' ')	      
      refwt=chrkey(ops,'reference-weight','legendre',' ')
      nwt=intkey(ops,'power-of-weight',1,' ')
      nq=intkey(ops,'reference-quadrature-size',n+1,' ')
      npt=intkey(ops,'number-of-points',n+1,' ')          
      testfn=chrkey(ops,'test-function','exponential',' ')
      call getmem(npt,pkplst,ndum,'keep',0)
      drop=logkey(ops,'drop',.false.,' ')
      do 10 i=1,npt
         kplst(i)=i
 10   continue   
      if(.not.drop) then
         ndrop=0
         nkeep=npt
      else
         ndrop=intkey(ops,'number-to-drop',2,' ')
         call getmem(ndrop,pdrlst,mdum,'drop',0)
         call intarr(ops,'drop',drlst,ndrop,' ')
         do 20 i=1,ndrop
            kplst(drlst(i))=0
 20      continue   
         nkeep=0
         do 30 i=1,npt
            if(kplst(i).ne.0) then
               nkeep=nkeep+1
               kplst(nkeep)=i
            endif
 30      continue   
         call getmem(-mdum,pdrlst,idum,'drop',idum)
      endif
      write(iout,6) (kplst(i),i=1,nkeep)
      if(refwt.eq.'jacobi'.or.refwt.eq.'laguerre') then
         refalf=fpkey(ops,'reference-alpha',0.d0,' ')
         refbet=fpkey(ops,'reference-beta',0.d0,' ')
      endif	                   
      write(iout,7) n, nq, refwt, refalf, refbet     
      need=wptoin(npt)
      call getmem(need,px,ngot(1),'x',0)
      wt=1
      junk=wt+npt
      need=wpadti(junk+npt)
      call getmem(need,pscr,njunk,'junk',0)
      call getqpt(x,scr(wt),edge(1),edge(2),refwt,'before',
     1            scr(junk),fix,npt,npt,1,.true.)
      call getmem(-njunk,pscr,idum,'junk',idum)

c
c     calculate the full set of orthogonal and coordinate functions
c     with no regard to necessarily using the entire set.
c

      a=1
      b=a+n+1
      pt=b+n+1
      wpt=pt+n+1
      cp=wpt+n+1
      dcp=cp+(n+1)*npt
      ddcp=dcp+(n+1)*npt
      wtt=ddcp+(n+1)*npt
      dwtt=wtt+npt
      ddwtt=dwtt+npt
      need=wpadti(ddwtt+npt)
      call getmem(need,pntply,ngot(2),'poly',0)  
#ifdef orthlib
      call orpoly(ply(a),ply(b),ply(pt),ply(wpt),
     1            ply(cp),ply(dcp),ply(ddcp),
     2            ply(wtt),ply(dwtt),ply(ddwtt),x,alpha,beta,
     3            refalf,refbet,edge,ptfix,fix,n,nq,npt,nfix,
     4            coord(1:lencrd),stndrd,parity,typwt,refwt,typint,prn)
#else
#ifdef coordlib
      call cordfn(ply(a),ply(b),ply(pt),ply(wpt),
     1            ply(cp),ply(dcp),ply(ddcp),
     2            ply(wtt),ply(dwtt),ply(ddwtt),x,alpha,beta,
     3            refalf,refbet,edge,ptfix,fix,n,nq,npt,nfix,
     4            coord(1:lencrd),stndrd,parity,typwt,refwt,typint,prn)
#endif coordlib
#endif orthlib

      call renorm(ply(cp),ply(dcp),ply(ddcp),ply(wpt),n+1)

c
c     using the new set see how it does on expanding some
c     test functions and their first and second derivatives.
c

      call cmpre(ply(pt),ply(wpt),ply(cp),ply(dcp),ply(ddcp),
     1           testfn,parity,n+1)

c
c     remove the unwanted coordinates and coordinate functions
c

      call basfun(ply(pt),ply(wpt),ply(cp),ply(dcp),ply(ddcp),
     1            ply(wtt),ply(dwtt),ply(ddwtt),kplst,n+1,nkeep,prn(5))


      need=wptoin(nkeep*nkeep)
      call getmem(need,pntkmt,ngot(3),'kmat',0)  
      call drvkmt(ply(pt),ply(wpt),ply(cp),ply(dcp),ply(ddcp),
     1            kmat,cordsys(1:lensys),coord(1:lencrd),mass,edge,
     2            parity,nkeep)
      call drvhmt(ply(pt),kmat,coord(1:lencrd),diag,nkeep)
      call chainx(0)               
      stop
 1    format(/,20x,'orthogonal polynomial basis function code',/,
     1       /,1x,'coordinate system           = ',a16,
     2       /,1x,'coordinate                  = ',a16,
     3       /,1x,'standard-interval type      = ',a16,
     4       /,1x,'weight function             = ',a16)
 2    format(/,1x,'actual left boundary point          = ',e15.8,
     1       /,1x,'actual right boundary point         = ',e15.8,
     2       /,1x,'number of fixed points = ',i1)
 3    format(/,1x,'fixed points = 'e15.8,1x,e15.8)
 4    format(/,1x,'alpha = ',e15.8,
     1       /,1x,'beta = ',e15.8) 
 5    format(/,1x,'rys alpha = ',e15.8)
 6    format(/,1x,'list of retained functions = ',(/,28x,10i4))
 7    format(/,1x,'polynomial n                 = ',i3,
     1       /,1x,'size of reference quadrature = ',i3,
     2       /,1x,'reference weight function    = ',a32,
     3       /,1x,'reference alpha              = ',e15.8,
     4       /,1x,'reference beta               = ',e15.8)
      end






















