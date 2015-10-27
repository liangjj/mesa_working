*deck h2h0.f
c***begin prologue     h2h0
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            expansion of eigenvectors of hamiltonian in terms of
c***                   non-interacting states.
c***                   
c***references         
c
c***routines called    
c***end prologue       h2h0
      subroutine h2h0(pham,sym,n,dim,prn)
      implicit integer (a-z)
      integer*8 pham
      real*8 ham, ov, hamx0, hamy0
      real*8 scr
      character*(*) sym
      logical prn
      dimension pham(4,2), n(4,2), prn(*)
      common/io/inp, iout
      pointer (ph,ham(1))
      pointer (pov,ov(1))
      pointer (phamx0,hamx0(1))
      pointer (phamy0,hamy0(1))
      pointer (pscr,scr(1))
c
c
      ph=pham(dim+1,1)
      phamx0=pham(1,2)
      hmtx0=1
      vx0=hmtx0+n(1,1)*n(1,1)
      eigvcx0=vx0+n(1,1)
      if(dim.gt.1) then
         phamy0=pham(2,2)
         hmty0=1
         vy0=hmty0+n(2,1)*n(2,1)
         eigvcy0=vy0+n(2,1)
      endif
      eigvc=1
      if(dim.eq.1) then
         need=wptoin(n(1,1)*n(1,1))
         call getmem(need,pov,ngot,'ovlp',0)
         call ov1d(ov,ham(eigvc),hamx0(eigvcx0),n(1,1),prn)
      else
c
         need=wpadti(1+n(1,1)*n(2,1)*n(dim+1,1))
         call getmem(need,pov,ngot,'ovlp',0)
         nbig=max(n(1,1),n(2,1))
         need=wpadti(1+nbig*nbig*n(dim+1,1))
         call getmem(need,pscr,njunk,'scr',0)
         if(sym.eq.'unsymmetric') then
         call ovusym(ov,scr,ham(eigvc),hamx0(eigvcx0),hamy0(eigvcy0),
     1               n(1,1),n(2,1),n(dim+1,1),prn)
         elseif(sym.eq.'symmetric') then
            call ovsym(ov,scr,ham(eigvc),hamx0(eigvcx0),n(1,1),
     1                 n(dim+1,1),prn)
         else
            call lnkerr('symmetry error')
         endif
         call getmem(-njunk,pscr,idum,'scr',idum)
      endif
      return      
      end       






