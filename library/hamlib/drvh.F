*deck drvh.f
c***begin prologue     drvh
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for dvr hamiltonian formation and 
c***                   diagonalization
c***                   
c***description        space hamiltonian in packed form.  
c***references         
c
c***routines called    
c***end prologue       drvh
      subroutine drvh(pham,px,pv,phamil,n,len,nonzro,nword,pack,
     1                drctv,type,ondisk,sym,dim,wdim,prn)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pham, px, pv, phone, phamil
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, px, pv, phone, phamil
#endif SGIPOINTER
      logical prn, incore, pack, ondisk
      character*80 title
      character*(*) sym, drctv, type
      dimension n(dim+1), nonzro(dim+1)
      dimension pham(dim+1), px(dim), pv(dim+1), phone(4)
      dimension len(dim+1), coord(dim)
      dimension nword(wdim,*), ngot(4)
      dimension prn(*)
      common/io/inp, iout
c
c     either calculate a packed or full form of the one body
c     matrices.
c
      do 10 i=1,dim
         call drvh1(pham(i),pv(i),phone(i),phamil,n(i),len(i),
     1              nonzro(i),ngot(i),pack,drctv,prn(4))
         if(pack) then
c
c           release memory for one body hamiltonian if 1d packing is used.
c
            call getmem(-nword(1,i),pham,idum,'hamil',idum)
c
         endif	       
 10   continue 
      if(dim.eq.1) then
         return
      elseif(dim.eq.2) then
         call drvh2(pham,pv,phone,px,phamil,n,len,nonzro,
     1              ngot,pack,drctv,type,ondisk,sym,dim,prn(2))
      elseif(dim.eq.3) then
         call drvh3(pham,pv,phone,px,phamil,n,len,nonzro,
     1              ngot,pack,drctv,type,ondisk,sym,dim,prn(2))
      else
         call lnkerr('error in dimension')
      endif               
      if(pack) then
c
c        release memory for 1d-grid, 1d-potentials, 1d-buffers,

c        at this point all that is left is the hamiltonian information
c        for the two-body matrices.
c
         do 20 i=1,dim
            call getmem(-nword(2,i),px(i),idum,'x',idum)
            call getmem(-nword(3,i),pv(i),idum,'v',idum)
            call getmem(-ngot(i),phone(i),idum,'one',idum)
 20      continue   
      else
         do 30 i=1,dim
               call getmem(-nword(1,i),pham(i),idum,'hamone',idum)
               call getmem(-nword(2,i),px(i),idum,'x',idum)
               call getmem(-nword(3,i),pv(i),idum,'v',idum)
 30         continue   
      endif
c  
      return
      end       



