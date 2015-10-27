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
      subroutine drvh(pham,n,len,nonzro,nword,pack,
     1                drctv,type,ondisk,sym,dim,wdim,prn)
      implicit integer (a-z)
      integer*8 phone, pham
      logical prn, incore, pack, ondisk
      character*80 title
      character*(*) sym, drctv, type
      dimension n(dim+1), nonzro(dim+1)
      dimension pham(dim+1), phone(4)
      dimension len(dim+1), coord(dim)
      dimension nword(wdim,*), ngot(4)
      dimension prn(*)
      common/io/inp, iout
c
c     either calculate a packed or full form of the one body
c     matrices.
c
      do 10 i=1,dim
         call drvh1(pham(i),phone(i),pham(dim+1),n(i),len(i),
     1              nonzro(i),ngot(i),pack,drctv,prn(4))
 10   continue 
      if(dim.eq.1) then
         return
      elseif(dim.eq.2) then
c         call drvh2(pham,phone,pham(dim+1),n,len,nonzro,
c     1              ngot,pack,drctv,type,ondisk,sym,dim,prn(2))
      elseif(dim.eq.3) then
c         call drvh3(pham,phone,pham(dim+1),n,len,nonzro,
c     1              ngot,pack,drctv,type,ondisk,sym,dim,prn(2))
      else
         call lnkerr('error in dimension')
      endif               
      if(pack) then
c
c
         do 20 i=1,dim
            call getmem(-ngot(i),phone(i),idum,'one',idum)
 20      continue   
      endif
c  
      return
      end       



