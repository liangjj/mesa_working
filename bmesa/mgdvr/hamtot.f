*deck hamtot.f
c***begin prologue     hamtot
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            packed hamiltonian
c***                   
c***description        space hamiltonian in packed form.  
c***references         
c
c***routines called    
c***end prologue       hamtot
      subroutine hamtot(qpnt,hampnt,indpnt,vpnt,dim,n,m,nreg,nfun,
     1                  maxreg,ngot,nonzro,len,addv,filnme,prn,
     2                  pack,sym,coord,region)
      implicit integer (a-z)
      integer*8 qpnt, hampnt, vpnt, pz, indpnt
      logical addv, prn, incore, pack
      character*80 title
      character*(*) filnme, sym, coord, region
      dimension nreg(3), nfun(maxreg,3)
      dimension qpnt(3), hampnt(4), vpnt(4), pz(4)
      dimension n(3), ngot(10,4), nonzro(4)
      dimension len(4), coord(3)
      dimension prn(*)
      common/io/inp, iout
      pointer (indpnt,index(1))      
c
c     calculate the interaction potential
c
      call vpert(qpnt(1),qpnt(2),qpnt(3),vpnt(1),vpnt(2),vpnt(3),
     1           vpnt(4),addv,m,n,dim,ngot(4,4),prn(1))
c
c     either calculate a packed or full form of the one body
c     matrices.
c
c
c     calculate the one dimensional non-zero matrix elements
c
      if(pack) then
         do 10 i=1,dim
            len(i)=n(i)*n(i)
            call pacone(hampnt(i),vpnt(i),pz(i),n(i),nreg(i),nfun(1,i),
     1                  len(i),nonzro(i),ngot(5,i),addv,
     2                  coord(i),prn(4),region)
            call memory(-ngot(3,i),hampnt(i),idum,'dum',idum)
            call memory(-ngot(4,i),vpnt(i),idum,'dum',idum)
 10      continue 
c  
         if(dim.eq.1) then
            return
         endif
c
         if(dim.eq.2) then
c
c           get the memory for the two dimensional index and then
c           calculate it
c
            call memory(n(1)*n(2),indpnt,ngot(1,4),'index',0)
            call setp2(index,m,n(1),n(2),sym,prn(2)) 
            call pac2(pz(4),index,hampnt(4),pz(1),pz(2),vpnt(4),len(4),
     1                len(1),len(2),nonzro(4),nonzro(1),nonzro(2),
     2                m,n(1),n(2),ngot(2,4),prn(2),sym,region)
c
         elseif(dim.eq.3) then
c
c           get the memory for the three dimensional index and then
c           calculate it
c
            call memory(n(1)*n(2)*n(3),indpnt,ngot(1,4),'index',0)
            call setp3(index,m,n(1),n(2),n(3),prn(2)) 
            call pac3(pz(4),index,hampnt(4),pz(1),pz(2),pz(3),vpnt(4),
     1                len(4),len(1),len(2),len(3),nonzro(4),
     2                nonzro(1),nonzro(2),nonzro(3),m,n(1),n(2),n(3),
     3                ngot(2,4),prn(2))
         endif
         call memory(-ngot(2,4),pz(4),idum,'hambuf',idum)
         call memory(-ngot(4,4),vpnt(4),idum,'vint',idum)
         do 20 i=1,dim
            call memory(-ngot(5,i),pz(i),idum,'bufone',idum)
 20      continue   
      else
         do 30 i=1,dim
            call one(hampnt(i),vpnt(i),n(i),addv,coord)
 30      continue
         if(dim.eq.1) then
            return
         endif
         if(dim.eq.2) then
            call memory(n(1)*n(2),indpnt,ngot(1,4),'index',0)
            call setp2(index,m,n(1),n(2),sym,prn(2)) 
            call ham2d(pz(4),index,vpnt(4),hampnt(1),hampnt(2),
     1                 hampnt(4),len(4),nonzro(4),m,n(1),n(2),
     2                 ngot(2,4),prn(2),sym,region)
         endif
         call memory(-ngot(2,4),pz(4),idum,'hambuf',idum)
         call memory(-ngot(4,4),vpnt(4),idum,'vint',idum)	 
      endif   
      return
      end       


