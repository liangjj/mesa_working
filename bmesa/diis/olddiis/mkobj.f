*deck mkobj.f
c***begin prologue     mkobj
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***description          
c***references         
c
c***routines called    
c***end prologue       mkobj
      subroutine mkobj(obj,ham,v1,v2,iter,n,type,prnt)
      implicit integer (a-z)
      real*8 obj, ham, v1, v2
      logical prnt
      character*80 title
      character*3 itoc
      character*(*) type
      dimension obj(*), ham(n,n), v1(n), v2(n)
      common/io/ inp, iout
      if(type.eq.'fock matrix') then
         ij=0
         do 10 i=1,n
            do 20 j=1,i
               ij = ij+1 
               obj(ij) = ham(i,j)  
 20         continue
 10      continue   
         nobj=ij
      elseif(type.eq.'wavefunction') then
         call copy(v1,obj,n)
         nobj=n
      elseif(type.eq.'hartree') then
         call copy(v2,obj,n)
         nobj=n
      endif                            
      if(prnt) then
         title='object vector iteration = '//itoc(iter)
         call prntrm(title,obj,nobj,1,nobj,1,iout)
      endif
      return
      end       

