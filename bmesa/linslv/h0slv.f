*deck h0slv.f
c***begin prologue     h0slv
c***date written       970503   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           precondition
c***author             schneider, barry (nsf)
c***source             
c***purpose            create a preconditioned trial vector based
c***                   on a soluble zeroth order hamiltonian.
c***references 
c                      ham(n,n) = lu factors of zeroth order hamiltonian
c                      v(n,n) = diagonal matrix of potential in the dvr
c                             representation.
c                      vecold(n,nvc) = vector from previous iteration in
c                                      dvr representation
c                      vecnew(n,nvc) = new vector in dvr rerpresentation
c***routines called    
c***end prologue       h0slv
      subroutine h0slv(ham,v,vecold,vecnew,ipvt,n,nvc,prnt)
      implicit integer (a-z)
      character*80 title
      logical prnt
      real*8 ham, v, vecold, vecnew
      dimension ham(n,n), v(n,n), vecold(n,nvc), vecnew(n,nvc)
      dimension ipvt(n)
      common/io/inp, iout
      do 10 i=1,n
         do 20 j=1,nvc
            vecnew(i,j) = v(i,i)*vecold(i,j)
   20    continue
   10 continue
      do 30 i=1,nvc
         call sgesl(ham,n,n,ipvt,vecnew(1,i),0)
 30   continue   
      if(prnt) then
         title='new trial vectors'
         call prntrm(title,vecnew,n,nvc,n,nvc,iout)
      endif
      return
      end       
