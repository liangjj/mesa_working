*deck seidel.f
c***begin prologue     seidel
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            gauss-seidel iteration.
c***references         
c
c***routines called    
c***end prologue       seidel
      subroutine seidel(ham,rhs,vecold,vecnew,n,nvc)
      implicit integer (a-z)
      character*80 title
      real*8 ham, rhs, vecold, vecnew, tmp
      logical incore
      dimension ham(n,n), rhs(n), vecold(n,nvc), vecnew(n,nvc)
      common/io/inp, iout
      call copy(vecold,vecnew,n*nvc)
      do 10 k=1,nvc
         do 20 i=1,n
            tmp = rhs(i) + ham(i,i)*vecnew(i,k)
            do 30 j=1,n
               tmp = tmp - ham(i,j)*vecnew(j,k)
 30         continue
            vecnew(i,k) = tmp/ham(i,i)
 20      continue 
 10   continue      
      return
      end       
