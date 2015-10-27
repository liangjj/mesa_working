*deck addvnl.f
c***begin prologue     addvnl
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c
c***routines called    
c***end prologue       addvnl
      subroutine addvnl(vecold,vecnew,vnl,n,nvc)
      implicit integer (a-z)
      real*8 vecold, vecnew, vnl
      dimension vecold(n,nvc), vecnew(n,nvc), vnl(n)
      common/io/inp, iout
      do 10 i=1,n
         do 20 j=1,nvc
            vecnew(i,j) = vecnew(i,j) + vnl(i)*vecold(i,j)
 20      continue
 10   continue   
      return
      end       
