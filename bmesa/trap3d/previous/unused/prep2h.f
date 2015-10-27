*deck prep2h.f
c***begin prologue     prep2h
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           simle one-dimensional hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       prep2h
      subroutine prep2h(ham01,ham02,diag,v,n01,n02,n)
      implicit integer (a-z)
      real*8 ham01, ham02, diag, v
      dimension ham01(n01,n01), ham02(n02,n02), diag(n), v(n) 
      common/io/inp, iout
      call rzero(diag,n)
      do 10 i=1,n
         i1=ind(i,1)
         j1=ind(i,2)
         diag(i) = diag(i) + ham01(i1,i1) + ham02(j1,j1) + v(i)
 10   continue
      return
      end       



