*deck preph.f
c***begin prologue     preph
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           simle one-dimensional hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       preph
      subroutine preph(ham,diag,v,n)
      implicit integer (a-z)
      real*8 ham, diag, v
      dimension ham(n,n), diag(n), v(n) 
      common/io/inp, iout
      do 10 i=1,n
         diag(i) = ham(i,i) + v(i)
         ham(i,i) = 0.d0
 10   continue
      return
      end       



