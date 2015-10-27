*deck honv.f
c***begin prologue     honv
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            multiply the hamiltonian matrix on a vector.
c***references         
c
c***routines called    
c***end prologue       honv
      subroutine honv(ham,vecold,vecnew,n,nvc)
      implicit integer (a-z)
      real*8 ham, vecold, vecnew
      dimension ham(n,n), vecold(n,nvc), vecnew(n,nvc)
      common/io/inp, iout
      call ebc(vecnew,ham,vecold,n,n,nvc)
      return
      end       
