*deck newsol.f
c***begin prologue     newsol
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           best current solution
c***author             schneider, barry (nsf)
c***source             
c***purpose            express best current davidson solution in original
c***                   basis.
c***                   
c***references         
c
c***routines called    
c***end prologue       newsol
      subroutine newsol(sol,vec,coef,n,m,nrhs,maxvec)
      implicit integer (a-z)
      real*8 sol, vec, coef
      dimension vec(n,maxvec), sol(n,nrhs), coef(maxvec,nrhs) 
      common/io/inp, iout
      call ebcxx(sol,vec,coef,n,m,nrhs,n,n,maxvec)
      return
      end       



