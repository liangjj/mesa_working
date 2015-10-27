*deck mkh1d.f
c***begin prologue     mkh1d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian for 1-dimensional spatial hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkh1d
      subroutine mkh1d(h1,vecout,vecin,tmp,n,nt,m)
      implicit integer (a-z)
      real*8 h1, vecout, vecin, tmp
      character*80 title
      dimension h1(n,n), vecin(n,m), vecout(n,nt,m)
      dimension tmp(n,m)
      common/io/inp, iout
      call ebc(tmp,h1,vecin,n,n,m)
      call filvt(vecout,tmp,n,nt,m)
      return
      end       
