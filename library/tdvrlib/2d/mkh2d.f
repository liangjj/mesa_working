*deck mkh2d.f
c***begin prologue     mkh2d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian for 2-dimensional spatial hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkh2d
      subroutine mkh2d(h1,h2,vecout,vecin,tmp,nd,nt,m)
      implicit integer (a-z)
      real*8 h1, h2, vecout, vecin, tmp
      dimension nd(*)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2))
      dimension vecout(nd(2),nd(1),nt,m)
      dimension vecin(nd(2),nd(1),m)
      dimension tmp(nd(2),nd(1),m)
      common/io/inp, iout
      call ebc(tmp,h2,vecin,nd(2),nd(2),nd(1)*m)
      call filvt(vecout,tmp,nd(1)*nd(2),nt,m)
      do 10 i=1,m
         call ebc(tmp(1,1,i),vecin(1,1,i),h1,nd(2),nd(1),nd(1))
 10   continue   
      call filvt(vecout,tmp,nd(1)*nd(2),nt,m)
      return
      end       


