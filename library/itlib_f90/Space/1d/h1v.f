*deck h1v.f
c***begin prologue     h1v
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            hamiltonian times space vector
c***                   
c***references         
c
c***routines called    
c***end prologue       h1v
      subroutine h1v(h1,vecout,vecin,n1,nvc)
      implicit integer (a-z)
      real*8 h1, vecout, vecin
      dimension h1(n1,n1), vecout(n1,nvc), vecin(n1,nvc)
      common/io/inp, iout
      call apbc(vecout,h1,vecin,n1,n1,nvc)
      return
      end       

