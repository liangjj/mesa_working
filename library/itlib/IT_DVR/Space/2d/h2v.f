*deck h2v.f
c***begin prologue     h2v
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            space hamiltonian times space*time vector
c***                   
c***references         
c
c***routines called    
c***end prologue       h2v
      subroutine h2v(h1,h2,vecout,vecin,n1,n2,nvc)
      implicit integer (a-z)
      real*8 h1, h2, vecout, vecin
      dimension h1(n1,n1), h2(n2,n2), vecout(n2,n1,nvc)
      dimension vecin(n2,n1,nvc)
      common/io/inp, iout
      call apbc(vecout,h2,vecin,n2,n2,n1*nvc)
      do 10 i=1,nvc
         call apbct(vecout(1,1,i),vecin(1,1,i),h1,n2,n1,n1)
 10   continue
      return
      end       
