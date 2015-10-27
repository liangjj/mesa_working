*deck h3v.f
c***begin prologue     h3v
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
c***end prologue       h3v
      subroutine h3v(h1,h2,h3,vecout,vecin,n1,n2,n3,nvc)
      implicit integer (a-z)
      real*8 h1, h2, h3, vecout, vecin
      dimension h1(n1,n1), h2(n2,n2), h3(n3,n3)
      dimension vecout(n3,n2,n1,nvc), vecin(n3,n2,n1,nvc)
      common/io/inp, iout
      call apbc(vecout,h3,vecin,n3,n3,n2*n1*nvc)
      do 10 i=1,n1
         do 20 j=1,nvc
            call apbct(vecout(1,1,i,j),vecin(1,1,i,j),
     1                 h2,n3,n2,n2)
 20      continue
 10   continue
      do 30 i=1,nvc
         call apbct(vecout(1,1,1,i),vecin(1,1,1,i),
     1              h1,n3*n2,n1,n1)
 30   continue
      return
      end       
