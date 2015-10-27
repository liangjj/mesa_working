*deck mkphi.f
c***begin prologue     mkphi
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       mkphi
      subroutine mkphi(phi,u,c,lst,ns,n,type)
      implicit integer (a-z)
      real*8 phi, u, c
      character*(*) type
      dimension phi(n), u(n,n), c(ns), lst(ns)
      common/io/inp, iout
      call rzero(phi,n)
      do 10 i=1,n
         do 20 j=1,ns
            phi(i) = phi(i) + u(i,lst(j))*c(j)
 20      continue
 10   continue   
      write(iout,1) type, ( lst(i),c(i), i=1,ns)
      return
 1    format(/,1x,'list of ',a1,' states and coefficients '//
     1            'superposed',(/,1x,5('state = ',i3,1x,
     2            'coefficient = ',e15.5)))
      end       

