*deck zfill3.f
c***begin prologue     zfill3
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill array
c***references         
c
c***routines called    
c***end prologue       zfill3
      subroutine zfill3(psi0,z1,z2,z3,n,m)
      implicit integer (a-z)
      real*8 psi0, rez, imagz
      complex*16 z1, z2, z3, prod
      dimension n(3)
      dimension psi0(m,2), z1(n(1)), z2(n(2)), z3(n(3))
      common/io/inp, iout
      count=0
      do 10 i=1,n(1)
         do 20 j=1,n(2)
    1       do 30 k=1,n(3)
               prod=z1(i)*z2(j)*z3(k)
               rez=real(prod) 
               imagz=imag(prod)  
               count=count+1
               psi0(count,1) = rez
               psi0(count,2) = imagz
 30         continue
 20      continue   
 10   continue   
      return
      end       
