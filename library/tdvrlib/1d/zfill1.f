*deck zfill1.f
c***begin prologue     zfill1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill array
c***references         
c
c***routines called    
c***end prologue       zfill1
      subroutine zfill1(psi0,z,n)
      implicit integer (a-z)
      real*8 psi0, rez, imagz
      complex*16 z
      dimension psi0(n,2), z(n)
      common/io/inp, iout
      do 10 i=1,n
         rez=real(z(i))
         imagz=imag(z(i))
         psi0(i,1)=rez
         psi0(i,2)=imagz 
 10   continue   
      return
      end       




