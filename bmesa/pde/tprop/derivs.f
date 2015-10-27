*deck derivs.f
c***begin prologue     derivs
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            evaluate the derivative of the time dependent
c***                   wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       derivs
      subroutine derivs(h01,h02,h03,vt,type,omega,efield,rpre,width,
     1                  shift,t,psi,dpsi,dim,nd,n)
      implicit integer (a-z)
      real*8 h01, h02, h03, vt
      real*8 psi, dpsi
      real*8 omega, efield, width, shift
      character*(*) rpre, type
      dimension nd(3)
      dimension h01(nd(1),nd(1)), h02(nd(2),nd(2)), h03(nd(3),nd(3))
      dimension vt(n), psi(n,2), dpsi(n,2)
      common/io/inp, iout
      data pi/3.1415926535897932384d0/
c      call rzero(vt,n)
c      if(type.eq.'t') then
c      elseif(type.eq.'cosine') then
c      elseif(type.eq.'gaussian-pulse') then
c      endif
      do 10 i=1,n
         dpsi(i,1) = psi(i,2)
         dpsi(i,2) = - psi(i,1)
 10   continue         
      return
      end       


