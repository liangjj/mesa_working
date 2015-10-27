*deck vadiab.f
c***begin prologue     vadiab
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            adiabatic potential matrix
c***                   
c***description        calculate the adiabatic hyperspherical potential
c***                   in the dvr representation.
c***                   
c***references         
c
c***routines called    
c***end prologue       vadiab
      subroutine vadiab(rho,phi,x,y,vad,type,a,b,
     1                  n,dim)
      implicit integer (a-z)
      real*8 rho, phi, x, y, vad
      real*8 a, b
      character*(*) type
      character*80 title
      dimension a(dim,dim), b(dim,dim), d(dim,dim)
      dimension phi(n), x(n), y(n), vad(n)
      common/io/inp, iout
      do 10 i=1,n
         x(i)=rho*sin(phi(i))
         y(i)=rho*cos(phi(i))
 10   continue   
      if(type.eq.'none') then
         return
      elseif(type.eq.'well') then
         call vwell(vad,x,b(1,1),a(1,1),n,.false.)
         call vwell(vad,y,b(2,2),a(1,1),n,.false.)
         call vwad12(vad,x,y,b(2,1),a(2,1),n,.false.)
      elseif(type.eq.'exponential') then
         call vexp(vad,x,a(1,1),b(1,1),n,.false.)
         call vexp(vad,y,a(2,2),b(2,2),n,.false.)
         call vadpr(vad,x,y,a(2,1),b(2,1),n,.false.)
      elseif(type.eq.'2d-model') then
         call vexp(vad,x,a(1,1),b(1,1),n,.false.)
         call vexp(vad,y,a(2,2),b(2,2),n,.false.)
         call vadex12(vad,x,y,a(2,1),b(2,1),n,.false.)
      elseif(type.eq.'coulomb') then
         call vcoul(vad,x,a(1,1),n,.false.)
         call vcoul(vad,y,a(2,2),n,.false.)
         call vcad12(vad,x,y,a(1,2),n,.false.)
      else
         call lnkerr('error in potential type')
      endif
c      title='potential'
c      call prntrm(title,vad,n,1,n,1,iout)
      return
      end       





