*deck tocord.f
c***begin prologue     tocord
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       tocord
      subroutine tocord(rhs,crhs,p,n,mattyp)
      implicit integer (a-z)
      real*8 rhs, p
      complex*16 crhs
      character*(*) mattyp
      dimension crhs(n), rhs(n), p(n,n)
      common/io/inp, iout
      if(mattyp.eq.'complex') then
         do 10 i=1,n
            crhs(i) = crhs(i)*p(i,i)
 10      continue            
      else
         do 20 i=1,n
            rhs(i) = rhs(i)*p(i,i)
 20      continue               
      endif            
      return
      end       
