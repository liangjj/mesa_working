*deck tfsol.f
c***begin prologue     tfsol
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            thomas-fermi solution for bec.
c***                   
c***references         
c
c***routines called    
c***end prologue       tfbec
      subroutine tfsol(atom,psi,p,mu,x,nat,scale,n)
      implicit integer (a-z)
      real*8 psi, p, x, nat, mass, u0, omega, hbar, pi
      real*8 mu, rtf, scale
      character*(*) atom
      dimension psi(n), p(n,n), x(n)
      common/io/inp, iout
      data call/ 0 /
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
      save call, mass, u0, omega
      if(call.eq.0) then
         call locase(atom,atom)
         if(atom.eq.'cs') then
            mass=2.2d-25
            u0=2.0d-51
c           default is 10 Hz for cs
            omega=2.d0*pi*10.d0
         elseif(atom.eq.'na') then
            mass=3.8176d-26
            u0=1.0d-50
c           default is 13.846 Hz for na
            omega=2.d0*pi*13.846d0
         else
            call lnkerr('error in atom database')
         endif
         call=1
      else
c    
c        thomas-fermi mu
c
         mu=( 15.d0*u0*nat/(8.d0*pi) )**(.4d0)
         mu=mu*( .5d0*mass*omega*omega )**(.6d0)
         rtf=(2.d0*mu/(mass*omega*omega))**.5d0
         do 10 i=1,n
            psi(i)=0.d0
            if(x(i).lt.rtf) then
               psi(i)=(mu-.5d0*mass*omega*omega*x(i)*x(i))/(nat*u0)
               psi(i)=sqrt(psi(i))
            endif
 10      continue
         do 20 i=1,n
            psi(i)=psi(i)*x(i)/(scale*p(i,i))
 20      continue     
      endif
      return
      end       
