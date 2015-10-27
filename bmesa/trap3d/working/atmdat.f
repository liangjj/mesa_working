*deck atmdat.f
c***begin prologue     atmdat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            atomic date for bec.
c***                   
c***references         
c
c***routines called    
c***end prologue       atmdat
      subroutine atmdat(atom,mass,scatl,u0,omega,natmax,
     1                  sfac,type)
      implicit integer (a-z)
      real*8 mass, u0, scatl, omega, natmax, sfac, hbar, pi
      real*8 rtf, mutf, avomeg, f1, f2
      dimension omega(3), sfac(3)
      character*(*) atom, type
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
c     hbar in joule-sec      
      data hbar/1.054592d-34/
c       frequency of trap in Radians/Sec           
c       mass in kilograms of atom.           
c       scattering length in meters                              
      call locase(atom,atom)
      if(atom.eq.'cs') then
         mass=2.2d-25
         u0=2.0d-51
c           default is 10 Hz for cs
         omega(1)=2.d0*pi*10.d0
         omega(2)=2.d0*pi*10.d0
         omega(3)=2.d0*pi*10.d0
      elseif(atom.eq.'na') then
         mass=3.8176d-26
         u0=1.0d-50
         if(type.eq.'spherical') then
c        default is 13.846 Hz for na
            omega(1)=2.d0*pi*13.846d0
            omega(2)=2.d0*pi*13.846d0
            omega(3)=2.d0*pi*13.846d0
         elseif(type.eq.'cylindrical') then
            omega(1)=2.d0*pi*16.93d0
            omega(2)=13.585d0*omega(1)            
            omega(3)=omega(2)
         endif   
      else
         call lnkerr('error in atom database')
      endif
      avomeg=( omega(1)*omega(2)*omega(3) )**(1.d0/3.d0)
      scatl=mass*u0/(4.d0*hbar*hbar*pi)
      sfac(1)=sqrt(hbar/(mass*omega(1)))
      sfac(2)=sqrt(hbar/(mass*omega(2)))
      sfac(3)=sqrt(hbar/(mass*omega(3)))
c    
c     make estimate of thomas-fermi radius
c
      mutf=( 15.d0*u0/(8.d0*pi) )**(.4d0)
      mutf=mutf*( .5d0*mass*avomeg*avomeg )**(.6d0)
      rtf=( 2.d0*mutf/(mass*avomeg*avomeg) )**.5d0
      f1=natmax**.4d0
      f2=natmax**.2d0
      rtf=rtf*f2
      write(iout,1) rtf, mutf*f1/(hbar*avomeg)
      return
 1    format(/,5x,'thomas-fermi estimate for trap size  = ',e15.8,/,5x,
     1            'largest value of thomas-fermi energy = ',e15.8)
      end       
