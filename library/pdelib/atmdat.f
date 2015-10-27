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
      subroutine atmdat(atom,mass,omega)
      implicit integer (a-z)
      real*8 mass, omega, pi
      dimension omega(3)
      character*(*) atom
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
c       frequency of trap in Radians/Sec           
c       mass in kilograms of atom.           
c       scattering length in meters                              
      call locase(atom,atom)
      if(atom.eq.'cs') then
         mass=2.2d-25
c           default is 10 Hz for cs
         omega(1)=2.d0*pi*10.d0
         omega(2)=2.d0*pi*10.d0
         omega(3)=2.d0*pi*10.d0
      elseif(atom.eq.'na') then
         mass=3.8176d-26
c        default is 13.846 Hz for na
         omega(1)=2.d0*pi*13.846d0
         omega(2)=2.d0*pi*13.846d0
         omega(3)=2.d0*pi*13.846d0
      elseif(atom.eq.'na-3d') then
         mass=3.8176d-26
         omega(1)=2.d0*pi*177.0d0
         omega(2)=sqrt(2.d0)*omega(1)
         omega(3)=2.d0*omega(1)
      else
         call lnkerr('error in atom database')
      endif
      return
      end       
