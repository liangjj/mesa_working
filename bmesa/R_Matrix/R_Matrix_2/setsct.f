*deck setsct.f
c***begin prologue     setsct
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            set up scattering information.
c***                   
c***references         
c
c***routines called    
c***end prologue       setsct
      subroutine setsct(pham,energy,sym,n,nc,nen,card,dim)
      implicit integer (a-z)
      integer*8 pham
      real*8 ham, energy
      character*80 chrkey
      character*16 eunit
      character*(*) card, sym
      dimension energy(nen)
      common/io/inp, iout
      pointer (pham,ham(1))
      nen=intkey(card,'number-of-energies',1,' ')
      eunit=chrkey(card,'units','energy',' ')
      call fparr(card,'energy',energy,nen,' ')
      if(eunit.eq.'wave-vector') then
         do 10 ene=1,nen
            energy(ene)=energy(ene)*energy(ene)*.5d0
 10      continue
      endif   
      if(dim.eq.1) then
         return
      else
c
c        define the incident coordinate in order to specify
c        the zero of energy.  this is a convenience in order
c        define the incident energies so that we are in the
c        continuum.
c
         offset=0
         tcoord=intkey(card,'target-coordinate',1,' ')            
         np=nc
         if(sym.eq.'unsymmetric') then
            call iosys('read integer "number of x channels" from ham',
     1                  1,nx,0,' ')
            call iosys('read integer "number of y channels" from ham',
     1                  1,ny,0,' ')
            np=nx
            if(tcoord.eq.2) then
               offset=nx
               np=ny
            endif
         endif         
         e=offset+1
         g=e+n
         e0=g+n*nc
         call relene(energy,ham(e0),nen,np)
      endif
      return      
      end       






