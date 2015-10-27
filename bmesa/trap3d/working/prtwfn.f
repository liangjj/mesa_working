*deck prtwfn.f
c***begin prologue     prtwfn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            print wavefunctions in format for easy plotting
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       prtwfn
      subroutine prtwfn(psi,vnl,eig,n,natmin,natmax,nstep,ndim,maxd)
      implicit integer (a-z)
      real*8 psi, vnl, eig, natoms, natmin, natmax, delat, energy 
      real*8 prob 
      character*3 itoc
      character*16 fptoc
      character*80 title
      dimension psi(n), vnl(n), eig(maxd,3)
      dimension nd(3)
      common/io/inp, iout
      do 10 i=1,ndim
         call iosys('read integer "number of points for '//
     1              'dimension = '//itoc(i)//'" from bec',1,
     2               nd(i),0,' ') 
         call iosys('read real "points for dimension = '//itoc(i)
     1              //'" from bec',nd(i),eig(1,i),0,' ')
 10   continue     
      natoms=natmin
      delat=(natmax - natmin)/nstep
      do 20 num=1,nstep
         call iosys('read real "wavefunction for '//
     1              'natoms = '//fptoc(natoms)//'" from bec',
     2               n,psi,0,' ')
         call iosys('read real "energy for natoms = '
     1               //fptoc(natoms)//'" from bec',1,
     2                 energy,0,' ')
         call iosys('read real "potential for '//
     1              'natoms = '//fptoc(natoms)//'" from bec',
     2               n,vnl,0,' ')
         title='wavefunction for energy = '//fptoc(energy)
     1         //' number of atoms = '//fptoc(natoms)
         write(iout,1) title
         do 30 i=1,n
            prob=psi(i)*psi(i)
            write(iout,2) eig(i,1), prob
 30      continue             
         natoms = natoms + delat
 20   continue   
      return
 1    format(//,a80)
 2    format(5x,' r = ',e15.8,2x,'psi**2 = ',e15.8)      
      end       
