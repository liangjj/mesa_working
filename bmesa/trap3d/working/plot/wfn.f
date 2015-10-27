*deck wfn.f
c***begin prologue     wfn
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             
c***purpose            print wavefunctions in format for easy plotting
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       wfn
      subroutine wfn(psi,vnl,eig,wt,arr,ind,n,nvals,maxd)
      implicit integer (a-z)
      real*8 psi, vnl, eig, wt, natoms, natmin, natmax, delat, energy 
      real*8 arr, prob 
      character*3 itoc
      character*16 fptoc
      character*80 title
      dimension psi(n), vnl(n), eig(maxd,3), wt(maxd,3)
      dimension nd(3), ind(*), arr(*)
      common/io/inp, iout
      call iosys('read integer "number of dimensions" from bec',
     1            1,ndim,0,' ')
      call iosys('read real "minimum number of atoms" from bec',
     1            1,natmin,0,' ')
      call iosys('read real "maximum number of atoms" from bec',
     1            1,natmax,0,' ')
      call iosys('read integer "number of steps" from bec',
     1            1,nstep,0,' ')
      do 10 i=1,ndim
         call iosys('read integer "number of points for '//
     1              'dimension = '//itoc(i)//'" from bec',1,
     2               nd(i),0,' ') 
         call iosys('read real "points for dimension = '//itoc(i)
     1              //'" from bec',nd(i),eig(1,i),0,' ')
         call iosys('read real "weights for dimension = '//itoc(i)
     1              //'" from bec',nd(i),wt(1,i),0,' ')
 10   continue     
      natoms=natmin
      delat=(natmax - natmin)/nstep
      do 20 num=1,nstep
         arr(num)=natoms 
         natoms = natoms + delat
 20   continue   
      do 30 num=1,nvals
         natoms=arr(ind(num))
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
         do 40 i=1,n
            prob=psi(i)/(eig(i,1)*sqrt(wt(i,1)))
            write(iout,2) eig(i,1), prob
 40      continue             

 30   continue   
      return
 1    format(//,a80)
 2    format(5x,' r = ',e15.8,2x,'psi = ',e15.8)      
      end       

