*deck rmtrx.f
c***begin prologue     rmtrx
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            calculate rmatrix and phase shifts.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       rmtrx
      subroutine rmtrx(p,q,eig,vec,energy,coul,n)
      implicit integer (a-z)
      real*8 p, q, eig, vec, energy, rmat, k, phase, fac
      real*8 sn, cn, dsn, dcn, tandel, cpse, dcpse, argm, dargm
      character*80 title
      logical coul
      dimension p(n,n), q(n), eig(n), vec(n,n)
      common/io/inp, iout
c
      write(iout,1)
      k=sqrt(2.d0*energy)
      cpse=0.d0
      dcpse=0.d0
      if(coul) then
        cpse=log(2.d0*k*q(n))/k
        dcpse=1.d0/(k*q(n))
        write(iout,2) cpse, dcpse
      else
        write(iout,3)
      endif
      argm=k*q(n)+cpse
      dargm=k+dcpse
      rmat=0.d0
      fac=p(n,n)
      do 10 i=1,n
         rmat = rmat + vec(n,i)*fac*fac*vec(n,i)/(eig(i)-energy)
   10 continue
      rmat=.5d0*rmat
      sn=sin(argm)
      cn=cos(argm)
      dsn=dargm*cn
      dcn=-dargm*sn
      tandel=( rmat*dsn - sn )/( cn - rmat*dcn )
      phase=atan(tandel)
      write(iout,4) energy, rmat, phase
      return
 1    format(/,25x,'scattering calculation',/)
 2    format(/,5x,'coulomb calculation:',/,5x,
     1            'coulomb phase               = ',e15.8,/,5x,
     2            'derivative of coulomb phase = ',e15.8) 
 3    format(/,5x,'short range potential')
 4    format(/,5x,'energy = ',e15.8,1x,'rmatrix = ',e15.8,/,5x,
     1            'phase shift = ',e15.8)
      end       
