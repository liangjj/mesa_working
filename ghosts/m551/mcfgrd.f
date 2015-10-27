*deck %W%  %G%
      subroutine mcfgrd(grad,f,fab,nco,nob,nbf,c,tv,temp,deg)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy fab,f,grad,c,temp,tv
cc
      dimension fab(nbf,2),f(nbf,2),grad(nbf,2),c(nbf,2),temp(nbf,2)
      dimension tv(2)
c
c----------------------------------------------------------c
c    core contributions to the lagrangian are obtained
c    from fock operators
c        f  ... core fock operator
c        fab .. sum d(ab)*(2*j(ab)-k(ab))
c               where a,b are active orbitals
c               d(ab) one-particle density matrix
c               j(ab) , k(ab) coulomb & exchange operators
c
c  the lagrangian is in the ao basis at this point
c----------------------------------------------------------c
c
      common / number / zero,pt5,one,two,four,eight
c
      do 10 i=1,nbf
         do 20 j=1,nbf
            f(i,j) = deg * ( f(i,j) + f(i,j) ) + fab(i,j) + fab(j,i)
 20      continue
 10   continue
c
      call apbc(grad,f,c,nbf,nbf,nco)
c
c
      return
      end
