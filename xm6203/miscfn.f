*deck @(#)miscfn.f	1.2  10/27/94
c***begin prologue     miscfn
c***date written       930806   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6203, link 6203, spherical harmonics
c***author             schneider, b. i.(lanl/nsf)
c***source             m6203
c***purpose            generate cos(m*phi) and sin(m*phi)
c***                   
c***references       
c
c***routines called
c***end prologue       miscfn
      subroutine miscfn(phi,cphi,sphi,nphi)
      implicit integer (a-z)
      real *8 phi, cphi, sphi
      common /io/ inp, iout
      dimension phi(nphi), cphi(nphi), sphi(nphi)
      do 10 pnt=1,nphi
         cphi(pnt)=cos(phi(pnt))
         sphi(pnt)=sin(phi(pnt))
   10 continue
      return
      end
