h65474
s 00029/00000/00000
d D 1.1 94/02/16 20:34:55 mesa 1 0
c date and time created 94/02/16 20:34:55 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
c***begin prologue     aij
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           size, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate atomic sizes based on bragg-slater radii
c***references         becke papers in jcp on lda and poisson equation.
c
c***routines called
c***end prologue       size
      function aij (za,zb)
      implicit integer (a-z)
      real*8 aij, za, zb, r, chi, uab
      dimension r(10)
      data r /10*1.d0 /
      common /io/ inp, iout
      na=za
      nb=zb
      chi=r(na)/r(nb)
      uab=(chi-1.d0)/(chi+1.d0)
      aij=uab/(uab*uab-1.d0)
      if(abs(aij).gt..5d0) then
         aij=sign(.5d0,aij)
      endif   
      return
      end

E 1
