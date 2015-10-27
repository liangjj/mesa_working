h64909
s 00038/00000/00000
d D 1.1 94/02/16 20:35:03 mesa 1 0
c date and time created 94/02/16 20:35:03 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
      program m6200
c***begin prologue     %M%
c***date written       930515   (yymmdd)  
c***revision date      %G%      
c
c***keywords           poisson
c***author             schneider, barry(nsf)
c***source             %W%   %G%
c***purpose            solves the generalized molecular poisson equation
c                      del**2 +k**2= rho 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- local variables ---
      integer ia,ioff,need,maxcor,intoff,wpadti
      real*8 z(1)
c
      common // ia(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m6200',0)
      intoff=wpadti(ioff)
      call pm6200(z(ioff),ia(intoff))
c
c
      stop
      end
E 1
