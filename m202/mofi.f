*deck @(#)mofi.f	5.1  11/6/94
      subroutine mofi(natoms,ian,c,atmass,pmom,eigvec)
c***begin prologue     mofi.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c***keywords           moment of inertia, atomic masses, center of mass
c***author             binkley, et al. (g82)
c***source             @(#)mofi.f	5.1   11/6/94
c***purpose            computes the principal moments of inertia.
c***description
c     call mofi(natoms,ian,c,atmass,pmom,eigvec)
c       pmom    vector(3) of principal moments of inertia (amu/bohr**2).
c
c***references
c***routines called    
c***end prologue       mofi.f
      implicit none
c     --- input variables -----
      integer natoms
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 c(3,natoms)
      real*8 atmass(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 pmom(3),eigvec(3,3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer iat,i,ierr
      real*8 zero,totwt,x,y,z,wt
      real*8 com(3),t(6),e(9),e2(9)
c
      data zero/0.d+00/
      save zero
c
c     --- compute the position of the center of mass and translate
c         it to the origin.
      com(1) = zero
      com(2) = zero
      com(3) = zero
      totwt = zero
      do 20 iat = 1, natoms
          wt = atmass(iat)
          totwt = totwt + wt
          com(1) = com(1) + wt*c(1,iat)
          com(2) = com(2) + wt*c(2,iat)
          com(3) = com(3) + wt*c(3,iat)
   20 continue
      com(1) = com(1) / totwt
      com(2) = com(2) / totwt
      com(3) = com(3) / totwt
c
c     --- compute the principal moments.
      do 60 i = 1, 6
   60     t(i) = zero
      do 80 iat = 1, natoms
          wt = atmass(iat)
          x  = c(1,iat) - com(1)
          y  = c(2,iat) - com(2)
          z  = c(3,iat) - com(3)
          t(1) = t(1) + wt * (y*y+z*z)
          t(3) = t(3) + wt * (x*x+z*z)
          t(6) = t(6) + wt * (x*x+y*y)
          t(2) = t(2) - wt * x * y
          t(4) = t(4) - wt * x * z
   80     t(5) = t(5) - wt * y * z
      call rsp(3,3,6,t,pmom,1,eigvec,e,e2,ierr)
c
c
      return
      end
