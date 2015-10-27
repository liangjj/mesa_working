*deck @(#)mofi.f	1.1  11/30/90
      subroutine mofi(natoms,ian,c,atmass,pmom,eigvec)
c***begin prologue     mofi
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           moment of inertia, atomic masses, center of mass
c***author             gauss82
c***source             @(#)mofi.f	1.1   11/30/90
c***purpose            computes the principal moments of inertia.
c***description
c     call mofi(natoms,ian,c,atmass,pmom,eigvec)
c       pmom    vector(3) of principal moments of inertia (amu/bohr**2).
c
c***references
c***routines called    sspev(math)
c***end prologue       mofi
      implicit real*8(a-h,o-z)
      dimension c(3,*), ian(*), atmass(*), pmom(*), eigvec(*), com(3),
     $          t(6), e(9), e2(9)
      data zero/0.d+00/
c
c     compute the position of the center of mass and translate
c     it to the origin.
c
      com(1) = zero
      com(2) = zero
      com(3) = zero
      totwt = zero
      do 20 iat = 1, natoms
          wt = atmass(iat)
          totwt = totwt + wt
          com(1) = com(1) + wt*c(1,iat)
          com(2) = com(2) + wt*c(2,iat)
   20     com(3) = com(3) + wt*c(3,iat)
      com(1) = com(1) / totwt
      com(2) = com(2) / totwt
      com(3) = com(3) / totwt
c
c     compute the principal moments.
c
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
