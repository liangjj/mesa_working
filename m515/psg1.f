*deck @(#)psg1.f	5.1  11/28/95
      subroutine psg1( nr, rpts, atnum, atrad, nrblk, angmx, nlmx,
     $     rblkord, rblkbgn, rblksiz,nlm)
c***begin prologue     psg1.f
c***date written       931216  
c***revision date      11/28/95      
c   march 13, 1994     rlm at lanl
c      fixing bug associated with number of lebedev points for l=13
c***keywords           numerical integration, standard grid
c                      angular quadrature 
c***author             RUSSO, thomas (lanl)
c***source             @(#)psg1.f	5.1   11/28/95
c***purpose            returns the points and weights associated with
c                      gaussian quadrature for a number of radial shells
c                      each shell is a lebedev grid scaled appropriately
c                      grid is pruned according to scheme in Gill,Johnson 
c                      and Pople
c***description
c
c***references
c                     Gill, Johnson, Pople, Chem. Phys. Lett 209, 506 (1993)
c
c***end prologue      psg1.f
      implicit none
c --- input variables
      integer nr,atnum
      real*8 atrad
c --- input arrays
      real*8 rpts(nr)
c --- output variables
      integer nrblk,angmx,nlmx
c --- output arrays (5 because always only 5 blocks in sg1
      integer rblkord(5), rblkbgn(5), rblksiz(5),nlm(5)

      if (atnum.le. 2) then
c        H-He
         a1=.25d0
         a2=.5d0
         a3=1.0d0
         a4=4.5d0
      else if (atnum .le. 10) then
c        Li-Ne
         a1=.1667d0
         a2=.5d0
         a3=.9d0
         a4=3.5d0
      else if (atnum .le. 18) then
c        Na-Ar
         a1=.1d0
         a2=.4d0
         a3=.8d0
         a4=2.5d0
      else
c        K-Kr; these were not in the original sg1 paper, but are
c              present in g92/dft. i have no idea if they're ok.
        a1=0.02d0
        a2=0.1d0
        a3=0.2d0
        a4=3.5d0
c        use bragg-slater radius for Z>18
      endif

C SG-1 always does these the same:
c rblkord is 2l+1 for l=1,4,7,11,7, nlm=(l+1)**2
      nrblk=5
      angmx=23
      nlmx=144
      rblkord(1)=3
      nlm(1)=4
      rblkord(2)=9
      nlm(2)=25
      rblkord(3)=15
      nlm(3)=81
      rblkord(4)=23
      nlm(4)=144
      rblkord(5)=15
      nlm(5)=81
      call izero(rblksiz,5)
      rblkbgn(1)=1
      rblkbgn(2)=0
      rblkbgn(3)=0
      rblkbgn(4)=0
      rblkbgn(5)=0

      do 100 ir=1,nr
         if (rpts(i).lt. a1*atrad) then
            rblksiz(1)=rblksiz(1)+1
         else if (rpts(i).lt.a2*atrad) then
            if (rblkbgn(2).eq.0) rblkbgn(2)=ir
            rblksiz(2)=rblksiz(2)+1
         else if (rpts(i).lt.a3*atrad) then
            if (rblkbgn(3).eq.0) rblkbgn(3)=ir
            rblksiz(3)=rblksiz(3)+1
         else if (rpts(i).lt.a4*atrad) then
            if (rblkbgn(4).eq.0) rblkbgn(4)=ir
            rblksiz(4)=rblksiz(4)+1
         else
            if (rblkbgn(5).eq.0) rblkbgn(5)=ir
            rblksiz(5)=rblksiz(5)+1
         endif
 100  continue 
      return
      end
