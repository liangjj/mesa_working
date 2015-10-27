*deck @(#)sg1.f	5.1  11/28/95
      subroutine sg1(xyzgrid,grdwts,mxr,tmprpts,tmprwts,
     $               mxang,tmppts,tmpwts,atrad,atcen,mxgrd,
     $               radpts,atnum,ngrid,ptrad)
c***begin prologue     sg1.f
c***date written       931216  
c***revision date      11/28/95      
c   march 13, 1994     rlm at lanl
c      fixing bug associated with number of lebedev points for l=13
c***keywords           numerical integration, standard grid
c                      angular quadrature 
c***author             RUSSO, thomas (lanl)
c***source             @(#)sg1.f	5.1   11/28/95
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
c***end prologue      sg1.f
      implicit none
c
c
      integer mxgrd,radpts,atnum,ngrid,i,j,lebord,nang,mxr,mxang
      integer ptrad(0:mxr)
      real*8 xyzgrid(mxgrd,3),grdwts(mxgrd),atcen(3)
      real*8 tmprwts(mxr),tmprpts(mxr)
      real*8 tmpwts(mxang),tmppts(mxang,3)
      real*8 atrad,a1,a2,a3,a4
      integer inp,iout
      integer angsiz
      integer stderr,ierr,nodeid
      common /io/inp,iout

      ierr=stderr()
c
c     --- for now always use 51 points for standard grids, which will
c         be cut down to 50 points by the euler scheme.
      radpts=51
c     --- determine the standard grid prunings
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
c         call lnkerr('Ack.  You called SG1 for something above krypton')
      endif
      ngrid=0
      call eulerq(tmprpts,tmprwts,radpts,atrad)
c
c     last rad point isn't really there. (zero weight)
      radpts=radpts-1
c
      ptrad(0)=1
      do 10 i=1,radpts
         if (tmprpts(i).lt. a1*atrad) then
            lebord=3
         else if (tmprpts(i).lt. a2*atrad) then
            lebord=9
         else if (tmprpts(i).lt. a3*atrad) then
            lebord=15
         else if (tmprpts(i).lt. a4*atrad) then
            lebord=23
         else 
            lebord=15
         endif
         nang=angsiz(lebord)
         call sphere(mxang,tmppts,tmpwts,lebord,nang)
         do 20 j=1,nang
            xyzgrid(j+ngrid,1)=tmppts(j,1)*tmprpts(i)+atcen(1)
            xyzgrid(j+ngrid,2)=tmppts(j,2)*tmprpts(i)+atcen(2)
            xyzgrid(j+ngrid,3)=tmppts(j,3)*tmprpts(i)+atcen(3)
            grdwts(j+ngrid)=tmpwts(j)*tmprwts(i)
 20      continue 
         ptrad(i)=ptrad(i-1)+nang
         ngrid=ngrid+nang
 10   continue 
      if (ngrid .gt. mxgrd) then
         write(iout,*) 'sg1: ngrid exceeds mxgrid; ngrid,mxgrid',
     $                  ngrid,mxgrd
         write(iout,*) 'resubmit with scf=(mxgrid=',ngrid,')'
         call plnkerr('sg1: ngrid exceeds mxgrd.',2100)
      endif
c
c
      return
      end
