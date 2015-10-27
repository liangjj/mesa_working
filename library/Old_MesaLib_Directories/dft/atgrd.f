*deck @(#)atgrd.f	5.2  2/5/95
      subroutine atgrd(xyzgrid,grdwts,lebord,nleb,radpts,ngrid,tmprpts,
     $                 tmprwts,mxang,tmppts,tmpwts,alpha,atcen,mxgrd,
     $                 ptrad)
c***begin prologue     atgrd.f
c***date written       930518  
c***revision date      2/5/95      
c
c***keywords           spherical harmonics, numerical integration,
c                      angular quadrature 
c***author             RUSSO, thomas (lanl)
c***source             @(#)atgrd.f	5.2   2/5/95
c***purpose            returns the points and weights associated with
c                      gaussian quadrature for a number of radial shells
c                      each shell is a lebedev grid scaled appropriately
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       atgrd.f
      implicit none
c     --- input variables -----
      integer mxgrd,mxang,radpts
      integer lebord,ngrid,nleb
c     --- input arrays (unmodified) ---
      real*8 atcen(3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ptrad(0:radpts)
      real*8 xyzgrid(mxgrd,3),grdwts(mxgrd)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 tmpwts(mxang),tmppts(mxang,3)
      real*8 tmprpts(radpts),tmprwts(radpts),alpha
c     --- local variables ---
      
      call eulerq(tmprpts,tmprwts,radpts,alpha)
c     --- last point not really used in Euler-McLaurin
c
      call sphere(mxang,tmppts,tmpwts,lebord,nleb)
      ngrid=nleb*(radpts-1)
c
      call catenat(radpts-1,mxgrd,mxang,nleb,tmprpts,tmprwts,tmppts,
     $             tmpwts,xyzgrid,grdwts,ptrad,atcen)
c 
c
      return
      end
