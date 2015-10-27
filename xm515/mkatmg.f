*deck @(#)mkatmg.f	1.1  4/25/95
      subroutine mkatmg(c,ian,grid,wts,rmax,lmax,nomega,nradial,natoms,
     $                  ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,radii,
     $                  akl,grdtyp,radshls,ptrad,dograd,gradwt,iatom)
c***begin prologue     mkatmg.f
c***date written       930518  
c***revision date      4/17/95      
c
c***keywords           
c***author             
c***source             @(#)mkatmg.f	1.1   4/25/95
c***purpose            return atomic integration grid for atom 'iatom'
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       mkatmg.f
      implicit none
c     --- input variables -----
      integer rmax,lmax,mxgrd
      integer nomega,nradial,natoms,iatom
      logical adjust,dograd
      character*(*) grdtyp
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 c(3,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grid(mxgrd,3),wts(mxgrd),
     $       gradwt(mxgrd,3,natoms)
      integer ptrad(*)
c     --- output variables ---
      integer ngrid,radshls
      real*8 radii(natoms)
c     --- scratch arrays ---
      real*8 rnuc(natoms,natoms),vwts(mxgrd)
      real*8 amu(natoms,natoms),pwtx(natoms),rr(natoms)
      real*8 akl(natoms,natoms)
c     --- local variables ---
      integer mxrpts,mxang
      parameter (mxrpts=1000,mxang=852)
      integer atom
      integer inp,iout
      real*8 pr(mxrpts),wr(mxrpts)
      real*8 pa(mxang,3),wa(mxang)
      real*8 rhomax,rhomx
      integer nodeid
      integer ierr,stderr

c
c
      common /io/inp,iout
      ierr=stderr()
      if (grdtyp.eq.'nogrid') then
         radshls=0
         ngrid=0
         return
      endif
c
c     --- generate an atomic grid 
      if (grdtyp.eq.'sg1') then
         adjust=.false.
         if (ian(iatom) .le. 18) then
            rhomx=rhomax('slater',ian(iatom))
         else 
            rhomx=rhomax('clementi',ian(iatom))
         endif
         call sg1(grid,wts,mxrpts,pr,wr,mxang,pa,wa,
     $            rhomx,c(1,iatom),mxgrd,nradial,ian(iatom),
     $            ngrid,ptrad)
         radshls=nradial
      else if(grdtyp(1:3).eq.'tag') then
         adjust=.true.
         do 20 atom=1,natoms
            radii(atom)=sqrt(rhomax('bragg',ian(atom)))
   20    continue
         call tagrid(grid,wts,ian(iatom),c(1,iatom),mxrpts,mxgrd,
     $               nradial,ngrid,ptrad,grdtyp,pr,wr,
     $               mxang,pa,wa)
         radshls=nradial
      else if(grdtyp.eq.'general') then
c        adjust=.false.
         rhomx=rhomax('bragg',ian(iatom))
         call atgrd(grid,wts,lmax,nomega,rmax,
     $              ngrid,pr,wr,mxang,pa,wa,rhomx,c(1,iatom),
     $              mxgrd,ptrad)
         radshls=rmax-1
      else
         if (nodeid().eq.0)
     $        write(iout,*) 'unrecognized grid for atom',iatom,':',
     $        grdtyp
         call plnkerr('bad grid in mkatmg',900)
      endif
c
c        --- generate the voronoi weightings --
      call vorgrad(natoms,c,grid,wts,mxgrd,
     $             ngrid,iatom,vwts,rnuc,amu,
     $             pwtx,rr,adjust,radii,akl,dograd,gradwt)
c
c
      return
      end
